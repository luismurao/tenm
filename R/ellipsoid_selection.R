#' ellipsoid_selection: Performs models selection for ellipsoid models
#'
#' @description
#' The function performs model selection for ellipsoid models
#' using three criteria: a) the omission rate, b) the significance of partial
#' ROC and binomial tests and c) the AUC value.
#'
#' @param env_train A data frame with the environmental training data.
#' @param env_test A data frame with the environmental testing data.
#' Default is NULL.
#' @param env_vars A vector with the names of environmental variables used in
#' the selection process. To help choosing which variables to use see
#' \code{\link[tenm]{correlation_finder}}.
#' @param nvarstest A vector indicating the number of variables to fit the
#' ellipsoids during model selection.
#' @param level Proportion of points to be included in the ellipsoids,
#' equivalent to the error (E) proposed by Peterson et al. (2008).
#' @param mve Logical. If \code{TRUE}, a minimum volume ellipsoid will be computed.
#' using \code{\link[MASS]{cov.rob}} from \pkg{MASS}. If \code{FALSE}, the covariance
#' matrix of the input data will be used.
#' @param omr_criteria Omission rate criteria: the allowable omission rate for
#' the selection process. Default is NULL (see details).
#' @param env_bg Environmental data to compute the approximated prevalence
#' of the model, should be a sample of the environmental layers of
#' the calibration area.
#' @param parallel Logical. If \code{TRUE}, computations will run in parallel.
#' Default is \code{F}.
#' @param ncores Number of cores to use for parallel processing. Default uses
#' all available cores minus one.
#' @param proc Logical. If \code{TRUE}, a partial ROC test will be run.
#' @param proc_iter Numeric. Total iterations for the partial ROC bootstrap.
#' @param rseed Logical. If \code{TRUE}, set a random seed for partial ROC bootstrap.
#' Default is \code{TRUE}.
#' @param comp_each Number of models to run in each job in parallel computation.
#'  Default is 100.
#' @return A data.frame with the following columns:
#'   - "fitted_vars": Names of variables that were fitted.
#'   - "nvars": Number of fitted variables
#'   - "om_rate_train": Omission rate of the training data.
#'   - "non_pred_train_ids": Row IDs of non-predicted training data.
#'   - "om_rate_test"': Omission rate of the testing data.
#'   - "non_pred_test_ids": Row IDs of non-predicted testing data.
#'   - "bg_prevalence": Approximated prevalence of the model (see details).
#'   - "pval_bin": p-value of the binomial test.
#'   - "pval_proc": p-value of the partial ROC test.
#'   - "env_bg_paucratio": Environmental partial AUC ratio value.
#'   - "env_bg_auc": Environmental AUC value.
#'   - "mean_omr_train_test": Mean value of omission rates (train and test).
#'   - "rank_by_omr_train_test": Rank value of importance in model selection
#'     by omission rate.
#'   - "rank_omr_aucratio": Rank value by AUC ratio.
#' @details
#' Model selection occurs in environmental space (E-space). For each variable
#' combination specified in nvarstest, the omission rate (omr) in E-space is
#' computed using \code{\link[tenm]{inEllipsoid}} function.
#' Results are ordered by omr of the testing data. If env_bg is provided,
#' an estimated prevalence is computed and results are additionally ordered
#' by partial AUC. Model selection can be run in parallel.
#' For more details and examples go to \code{\link[tenm]{ellipsoid_omr}} help.
#' @export
#' @import future
#' @author Luis Osorio-Olvera <luismurao@gmail.com>
#' @references Peterson, A.T. et al. (2008) Rethinking receiver operating
#' characteristic analysis applications in ecological niche modeling. Ecol.
#' Modell. 213, 63–72. \doi{10.1016/j.ecolmodel.2007.11.008}
#'
#' @examples
#' \donttest{
#' library(tenm)
#' data("abronia")
#' tempora_layers_dir <- system.file("extdata/bio",package = "tenm")
#' abt <- tenm::sp_temporal_data(occs = abronia,
#'                               longitude = "decimalLongitude",
#'                               latitude = "decimalLatitude",
#'                               sp_date_var = "year",
#'                               occ_date_format="y",
#'                               layers_date_format= "y",
#'                               layers_by_date_dir = tempora_layers_dir,
#'                               layers_ext="*.tif$")
#' abtc <- tenm::clean_dup_by_date(abt,threshold = 10/60)
#' future::plan("multisession",workers=2)
#' abex <- tenm::ex_by_date(this_species = abtc,train_prop=0.7)
#' abbg <- tenm::bg_by_date(this_species = abex,
#'                          buffer_ngbs=10,n_bg=50000)
#' future::plan("sequential")
#' varcorrs <- tenm::correlation_finder(environmental_data =
#'                                      abex$env_data[,-ncol(abex$env_data)],
#'                                      method = "spearman",
#'                                      threshold = 0.8,
#'                                      verbose = FALSE)
#' edata <- abex$env_data
#' etrain <- edata[edata$trian_test=="Train",] |> data.frame()
#' etest <- edata[edata$trian_test=="Test",] |> data.frame()
#' bg <- abbg$env_bg
#' res1 <- tenm::ellipsoid_selection(env_train = etrain,
#'                                   env_test = etest,
#'                                   env_vars = varcorrs$descriptors,
#'                                   nvarstest = 3,
#'                                   level = 0.975,
#'                                   mve = TRUE,
#'                                   env_bg = bg,
#'                                   omr_criteria = 0.1,
#'                                   parallel = FALSE,proc = TRUE)
#' head(res1)

#' }
#'

ellipsoid_selection <- function(env_train,env_test=NULL,env_vars,nvarstest,
                                level=0.95,
                                mve=TRUE,env_bg=NULL,omr_criteria,
                                parallel=FALSE,
                                ncores=NULL,
                                comp_each=100,proc=FALSE,
                                proc_iter=100,rseed=TRUE){
  n_vars <- length(env_vars)
  ntest <- sapply(nvarstest, function(x) choose(n_vars,x))
  nmodels <- sum(ntest)
  cat("-------------------------------------------------------------------\n")

  cat("\t\t**** Starting model selection process ****\n")
  cat("-------------------------------------------------------------------\n\n")
  for(i in 1:length(ntest)){
    cat("A total number of",ntest[i] ,"models will be created for combinations",
        "of",n_vars, "variables taken by",nvarstest[i],"\n\n")
  }
  cat("-------------------------------------------------------------------\n")
  cat("\t **A total number of",nmodels ,"models will be tested **\n\n")
  cat("-------------------------------------------------------------------\n")

  if(nmodels >100 && parallel){
    options(future.rng.onMisuse="ignore")
    max_var <- max(nvarstest)
    cvars <- lapply(nvarstest, function(x) {

      cb <- utils::combn(env_vars,x)
      if(x < max_var){
        nrowNA <-max_var-nrow(cb)
        na_mat <- matrix(nrow = nrowNA,ncol=ncol(cb))
        cb <- rbind(cb,na_mat)
      }
      return(cb)
    })
    big_vars <- do.call(cbind,cvars)

    n_cores <- future::availableCores() -1
    if(ncores>n_cores || is.null(ncores)){
      n_cores <- n_cores
    } else{
      n_cores <- ncores
    }
    niter_big <- floor(nmodels/n_cores)
    if(niter_big>comp_each)
      niter_big <- comp_each
    steps <- seq(1, nmodels, niter_big)
    nsteps <- length(steps)
    if(steps[nsteps]<nmodels){
      kkk <- c(steps,  nmodels + 1)
    } else {
      kkk <- steps
      kkk[nsteps] <- kkk[nsteps] + 1
    }

    long_k <- length(kkk)
    pasos <- 1:(length(kkk) - 1)
    pasosChar <- paste0(pasos)
    globs <- c("env_train",
               "env_test",
               "env_bg")
    #furrr::furrr_options(globals = c("env_train",
    #                                 "env_test",
    #                                 "env_bg",
    #                                 "rseed","level"),
    #                     packages = c("tenm"))
    oplan <- plan(tweak(multisession, workers = n_cores))
    #plan(multisession,workers=n_cores)
    options(future.globals.maxSize= 8500*1024^2)
    model_select <- new.env()
    on.exit(plan(oplan), add = TRUE)
    for (paso in pasosChar) {
      x <- as.numeric(paso)
      #fname <- file.path(dir1,paste0("eselection_",x,".txt"))
      #if(x>n_cores) core <- 1

      cat("Doing calibration from model ",
          kkk[x],"to ",kkk[x + 1] - 1,
          "in process ",x,"\n\n")
      model_select[[paso]] %<-% {

        seq_model <- kkk[x]:(kkk[x + 1] - 1)
        combs_v <- as.matrix(big_vars[,seq_model])

        results_L <- lapply(1:ncol(combs_v),function(x_comb) {
          var_comb <- stats::na.omit(combs_v[,x_comb])
          env_data0 <- stats::na.omit(env_train[,var_comb])
          env_test0 <- stats::na.omit(env_test[,var_comb])
          env_bg0 <-   stats::na.omit(env_bg[,var_comb])
          r1 <- tenm::ellipsoid_omr(env_data = env_data0,
                                    env_test = env_test0,
                                    env_bg = env_bg0,
                                    cf_level = level,
                                    proc = proc,
                                    proc_iter,
                                    rseed=rseed)

          return(r1)
        })
        results_df <- do.call("rbind.data.frame",results_L)
        cat("Finishing calibration of models ",kkk[x],"to ",kkk[x + 1] - 1,
            "\n\n")
        return(results_df)
      }

    }
    mres <- as.list(model_select)

    cat("Finishing...\n\n")
    cat("-------------------------------------------------------------------\n")


    rfinal <- do.call("rbind.data.frame", mres )

    future::plan(sequential)
  }
  else{
    cvars <- lapply(nvarstest, function(x) utils::combn(env_vars,x))

    results_L <- lapply(1:length(cvars), function(x) {
      combs_v <- cvars[[x]]
      results_L <- lapply(1:ncol(combs_v),function(x_comb) {
        var_comb <- stats::na.omit(combs_v[,x_comb])
        env_data <- stats::na.omit(env_train[,var_comb])
        env_test <- stats::na.omit(env_test[,var_comb])
        env_bg <-   stats::na.omit(env_bg[,var_comb])
        r1 <- tenm::ellipsoid_omr(env_data = env_data,
                                   env_test = env_test,
                                   env_bg = env_bg,
                                   cf_level = level,
                                   proc = proc,
                                   proc_iter,rseed=rseed)
        return(r1)
      })
      results_df <- do.call("rbind.data.frame",results_L)
      return(results_df)
    })
    rfinal <- do.call("rbind.data.frame",results_L)
  }
  bg_omr <- c("bg_prevalence","om_rate_test") %in% names(rfinal)
  bg_omr_in <- all(bg_omr)
  if( bg_omr_in){
    rfinal[["om_rate_test"]] <- ifelse(is.na(rfinal[["om_rate_test"]]),0,
                                       rfinal[["om_rate_test"]])
    rfinal[["om_rate_train"]] <- ifelse(is.na(rfinal[["om_rate_train"]]),0,
                                       rfinal[["om_rate_train"]])
    mean_omr <- rowMeans(rfinal[,c("om_rate_train",
                                   "om_rate_test")],na.rm = TRUE)
    #mean_omr <- ifelse(is.na(mean_omr),0,mean_omr)
    rfinal$mean_omr_train_test <- mean_omr
    rfinal <- rfinal[order(rfinal$mean_omr_train_tes,
                           rfinal$bg_prevalence,
                           decreasing = FALSE),]

    rfinal <- data.frame(rfinal,rank_by_omr_train_test=1:nrow(rfinal))
    met_criteriaID_train <- which(rfinal$om_rate_train <= omr_criteria)
    met_criteriaID_test <- which(rfinal$om_rate_test <= omr_criteria)
    met_criteriaID_both <- intersect(met_criteriaID_train,
                                     met_criteriaID_test)

    if(length(met_criteriaID_train) > 0L){
      cat("\t",length(met_criteriaID_train),
          "models passed omr_criteria for train data\n")
    }
    if(length(met_criteriaID_test) > 0L){
      cat("\t",length(met_criteriaID_test),
          "models passed omr_criteria for test data\n")

    }
    if(length(met_criteriaID_both) > 0L){
      cat("\t",length(met_criteriaID_both),
          "models passed omr_criteria for train and test data\n")
    }
    else{
      cat("\tNo model passed the omission criteria ranking by mean omission rates\n")
      return(rfinal)
    }
    best_r <- rfinal[met_criteriaID_both,]
    if(proc){
      best_r <- best_r[order(best_r$env_bg_paucratio,
                             decreasing = TRUE),]
    }

    rfinal <- rbind(best_r,
                    rfinal[-met_criteriaID_both,])
    if(proc){
      rfinal <- data.frame(rfinal,
                           rank_omr_aucratio=1:nrow(rfinal))
    }
  }
  else
    rfinal <- rfinal[order(rfinal$om_rate_train,
                           decreasing = FALSE),]
  rownames(rfinal) <- NULL
  return(rfinal)
}
