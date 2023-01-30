#' ellipsoid_selection: Performs variable selection for ellipsoid models
#'
#' @description Performs variable selection for ellipsoid models according to
#' omission rates in the environmental space.
#' @param env_train A data frame with the environmental training data.
#' @param env_test A data frame with the environmental testing data. The default
#' is NULL if given the selection process will show the p-value of a binomial
#' test.
#' @param env_vars A vector with the names of environmental variables to be
#' used in the selection process.
#' @param nvarstest A vector indicating the number of variables to fit the
#' ellipsoids during model selection. It is allowed to test models with a
#' different number of variables (i.e. nvarstest=c(3,6)).
#' @param level Proportion of points to be included in the ellipsoids. This
#' parameter is equivalent to the error (E) proposed by Peterson et al. (2008).
#' @param mve A logical value. If TRUE a minimum volume ellipsoid will be
#' computed using the function \code{\link[MASS]{cov.rob}} of the \pkg{MASS}
#' package. If False the covariance matrix of the input data will be used.
#' @param omr_criteria Omission rate criteria. Value of the omission rate
#' allowed for the selection process. Default NULL see details.
#' @param env_bg Environmental data to compute the approximated prevalence of
#' the model. The data should be a sample of the environmental layers of the
#' calibration area.
#' @param parallel The computations will be run in parallel. Default FALSE
#' @param ncores The number of cores that will be used for the parallel process.
#' By default tenm will use the total number of available cores less one.
#' @param proc Logical if TRUE a partial roc test will be run.
#' @param proc_iter Numeric. The total number of iterations for the partial ROC
#' bootstrap.
#' @param rseed Logical. Whether or not to set a random seed for partial roc
#' bootstrap. Default TRUE.
#' @param comp_each Number of models to run in each job in the parallel
#' computation. Default 100
#' @return A data.frame with 5 columns: i) "fitted_vars" the names of variables
#' that were fitted; ii) "om_rate" omission rates of the model; iii)
#' "bg_prevalence" approximated prevalence of the model see details section;
#' iv) The rank value of importance in model selection by omission rate; v)
#' The rank value by prevalence after if the value of omr_criteria is passed.
#' @details Model selection occurs in environmental space (E-space). For each
#' variable combination the omission rate (omr) in E-space is computed using the
#' function \code{\link[tenm]{inEllipsoid}}. The results will be ordered by omr
#' and if the user-specified the environmental background "env_bg" an
#' estimated prevalence will be computed and the results will be ordered also
#' by "bg_prevalence". The number of variables to construct candidate models
#' can be specified by the user in the parameter "nvarstest". Model selection
#' will be run in parallel if the user-specified more than one set of
#' combinations and the total number of models to be tested is greater than 500.
#' If given"omr_criteria" and "bg_prevalence", the results will be shown
#' pondering those models that met the "omr_criteria" by the value of
#' "bg_prevalence".
#' For more details and examples go to \code{\link[tenm]{ellipsoid_omr}} help.
#' @export
#' @import future
#' @author Luis Osorio-Olvera <luismurao@gmail.com>
#' @references Peterson, A.T. et al. (2008) Rethinking receiver operating
#' characteristic analysis applications in ecological niche modeling. Ecol.
#' Modell., 213, 63–72.


ellipsoid_selection <- function(env_train,env_test=NULL,env_vars,nvarstest,
                                level=0.95,
                                mve=TRUE,env_bg=NULL,omr_criteria,parallel=F,
                                ncores=NULL,
                                comp_each=100,proc=FALSE,
                                proc_iter=100,rseed=TRUE){
  n_vars <- length(env_vars)
  ntest <- sapply(nvarstest, function(x) choose(n_vars,x))
  nmodels <- sum(ntest)
  cat("-----------------------------------------------------------------------------------------\n")

  cat("\t\t**** Starting model selection process ****\n")
  cat("-----------------------------------------------------------------------------------------\n\n")
  for(i in 1:length(ntest)){
    cat("A total number of",ntest[i] ,"models will be created for combinations",
        "of",n_vars, "variables taken by",nvarstest[i],"\n\n")
  }
  cat("-----------------------------------------------------------------------------------------\n")
  cat("\t **A total number of",nmodels ,"models will be tested **\n\n")
  cat("-----------------------------------------------------------------------------------------\n")
  options(future.rng.onMisuse="ignore")
  if(nmodels >100 && parallel){
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
    furrr::furrr_options(globals = c("env_train",
                                     "env_test",
                                     "env_bg",
                                     "rseed","level"),
                         packages = c("tenm"))
    plan(multisession,workers=n_cores)
    options(future.globals.maxSize= 8500*1024^2)
    model_select <- new.env()
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
    cat("-----------------------------------------------------------------------------------------\n")


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
                                   "om_rate_test")],na.rm = T)
    #mean_omr <- ifelse(is.na(mean_omr),0,mean_omr)
    rfinal$mean_omr_train_test <- mean_omr
    rfinal <- rfinal[order(rfinal$mean_omr_train_tes,
                           rfinal$bg_prevalence,
                           decreasing = F),]

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
                           decreasing = F),]
  rownames(rfinal) <- NULL
  return(rfinal)
}
