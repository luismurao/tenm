#' ellipsoid_omr
#'
#' @description Compute the omission rate of ellipsoid models
#' @param env_data A data frame with the environmental data.
#' @param env_test A data frame with the environmental testing data. The default
#' is NULL if given the selection process will show the p-value of a binomial
#' test.
#' @param env_bg Environmental data to compute the approximated prevalence of
#' the model. The data should be a sample of the environmental layers of the
#' calibration area.
#' @param cf_level Proportion of points to be included in the ellipsoids. This
#' parameter is equivalent to the error (E) proposed by Peterson et al. (2008).
#' @param mve A logical value. If TRUE a minimum volume ellipsoid will be
#' computed using
#' the function \code{\link[MASS]{cov.rob}} of the \pkg{MASS} package. If False
#' the covariance matrix of the input data will be used.
#' @param proc Logical if TRUE a partial roc test will be run.
#' @param proc_iter Numeric. The total number of iterations for the partial ROC
#' bootstrap.
#' @param rseed Logical. Whether or not to set a random seed for partial roc
#' bootstrap. Default TRUE.
#' @return A data.frame with 5 columns: i) "fitted_vars" the names of variables
#' that were fitted; ii) "om_rate" omission rates of the model; iii)
#' "bg_prevalence" approximated prevalence of the model see details section.
#' @export

ellipsoid_omr <- function(env_data,env_test=NULL,env_bg,cf_level,mve=TRUE,
                          proc=FALSE,proc_iter=100,rseed=TRUE){
  emd <- try(tenm::cov_center(data = env_data,
                               mve = mve,
                               level = cf_level,
                               vars = 1:ncol(env_data)),
             silent = TRUE)

  message1 <- attr(emd,"class")== "try-error"
  if(length(message1)>0L)
    return()

  in_e <-  inEllipsoid(centroid = emd$centroid,
                       eShape = emd$covariance,
                       env_data = env_data,
                       level = cf_level)

  fails_train_ids <- which(in_e$in_Ellipsoid== 0)

  if(length(fails_train_ids)>0){
    fails_train_ids <- paste0(fails_train_ids,collapse = ",")
  } else {
    fails_train_ids <- NA
  }

  occs_table <- table( in_e$in_Ellipsoid)

  succsID <- which(names(occs_table) %in% "1")
  failsID <- which(names(occs_table) %in% "0")

  occs_succs <-  if(length(succsID)>0L){
    occs_table[[succsID]]
  } else{
    0
  }
  occs_fail <-  if(length(failsID)>0L){
    occs_table[[failsID]]
  } else{
    0
  }

  a_train <-  occs_fail
  omrate_train <- a_train /nrow( in_e)

  d_results <- data.frame(fitted_vars =paste(names(emd$centroid),
                                             collapse =  ","),
                          nvars=length(emd$centroid),
                          om_rate_train= omrate_train,
                          non_pred_train_ids = fails_train_ids)
  if(is.data.frame(env_test) || is.matrix(env_test)){
    in_etest <-  tenm::inEllipsoid(centroid = emd$centroid,
                                    eShape = emd$covariance,
                                    env_data = env_test,
                                    level = cf_level)

    fails_test_ids <- which(in_etest$in_Ellipsoid== 0)

    if(length(fails_train_ids)>0){
      fails_test_ids <- paste0(fails_test_ids,collapse = ",")
    } else {
      fails_test_ids <- NA
    }

    suits_val <- exp(-0.5*( in_etest$mh_dist))

    occs_table_test <- table(in_etest$in_Ellipsoid)

    succsID <- which(names(occs_table_test) %in% "1")
    failsID <- which(names(occs_table_test) %in% "0")

    occs_succs_test <-  if(length(succsID)>0L){
      occs_table_test[[succsID]]
    } else{
      0
    }
    occs_fail_test <-  if(length(failsID)>0L){
      occs_table_test[[failsID]]
    } else{
      0
    }
    a_test <-  occs_fail_test
    omrate_test <- a_test /nrow( in_etest)
    d_results <- data.frame(d_results,
                            om_rate_test=omrate_test,
                            non_pred_test_ids=fails_test_ids)
  }

  if(!is.null(env_bg)){

    env_bg <- data.frame(env_bg)
    in_ebg <-  tenm::inEllipsoid(centroid = emd$centroid,
                                  eShape = emd$covariance,
                                  env_data = env_bg,
                                  level = cf_level)
    suits_bg <- exp(-0.5*in_ebg$mh_dist)

    bg_table <- table(c(in_ebg$in_Ellipsoid,in_e$in_Ellipsoid))
    succs_bg_ID <- which(names(bg_table) %in% "1")
    fails_bg_ID <- which(names(bg_table) %in% "0")

    bg_succs <-  if(length(succs_bg_ID)>0L){
      bg_table[[succs_bg_ID]]
    } else{
      0
    }

    bg_fails <-  if(length(fails_bg_ID)>0L){
      bg_table[[fails_bg_ID]]
    } else{
      0
    }
    prevBG <- bg_succs/(bg_fails+bg_succs)
    d_results <-data.frame( d_results,
                            bg_prevalence= prevBG)

    if(exists("in_etest")){
      #bin_table <- table(c(in_ebg$in_Ellipsoid,
      #                     in_etest$in_Ellipsoid))
      #binBG <- bin_table[[2]]/(bin_table[[1]]+bin_table[[2]])
      test_fail <-  occs_fail_test
      test_succs <- occs_succs_test
      p_bin <- 1 - stats::pbinom(test_succs,
                                 size=test_succs+test_fail,
                                 prob = prevBG)
      d_results <-data.frame( d_results,
                              pval_bin=p_bin)
      if(proc){
        proc <- try(tenm::pROC(suits_bg,test_data = suits_val,
                            n_iter = proc_iter,rseed = rseed))
        message1 <- attr(proc,"class")== "try-error"
        if(length(message1)>0L){
          pval_proc <- NA
          mean_aucratio <- NA
          mean_auc <- NA
        } else{
          pval_proc <- proc$pROC_summary[3]
          mean_aucratio <- proc$pROC_summary[2]
          mean_auc <- proc$pROC_summary[1]
        }

        d_results <-data.frame( d_results,
                                pval_proc,
                                env_bg_paucratio= mean_aucratio,
                                env_bg_auc = mean_auc)
      }


    }

  }
  return(d_results)
}