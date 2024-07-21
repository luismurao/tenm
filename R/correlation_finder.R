#' Function to find strong correlations within environmental predictors
#' @description
#' This function identifies variables with strong correlations based on a
#' specified threshold.
#' @param environmental_data A matrix or data.frame containing
#' environmental data.
#' @param method Method used to estimate the correlation matrix. Possible
#' options include "spearman" (Spearman's rank correlation),
#' "pearson" (Pearson's correlation),
#' or "kendall" (Kendall's tau correlation).
#' @param threshold Correlation threshold value. Variables with absolute
#' correlation values greater than or equal to this threshold are considered
#' strongly correlated.
#' @param verbose Logical. If \code{TRUE}, prints verbose output detailing
#' correlations.
#' @return A list with two elements:
#'   - `not_correlated_vars`: A vector containing names of variables that are
#'      not strongly correlated.
#'   - `correlation_values`: A list with correlation values for all pairs of
#'      variables.

#' @export
#' @examples
#' \donttest{
#'
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
#' abex <- tenm::ex_by_date(abtc,train_prop=0.7)
#' future::plan("sequential")
#' envdata <- abex$env_data[,-ncol(abex$env_data)]
#' ecors <- tenm::correlation_finder(environmental_data =envdata,
#'                                   method="spearman",
#'                                   threshold = 0.7 )
#' }
#'

correlation_finder <- function(environmental_data,method="spearman",threshold,
                               verbose=TRUE){
  if(is.matrix(environmental_data) || is.data.frame(environmental_data)){
    environmental_data <- stats::na.omit(environmental_data)
    cor_mat <- data.frame(stats::cor(environmental_data,method = method))
    variables_cor <- function(x,threshold){
      x <- as.numeric(x)
      vars_pos <- which(x > threshold)
      vars_neg <- which(x < (-1)*threshold)
      vars <- c(vars_pos,vars_neg)
      cors <- x[vars]
      names(cors) <- names(cor_mat)[vars]
      return(cors)
    }

    list_cor <- mapply(variables_cor,x=cor_mat,
                       threshold=threshold,SIMPLIFY = FALSE)

    nomvar <- names(cor_mat)
    list_cor2 <- sapply(nomvar, function(x)
      list_cor[[x]][which(list_cor[[x]]!=1)])


    nomvar2 <- nomvar
    descriptors <- NULL
    for(i in 1:length(list_cor2)){
      descriptors <- c(descriptors,names(list_cor2[i]))
      if(names(list_cor2[i]) %in% nomvar2){
        trash <- sort(unlist(sapply(c(descriptors,
                                      names(list_cor2[[i]])),
                                    function(x) which(x==nomvar2))))
        nomvar2 <- nomvar2[-trash]
      }
      else{
        descriptors <- descriptors[-length(descriptors)]
      }
    }


    if(verbose){
      f1 <- '********************************'
      f1 <- paste0(f1,'*********************************\n\n')
      f2 <- '---------------------------------'
      f2 <- paste0(f2,'-------------------------------\n\n')
      cat(f1)
      cat(' Here is a list of variables that can summarize your niche\n')
      cat(' information, according to the threshold of',threshold,":\n\n")
      cat(' ',descriptors,'\n\n')
      cat(f1)
      cat(f2)
      cat('Correlation list:\n\n')

      for(i in 1:dim(cor_mat)[1]){
        cat("Variable",names(list_cor)[i],"is strongly correlated with:\n\n")
        print(list_cor[[i]])
        cat('----------------------------------------------------------------\n\n')
      }
    }
    return(list(descriptors=descriptors,list_cor=list_cor))
  }
  else
    stop("cor_mat must be a matrix or a data.frame")
}

