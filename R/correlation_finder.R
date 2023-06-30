#' Function to find out strong correlations in a correlation matrix
#' @description The function finds out which variables have strong
#' correlations according to a correlation threshold. The output
#' returns a list of variables names that can summarize the information
#' and removes the variables that are redundant.
#' @param environmental_data A matrix or a data.frame of environmental data
#' @param method A method to estimate correlation matrix. Posible options are "spearman", "pearson" or
#' "kendall".
#' @param threshold Threshold value from which it is considered that the correlation is high.
#' @param verbose Verbose output.
#' @return Returns a vector with variable names that can summarize the information.
#' @export

correlation_finder <- function(environmental_data,method="spearman",threshold,verbose=TRUE){
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
      return()
    }
    return(list(descriptors=descriptors,list_cor=list_cor))
  }
  else
    stop("cor_mat must be a matrix or a data.frame")
}
