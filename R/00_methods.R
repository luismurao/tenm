
#' Predict potential distribution for \pkg{tenm} regarding to a specific time
#' period of environmental conditions or averages.
#' @importFrom methods new
#' @docType methods
#' @param object An object of class sp.temporal.selection
#' @param model_variables A vector with variable names
#' @param layers A SpatRaster object or a list where each element is a SpatRaster.
#' @param layers_path Path to layers
#' @param layers_ext Layers extension
#' @param mve If the projection will use the minimum volume ellipsoid algorithm
#' @param level Proportion of data used to fit the minimum volume ellipsoid
#' @param output A character to indicate if model uses suitability values
#' or Mahalanobis distances. Possible values are "suitability" and "mahalanobis"
#' @param ... Additional parameters passed to
#' \code{\link[tenm]{ellipsoid_projection}}
#' @details Note that the SpatRaster in layers parameter should have the same
#' number of elements (layers) than model_variables. The predict method will
#' assume that variables in each SpatRaster are the same as the ones in
#' model_variables.
#' @rdname predict
#' @aliases predict
#' @export

methods::setMethod('predict', signature(object="sp.temporal.selection"),
                   function(object, model_variables=NULL,layers = NULL,
                            layers_path=NULL,layers_ext=NULL,
                            mve = TRUE, level=0.975,output = "suitability",...){

                     out_check <- match.arg(output,
                                            choices = c("suitability",
                                                        "mahalanobis"))


                     mod_table <- object$mods_table
                     model_vars <- stringr::str_split(mod_table$fitted_vars,",")

                     if(is.null(model_variables)){
                       message(paste0("No selected variables. Using the first model in mods_table"))
                       mod_vars <- model_vars[[1]]
                     } else{
                       mod_vars <- model_variables
                       idvars <- which(!mod_vars %in% names(object$env_bg))
                       if(length(idvars)>0L){
                         stop(paste("Not valid variable names:",model_vars[idvars],
                                     "please provide valid variable names"))
                       }
                     }
                     #----------------------------------------------------------
                     # Fit ellipsoid model in E-space
                     trian_ids <- object$temporal_df$trian_test== "Train"
                     env_data <- object$temporal_df[trian_ids,mod_vars]
                     env_data <- stats::na.omit(env_data)

                     mod <- tenm::cov_center(env_data,mve = mve,
                                             level = level,vars = mod_vars)
                     #----------------------------------------------------------

                     if(methods::is(layers,"SpatRaster")){
                       projmods <- tenm::ellipsoid_projection(envlayers = terra::rast(layers),
                                                             centroid = mod$centroid,
                                                             covar = mod$covariance,
                                                             level = 0.9999,
                                                             plot = TRUE,size = 2,
                                                             output = output,
                                                             ...)
                       names(projmods) <- paste0("suitability_proj")


                     } else if(methods::is(layers,"list")){
                       pb <- utils::txtProgressBar(min = 0,max = length(layers),style = 3)
                       checK_if_stack <- seq_along(layers) |> purrr::map(function(x){
                         if(methods::is(layers[[x]][1],"SpatRaster")){
                           if(terra::nlyr(layers[[x]]) != length(model_vars)){
                             stop(paste0("The number of layers in 'layers[[",
                                         x,"]]' should be the same as 'model_variables' object" ))
                           }
                           return(TRUE)
                         } else{
                           stop(paste0("Objet 'layers[[",
                                       x,"]]' should be a 'SpatRaster'" ))
                         }
                       })

                       projmods <- seq_along(layers) |> purrr::map(function(x){
                         suitmod <- tenm::ellipsoid_projection(envlayers = terra::rast(layers[[x]]),
                                                               centroid = mod$centroid,
                                                               covar = mod$covariance,
                                                               level = 0.9999,
                                                               plot = TRUE,
                                                               size = 2,
                                                               output = output,
                                                               ...)
                         utils::setTxtProgressBar(pb, x, title = NULL, label = NULL)

                         return(suitmod)

                         })
                       names(projmods) <- paste0("suitability_proj_",1:length(layers))
                     } else{
                       pb <- utils::txtProgressBar(min = 0,max = length(layers_path),style = 3)
                       projmods <-    seq_along(layers_path) |> purrr::map(function(x){
                         lnames <- list.files(layers_path[x],pattern = layers_ext)
                         lanames <- gsub(layers_ext,"",lnames)
                         var_ids <- which(lanames %in% mod_vars)
                         layer2proj <- list.files(layers_path[x],full.names = T,
                                                  pattern = layers_ext)[var_ids]

                         slayers <- terra::rast(layer2proj)

                         suitmod <- tenm::ellipsoid_projection(envlayers = slayers,
                                                               centroid = mod$centroid,
                                                               covar = mod$covariance,
                                                               level = 0.9999,
                                                               plot = TRUE,
                                                               size = 2,
                                                               output = output,
                                                               ...)

                         utils::setTxtProgressBar(pb, x, title = NULL, label = NULL)

                         return(suitmod)
                       })
                       if(length(projmods) == 1){
                         projmods <- projmods[[1]]
                       } else{
                         projmods <- terra::rast(projmods)
                       }
                       names(projmods) <- layers_path
                     }
                     gc()
                     return(projmods)
                   })
