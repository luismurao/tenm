
#' Predict the potential distribution of species based on environmental
#' conditions
#' @importFrom methods new
#' @docType methods
#' @param object An object of class sp.temporal.selection
#' @param model_variables A character vector specifying the variable names used
#' to build the model.
#' @param layers A SpatRaster object or a list where each element is a
#' SpatRaster.
#' @param layers_path Path to the directory containing raster layers.
#' @param layers_ext File extension of the raster layers.
#' @param mve Logical indicating whether to use the minimum volume
#' ellipsoid algorithm.
#' @param level Proportion of data to include inside the ellipsoid
#' if mve is \code{TRUE}.
#' @param output Character indicating if the model outputs "suitability" values
#'  or "mahalanobis" distances.
#' @param ... Additional parameters passed to
#' \code{\link[tenm]{ellipsoid_projection}}.
#' @return A SpatRaster object representing predicted suitability values or
#' Mahalanobis distances to niche center.
#' @details
#' This function predicts the potential distribution of a species based on
#' environmental conditions represented by raster layers. The prediction is
#' based on the model statistics and environmental variables specified in
#' 'model_variables'. If 'mve' is \code{TRUE}, the minimum volume ellipsoid algorithm
#' is used to model the niche space. The output can be either "suitability",
#' or "mahalanobis", indicating distance to the niche center.
#' Note that each SpatRaster in the 'layers' parameter should have the
#' same number of elements (layers) as 'model_variables'. The predict method
#' assumes that variables in each SpatRaster correspond to those in
#' 'model_variables'. If layers in the 'layers' parameter are given as a
#' list of objects of class SpatRaster, then the number of prediction layers
#'  will have the same number of elements in the list.
#' @rdname predict
#' @aliases predict
#' @export
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
#'                          buffer_ngbs=NULL,n_bg=50000)
#' abbg <- tenm::bg_by_date(this_species = abex,
#'                          buffer_ngbs=10,n_bg=50000)

#' future::plan("sequential")
#' varcorrs <- tenm::correlation_finder(environmental_data =
#'                                        abex$env_data[,-ncol(abex$env_data)],
#'                                      method = "spearman",
#'                                      threshold = 0.8,
#'                                      verbose = FALSE)
#' mod_sel <- tenm::tenm_selection(this_species = abbg,
#'                                 omr_criteria =0.1,
#'                                 ellipsoid_level=0.975,
#'                                 vars2fit = varcorrs$descriptors,
#'                                 nvars_to_fit=c(3,4),
#'                                 proc = TRUE,
#'                                 RandomPercent = 50,
#'                                 NoOfIteration=1000,
#'                                 parallel=TRUE,
#'                                 n_cores=2)
#' # Prediction using variables path
#' layers_70_00_dir <- system.file("extdata/bio_1970_2000",package = "tenm")
#' # The if the 'model_variables' parameter is set to NULL, the method uses
#' # the first model in the results table (mod_sel$mods_table)
#' suit_1970_2000 <- predict(mod_sel,
#'                           model_variables = NULL,
#'                           layers_path = layers_70_00_dir,
#'                           layers_ext = ".tif$")
#' # You can select the modeling variables used to project the model
#' suit_1970_2000 <- predict(mod_sel,
#'                           model_variables = c("bio_01","bio_04",
#'                                               "bio_07","bio_12"),
#'                           layers_path = layers_70_00_dir,
#'                           layers_ext = ".tif$")
#'
#' # Pass a list containing the paths of the modeling layers
#' layers_1939_2016 <- file.path(tempora_layers_dir,c("1939","2016"))
#' suit_1939_2016 <- predict(mod_sel,model_variables = NULL,
#'                           layers_path = layers_1939_2016,
#'                           layers_ext = ".tif$")
#' # Pass a list of raster layers
#' layers_1939 <- terra::rast(list.files(layers_1939_2016[1],
#'                                       pattern = ".tif$",full.names = TRUE))
#' layers_2016 <- terra::rast(list.files(layers_1939_2016[2],
#'                                       pattern = ".tif$",full.names = TRUE))
#' layers_1939 <- layers_1939[[c("bio_01","bio_04","bio_07")]]
#' layers_2016 <- layers_2016[[c("bio_01","bio_04","bio_07")]]
#' layers_list <- list(layers_1939,layers_2016)
#' suit_1939_2016 <- predict(object = mod_sel,
#'                           model_variables = c("bio_01","bio_04","bio_07"),
#'                           layers_path = NULL,
#'                           layers = layers_list,
#'                           layers_ext = ".tif$")
#'
#' }
#'

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
                       projmods <- tenm::ellipsoid_projection(envlayers = layers,
                                                             centroid = mod$centroid,
                                                             covar = mod$covariance,
                                                             level = 0.9999,
                                                             plot = TRUE,size = 2,
                                                             output = output,
                                                             ...)
                       names(projmods) <- output


                     } else if(methods::is(layers,"list")){
                       pb <- utils::txtProgressBar(min = 0,max = length(layers),style = 3)
                       checK_if_stack <- seq_along(layers) |> purrr::map(function(x){
                         if(methods::is(layers[[x]],"SpatRaster")){
                           if(terra::nlyr(layers[[x]]) != length(mod_vars)){
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
                         suitmod <- tenm::ellipsoid_projection(envlayers = layers[[x]],
                                                               centroid = mod$centroid,
                                                               covar = mod$covariance,
                                                               level = 0.9999,
                                                               plot = TRUE,
                                                               size = 2,
                                                               output = output,
                                                               ...)
                         utils::setTxtProgressBar(pb, x, title = NULL, label = NULL)
                         names(suitmod) <- output

                         return(suitmod)

                         })
                       names(projmods) <- paste0(output,"_",1:length(layers))
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
                         names(suitmod) <- paste0(output,"_",x)

                         return(suitmod)
                       })
                       if(length(projmods) == 1){
                         projmods <- projmods[[1]]
                       } else{
                         projmods <- terra::rast(projmods)
                       }
                       #names(projmods) <- layers_path
                     }
                     gc()
                     return(projmods)
                   })

