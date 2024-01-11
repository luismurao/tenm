#' Extract environmental data by date
#' @description Function to extract environmental data by date. This function
#' generates training and testing data sets using a random partition with a
#' definite proportion.
#' @param this_species Species Temporal Data object see
#' \code{\link[tenm]{sp_temporal_data}}.
#' @param train_prop Numeric. Training proportion. This is the proportion at
#' which each record (occurrences points) will be used as training or testing.
#' A 0.7 train_prop indicates that 70 of the data will be used to train the
#' model, and the other 30 to model itself.
#' @importFrom future plan tweak sequential
#' @return An object of class sp.temporal.env that contains a temporal
#' data.frame (temporal_df) with the following columns: latitude, longitude,
#' year, layer_dates, layers_path, cell_ids_year, and environmental data.
#' @examples
#' \dontrun{
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
#' abex <- tenm::ex_by_date(this_species = abtc,
#'                          train_prop=0.7)
#' future::plan("sequential")
#' }
#' @export

ex_by_date <- function(this_species,train_prop=0.7){

  tdf <- this_species$temporal_df
  layers_path <- var_name <- layer_val <- NULL

  layer_mask <-   list.files(tdf$layers_path[1],
                             pattern = this_species$layers_ext,
                             full.names = T)[1]

  cell_ids_year <- terra::cellFromXY(terra::rast(layer_mask),
                                      tdf[,this_species$lon_lat_vars])


  tdf$cell_ids_year <- cell_ids_year

  capasDatePath <- list.files(unique(tdf$layers_path),
                             pattern = this_species$layers_ext,
                             full.names = T,recursive = T) |>
    normalizePath(winslash = "/")


  unicos <-  paste0("/",unique(base::basename(capasDatePath)),
                    collapse = "|")
  lpaths <- gsub(pattern = unicos,replacement = "",capasDatePath)
  capasByResDF <- data.frame(layers_path=lpaths,
                             capasDatePath)

  ex_time <- seq_len(nrow(capasByResDF)) |> furrr::future_map_dfr(function(x){
    time_obs <- tdf  |> dplyr::filter(layers_path ==
                                         capasByResDF$layers_path[!!x])
    env_layers <- terra::rast(capasByResDF$capasDatePath[x])
    #layer_val <- env_layers[time_obs$cell_ids_year]
    layer_val <- terra::extract(env_layers,
                                time_obs[,this_species$lon_lat_vars])
    df1 <- data.frame(time_obs[,c(1:6)],
                      layer_val = layer_val[[2]],
                      var_name = env_layers@cpp$get_sourcenames())
    return(df1)
  },.progress = TRUE,.options = furrr::furrr_options(seed = NULL))
  gc()

  years_env <- tidyr::pivot_wider(ex_time,
                                  names_from = "var_name",
                                  values_from = "layer_val")
  years_envL <- split(years_env,years_env$layer_dates)

  trian_test <- seq_along(years_envL) |> purrr::map(function(x){
    ndata <- nrow(years_envL[[x]])
    if(ndata==1) train_test <- "Train"
    if(ndata==2) train_test <- c("Train","Test")
    if(ndata==3) {
      train_test <- sample(c("Train","Train","Test"),size = 3)
    } else if(ndata >3){
      train_test <- rep("Test",ndata)
      ids_train <- sample(ndata,size = ceiling(ndata*train_prop))
      train_test[ids_train] <- "Train"

    }
    return(train_test)
  }) |> unlist()
  years_env$trian_test <- trian_test

  #train_data <- which(years_env$trian_test=="Train")

  sp.temp.data.env <- list(temporal_df = years_env,
                           sp_date_var = this_species$sp_date_var,
                           lon_lat_vars =this_species$lon_lat_vars ,
                           layers_ext= this_species$layers_ext,
                           env_data = years_env[,-(1:6)]
                           )
  class(sp.temp.data.env) <- c("sp.temporal.env")

  return(sp.temp.data.env)

}
