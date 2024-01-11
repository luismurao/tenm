#' Function to clean duplicated coordinates data
#' @description Clean up duplicated longitude and latitude data by year using a
#' threshold distance. The main propose of this function is eliminate data
#' occurrence points spliced in each different year or date
#' by a given distance threshold. This function also allows you to clean
#' duplicate data by pixel at a given resolution of a raster mask.
#' @param this_species Species Temporal Data object
#' see \code{\link[tenm]{sp_temporal_data}}.
#' @param threshold A numeric value representing the distance
#' between coordinates to be considered as a duplicate.
#' @param by_mask Logical. If TRUE the elimination of duplicates will be done
#' using a raster layer as a mask; If False the elimination of duplicates will
#' be done by the distance threshold.
#' @param raster_mask An object of class SpatRaster that will be used to clean
#' duplicates that are present in the same ID pixel.
#' @param n_ngbs Number of pixel neighbors. Remove duplicates depending on how
#' many pixels range you want,"0" is use in order to remove all the duplicated
#' occurrence points present in a single pixel (1x1), "1" is to eliminate
#' duplicates in a 1 pixel-long neighborhood of a 9 pixel area (3x3), "2"
#' correspond to a neighborhood of 2 adjacent pixels of a 25 pixel area (5x5)
#' and so on depending on how much area you want to cover.
#' @return A sp.temporal.modeling object that contains a temporal data.frame.
#' This table has five columns: longitude, latitude, year, layers_dates and
#' layers_path
#'
#' @details This function is build on the basis of
#' \code{\link[tenm]{clean_dup}}. See the help of the function for more examples
#' @examples
#' library(tenm)
#' data("abronia")
#' tempora_layers_dir <- system.file("extdata/bio",package = "tenm")
#' tenm_mask <- terra::rast(file.path(tempora_layers_dir,"1939/bio_01.tif"))
#' # Clean duplicates without raster mask (just by distance threshold)
#' abt <- tenm::sp_temporal_data(occs = abronia,
#'                               longitude = "decimalLongitude",
#'                               latitude = "decimalLatitude",
#'                               sp_date_var = "year",
#'                               occ_date_format="y",
#'                               layers_date_format= "y",
#'                               layers_by_date_dir = tempora_layers_dir,
#'                               layers_ext="*.tif$")
#' abtc1 <- tenm::clean_dup_by_date(abt,threshold = terra::res(tenm_mask)[1])
#' # Check number of records
#' print(nrow(abtc1$temporal_df))
#' # Clean duplicates using a raster mask
#' abtc2 <- tenm::clean_dup_by_date(this_species = abt,
#'                                 by_mask = TRUE,
#'                                 threshold = terra::res(tenm_mask)[1],
#'                                 raster_mask = tenm_mask[1],
#'                                 n_ngbs = 0)
#' # Check number of records
#' print(nrow(abtc2$temporal_df))
#'
#' @export
#'

clean_dup_by_date <- function(this_species,threshold,by_mask = FALSE,
                              raster_mask = NULL, n_ngbs = 0){
  stopifnot(inherits(this_species, "sp.temporal.modeling"))
  df_occs_date <- this_species$temporal_df
  df_occs_dateL <- split(df_occs_date,df_occs_date$layers_path,drop=T)
  clean_by_date <- seq_along(df_occs_dateL) |>
    purrr::map_df(function(x){
      dd <- tenm::clean_dup(data = df_occs_dateL[[x]],
                            longitude = this_species$lon_lat_vars[1],
                            latitude = this_species$lon_lat_vars[2],
                            threshold = threshold,
                            by_mask = by_mask,
                            raster_mask = raster_mask,
                            n_ngbs = n_ngbs)
      return(dd)
    })

  sp.temp.data.clean <- list(temporal_df = clean_by_date,
                             sp_date_var = this_species$sp_date_var,
                             lon_lat_vars =this_species$lon_lat_vars ,
                             #layers_path = layers_all[dates_ids],
                             layers_ext= this_species$layers_ext)
  class(sp.temp.data.clean) <- c("sp.temporal.modeling")


  return(sp.temp.data.clean)
}
