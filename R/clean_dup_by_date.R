#' Function to thin occurrence data
#' Cleans up duplicated longitude and latitude data by year using a specified
#' distance threshold. The distance can be specified as a geographic distance
#' or, if a raster_mask is provided, as a pixel distance.
#' @param this_species An object of class sp.temporal.modeling representing
#' species occurrence data organized by date.
#' See \code{\link[tenm]{sp_temporal_data}}.
#' @param threshold A numeric value representing the distance threshold between
#' coordinates to be considered duplicates. Units depend on whether
#' `by_mask` is \code{TRUE} or \code{FALSE}. If \code{TRUE}, the user needs
#' to specify the number of pixels that define the neighborhood of duplicates
#' (see n_ngbs parameter).
#' @param by_mask Logical. If \code{TRUE}, the thinning process will use a
#' raster layer as a mask for defining distance in pixel units.
#' @param raster_mask An object of class SpatRaster that serves as a reference
#' to thin the occurrence data. Required if `by_mask` is \code{TRUE}.
#' @param n_ngbs Number of pixels used to define the neighborhood matrix that
#' helps determine which occurrences are duplicates:
#'   - 0 removes occurrences within the same pixel, keeping one.
#'   - 1 considers duplicates all occurrences within a distance of one pixel.
#'   - n considers duplicates all occurrences within a distance of n pixels.
#' @return An object of class sp.temporal.modeling containing a temporal
#' data.frame with cleaned occurrence data, including columns for
#'  longitude, latitude, date variable, layers_dates, and layers_path.
#'
#' @details
#' This function is based on \code{\link[tenm]{clean_dup}}. It cleans up
#' duplicated occurrences based on the specified threshold. If `by_mask`
#' is \code{TRUE}, the distance is interpreted as pixel distance using the provided
#' raster_mask; otherwise, it is interpreted as geographic distance.

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
#'                                 raster_mask = tenm_mask,
#'                                 n_ngbs = 0)
#' # Check number of records
#' print(nrow(abtc2$temporal_df))
#'
#' abtc3 <- tenm::clean_dup_by_date(this_species = abt,
#'                                 by_mask = TRUE,
#'                                 threshold = terra::res(tenm_mask)[1],
#'                                 raster_mask = tenm_mask,
#'                                 n_ngbs = 2)
#' # Check number of records
#' print(nrow(abtc3$temporal_df))
#' @export
#'

clean_dup_by_date <- function(this_species,threshold,by_mask = FALSE,
                              raster_mask = NULL, n_ngbs = 0){
  stopifnot(inherits(this_species, "sp.temporal.modeling"))
  df_occs_date <- this_species$temporal_df
  df_occs_dateL <- split(df_occs_date,df_occs_date$layers_path,drop=T)
  clean_by_date <- seq_along(df_occs_dateL) |>
    furrr::future_map_dfr(function(x){
      dd <- tenm::clean_dup(data = df_occs_dateL[[x]],
                            longitude = this_species$lon_lat_vars[1],
                            latitude = this_species$lon_lat_vars[2],
                            threshold = threshold,
                            by_mask = by_mask,
                            raster_mask = raster_mask,
                            n_ngbs = n_ngbs)
      return(dd)
    },.progress = TRUE,
    .options = furrr::furrr_options(globals = c("this_species",
                                                "df_occs_date",
                                                "df_occs_dateL"),
                                    seed = NULL))

  sp.temp.data.clean <- list(temporal_df = clean_by_date,
                             sp_date_var = this_species$sp_date_var,
                             lon_lat_vars =this_species$lon_lat_vars ,
                             #layers_path = layers_all[dates_ids],
                             layers_ext= this_species$layers_ext)
  class(sp.temp.data.clean) <- c("sp.temporal.modeling")


  return(sp.temp.data.clean)
}
