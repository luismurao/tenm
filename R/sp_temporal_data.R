#' Function to create a Species Temporal Data object (STD object).
#' @description
#' Creates an object of class sp.temporal.modeling that contains
#' a list with four attributes:
#'   - temporal_df: A data frame with the following columns:
#'     - Longitude: Longitude coordinates of occurrence records.
#'     - Latitude: Latitude coordinates of occurrence records.
#'     - Date: Date variable indicating when the species were observed.
#'     - Layer Dates: Format of dates for each layer of environmental data.
#'     - Layers Path: Path to the bioclimatic layer corresponding to each year.
#'   - sp_date_var: Name of the date variable column in the occurrence records.
#'   - lon_lat_vars: Names of the longitude and latitude columns.
#'   - layers_ext: Final extension format of the environmental information
#'     (e.g., ".tif").
#' @param occs A data.frame with information about occurrence records of the
#' species being modeled.
#' @param longitude Name of the variable in 'occs' containing longitude data.
#' @param latitude Name of the variable in 'occs' containing latitude data.
#' @param sp_date_var Name of the date variable.
#' @param occ_date_format Format of dates in occurrence records.
#' Options: "y", "ym", "ymd", "mdy", "my", "dmy".
#' @param layers_date_format Format of dates in raster layers. Options:
#' "y", "ym", "ymd", "mdy", "my", "dmy".
#' @param layers_by_date_dir Directory containing folders organized by date
#' with raster layers of environmental information.
#' @param layers_ext Extension or path of each raster layer archive
#' (e.g., ".tif").
#' @return Returns a  sp.temporal.modeling object (list) with the coordinates
#' of each occurrences points, the years of observation and the path to the
#' temporal layers.
#' @details
#'  The format of dates for each layer can be organized in a particular pattern,
#'  for example year/month/day ("ymd"), year/month ("ym"), just year ("y") or
#'  some other arrangement like month/year ("my"), month/year/day ("myd"),
#'  day/month/year ("dmy").
#'
#' @examples
#' library(tenm)
#' #A data.frame with occurrences points information of Abronia graminea.
#' # See help(abronia)
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

#' @export

sp_temporal_data <- function(occs,longitude,
                             latitude,sp_date_var,
                             occ_date_format="y",
                             layers_date_format= "y",
                             layers_by_date_dir,layers_ext="*.tif$"){
  classes_sp_data <- c("data.frame")
  if(class(occs) %in% classes_sp_data){

    if(methods::is(occs, "data.frame")) {
      lon_lat_vars <- c(longitude,latitude)
      if(!all(lon_lat_vars %in% names(occs)))
        stop("\n Please provide species a valid 'longitude' and 'latitude'")
    }
    if(sp_date_var %in% names(occs)){

      nearest_date_id <- function(occ_dates,layers_dates){
        r1 <-  seq_along(occ_dates) |> purrr::map_int(function(x){
          id_layer_dir <- which.min(abs(occ_dates[x] - layers_dates))
          return(id_layer_dir)
        }) #|> do.call('c', .)
        return(r1)
      }
      date_formats <- c("y","ym","ymd","my","myd","dmy")
      if(!occ_date_format %in% date_formats)
        stop("occ_date_format should have one of these values ",
             paste(date_formats,collapse = " "))
      if(occ_date_format == date_formats[1]){
        date_in_occs <- lubridate::ymd(occs[,sp_date_var], truncated = 2L)
      } else if(occ_date_format == date_formats[2]){
        date_in_occs <- lubridate::ym(occs[,sp_date_var])
      } else if(occ_date_format == date_formats[3]){
        date_in_occs <- lubridate::ymd(occs[,sp_date_var])
      } else if(occ_date_format == date_formats[4]){
        date_in_occs <- lubridate::my(occs[,sp_date_var])
      } else if(occ_date_format == date_formats[5]){
        date_in_occs <- lubridate::myd(occs[,sp_date_var])
      } else if(occ_date_format == date_formats[6]){
        date_in_occs <- lubridate::dmy(occs[,sp_date_var])
      }

      layers_names <- list.dirs(path = layers_by_date_dir,recursive = F,
                                full.names = F)

      if(!layers_date_format %in% date_formats)
        stop("layers_date_format should have one of these values ",
             paste(date_formats,collapse = " "))
      if(layers_date_format == date_formats[1]){
        layers_dates <- lubridate::ymd(layers_names, truncated = 2L)
      } else if(layers_date_format == date_formats[2]){
        layers_dates <- lubridate::ym(layers_names)
      } else if(layers_date_format == date_formats[3]){
        layers_dates <- lubridate::ymd(layers_names)
      } else if(layers_date_format == date_formats[4]){
        layers_dates <- lubridate::my(layers_names)
      } else if(layers_date_format == date_formats[5]){
        layers_dates <- lubridate::myd(layers_names)
      } else if(layers_date_format == date_formats[6]){
        layers_dates <- lubridate::dmy(layers_names)
      }

      date_in_occs_index <- which(!is.na(date_in_occs))
      date_sps <- date_in_occs[date_in_occs_index]
      # environmental layers by date

      dates_ids <- nearest_date_id(occ_dates = date_sps,
                                   layers_dates = layers_dates)

      layers_all <- list.dirs(path = layers_by_date_dir,recursive = FALSE,
                              full.names = TRUE)
      layers_all <- base::normalizePath(layers_all,winslash = "/")
      temporal_df <- data.frame(occs[date_in_occs_index,c(lon_lat_vars,sp_date_var)],
                                layer_dates = layers_dates[dates_ids],
                                layers_path=layers_all[dates_ids])

      sp_temp_data <- list(temporal_df = temporal_df,
                           sp_date_var = sp_date_var,
                           lon_lat_vars =lon_lat_vars ,
                           #layers_path = layers_all[dates_ids],
                           layers_ext = layers_ext
                           #multiband = multiband
                           )

      class(sp_temp_data) <- c("sp.temporal.modeling")

      return(sp_temp_data)
    }
    else
      stop("\n Please provide species a valid 'date_var'")
  }
  else
    stop("\n Please provide species occurrence data as a 'data.frame'")

}
