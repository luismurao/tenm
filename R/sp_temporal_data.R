#' Function to create a Species Temporal Data object (STD object).
#' @description This function creates an object of class sp.temporal.modeling
#'  that contains a list with four attributes. Inside the list, there is a data
#'  frame called temporal_df containing five columns: the first two columns are
#'  coordinates of longitude and latitude that came from the occurrences
#'  records, the third is the date variable at which the species were
#'  registered, the fourth is the format of dates for each layer, organized in
#'  a particular pattern, for example year/month/day, year/month, just year or
#'  some other arrangement like month/year, month/year/day, day/month/year and
#'  the fifth one is the path where is stored the bioclimatic layer
#'  corresponding to each year. The other three variables present in the list
#'  are character class objects: sp_date_var [1] is the name of the date
#'  variable column available in the occurrences records, lon_lat_vars [2] are
#'  both of the columns with the coordinates of longitude and latitude
#'  correspondingly, layers_ext [1] is the final extension format of the
#'  environmental information (“.tif$”).
#' @param occs A DataFrame or a SpatialPointsDataFrame with information about
#' the occurrence records of the specie that is being modeled. It is fundamental
#' to count with exact geographical coordinates of longitude and latitude where
#' the specie was detected or at least the nearest, also a temporal column
#'   indicating the time at which was the record.
#' @param longitude If occs is a data.frame the user must indicate the variable
#' name of longitude data.
#' @param latitude If occs is a data.frame the user must indicate the variable
#' name of latitude data.
#' @param sp_date_var A date variable indicating the date of each observation.
#' The name of the variable where is stored the date of each observation.
#' @param occ_date_format Occurrences date format. It is the format with which
#' the dates of the occurrence points are organized. The possible options are
#' "y" for years; "ym" for years and months; "ymd" for year, month and day;
#' "mdy" for month, day and year; "my" for month and year; "dmy" for day, month
#' and year.
#' @param layers_date_format Raster layers of environmental information data
#' format. The possible options are "y" for years; "ym" for years and months;
#' "ymd" for year, month and day; "mdy" for month, day and year; "my" for month
#' and year; "dmy" for day, month and year.
#' @param layers_by_date_dir A directory which has contain inside other folders
#'  organized by date with the raster layers of environmental information.
#' @param layers_ext This is the extension or path of each raster layer archive.
#'  In other words, this is the object where is stored each one of the paths
#'  that leads to the location of the environmental raster layers into the inner
#'  memory.
#' @return Returns a  sp.temporal.modeling object (list) with the coordinates
#' of each occurrences points, the years of observation and the path to the
#' temporal layers.
#' @importFrom rgdal readOGR
#' @importFrom raster raster stack
#' @importFrom magrittr %>%
#'
#' @examples
#' library(tenm)
#' #A data.frame with occurrences points information of Abronia graminea.
#' See help(abronia)
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
  classes_sp_data <- c("SpatialPointsDataFrame","data.frame")
  if(class(occs) %in% classes_sp_data){
    if(methods::is(occs, "SpatialPointsDataFrame")){
      lon_lat_vars <- colnames(occs@coords)
      occs <- data.frame(occs@coords,occs@data)
    }

    if(methods::is(occs, "data.frame")) {
      lon_lat_vars <- c(longitude,latitude)
      if(!all(lon_lat_vars %in% names(occs)))
        stop("\n Please provide species a valid 'longitude' and 'latitude'")
    }
    if(sp_date_var %in% names(occs)){

      nearest_date_id <- function(occ_dates,layers_dates){
        r1 <-  seq_along(occ_dates) %>% purrr::map_int(function(x){
          id_layer_dir <- which.min(abs(occ_dates[x] - layers_dates))
          return(id_layer_dir)
        }) #%>% do.call('c', .)
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

      layers_all <- list.dirs(path = layers_by_date_dir,recursive = F,
                              full.names = T)
      temporal_df <- data.frame(occs[date_in_occs_index,c(lon_lat_vars,sp_date_var)],
                                layer_dates = layers_dates[dates_ids],
                                layers_path=layers_all[dates_ids])

      sp_temp_data <- list(temporal_df = temporal_df,
                           sp_date_var = sp_date_var,
                           lon_lat_vars =lon_lat_vars ,
                           #layers_path = layers_all[dates_ids],
                           layers_ext= layers_ext)

      class(sp_temp_data) <- c("sp.temporal.modeling")

      return(sp_temp_data)
    }
    else
      stop("\n Please provide species a valid 'date_var'")
  }
  else
    stop("\n Please provide species occurrence data as a 'data.frame' or as a
         'SpatialPointsDataFrame'")

}
