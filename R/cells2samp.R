#' Helper function to generate cell IDs to be used for choosing  the
#' environmental background data
#' @description Returns pixel IDs to be sample for generating
#' environmental background form modeling layers.
#' @param data A data.frame with longitude and latitude data
#' @param longitude A character vector of the column name of longitude.
#' @param latitude A character vector of the column name of latitude.
#' @param cell_ids A numeric vector with pixel ids. The default values NULL.
#' Use this parameter if you have obtained this information using the function
#' cellFromXY.
#' @param buffer_ngbs Number of pixel neighbors around occurrences to be used
#' to build the buffer.
#' @param raster_mask An object of class RasterLayer that will be used to
#' obtain pixel IDs.
#' @param process_ngbs_by Numeric. Estimates neighbor cells each x cells. This
#' is for memory management.
#' @param n_bg Number of background pixels.
#' @param progress Logical. Show computation progress
#' @examples
#' # cells to sample
#' data(abronia)
#' temporal_layer <- system.file("extdata/bio/2016/bio_01.tif",package = "tenm")
#' raster_mask <- terra::rast(temporal_layer)
#' set.seed(123)
#' samp_01 <- tenm::cells2samp(data = abronia,
#'                             longitude = "decimalLongitude",
#'                             latitude = "decimalLatitude",
#'                             cell_ids = NULL,
#'                             buffer_ngbs = 4,
#'                             raster_mask = raster_mask,
#'                             process_ngbs_by = 10,
#'                             n_bg = 50000,
#'                             progress =TRUE)
#'
#' @export
cells2samp <- function(data,longitude,latitude,cell_ids = NULL,buffer_ngbs = 2,
                       raster_mask,process_ngbs_by = 10,n_bg = 50000,
                       progress =TRUE){
  if(is.null(cell_ids)){
    data <- data.frame(data[,c(longitude,latitude)])
    cell_ids <- terra::cellFromXY(raster_mask,data[,c(longitude,latitude)])
  } else if(!is.numeric(cell_ids)){
    stop("Provide valid cell numbers")
  }
  cell_ids <- base::sort(cell_ids)
  n_cells <- length(cell_ids)
  cut_offs <- ceiling(n_cells/process_ngbs_by)
  if(cut_offs > 1){
    pro_cells <- base::cut(x = seq_len(n_cells),cut_offs)
    cell_idsL <- base::split(cell_ids,pro_cells)

  } else{
    cell_idsL <- list(cell_ids)

  }

  nbase <- 2 * buffer_ngbs + 1
  ngMat <- base::matrix(rep(1, nbase * nbase), ncol = nbase, byrow = TRUE)
  ngMat[buffer_ngbs + 1, buffer_ngbs + 1] <- 0
  rcellsL <- seq_along(cell_idsL) |> purrr::map(function(x){
    adj_cells <- terra::adjacent(x = raster_mask,
                                 cells = cell_idsL[[x]],pairs=TRUE,
                                 directions = ngMat,include=FALSE)
    rcells <- unique(adj_cells[, 2])
    return(rcells)

  },.progress = progress)

  rcells <- unique(unlist(rcellsL))
  if (length(rcells) < n_bg) {
    n_bg  <- length(rcells)
  }
  sam <- sample(rcells, size = n_bg)
  return(sam)
}
