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
#' @param n_bg Number of background pixels.
#' @export
cells2samp <- function(data,longitude,latitude,cell_ids = NULL,buffer_ngbs = 2,
                       raster_mask,n_bg = 50000){
  if(is.null(cell_ids)){
    data <- data.frame(data[,c(longitude,latitude)])
    cell_ids <- raster::cellFromXY(raster_mask,data[,c(longitude,latitude)])
  } else if(!is.numeric(cell_ids)){
    stop("Provide valid cell numbers")
  }

  nbase <- 2 * buffer_ngbs + 1
  ngMat <- base::matrix(rep(1, nbase * nbase), ncol = nbase, byrow = TRUE)
  ngMat[buffer_ngbs + 1, buffer_ngbs + 1] <- 0
  adj_cells <- raster::adjacent(x = raster_mask, cells = cell_ids,
                                directions = ngMat)
  rcells <- unique(adj_cells[, 2])
  if (length(rcells) < n_bg) {
    n_bg  <- length(rcells)
  }
  sam <- sample(rcells, size = n_bg)
  return(sam)
}
