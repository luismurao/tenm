#' Helper function to randomly select cell IDs for generating
#' environmental background data.
#' @description
#' This function returns pixel IDs to be sampled for generating environmental
#' background data around species occurrence points.
#' @param data A data.frame containing longitude and latitude data of occurrence
#' points.
#' @param longitude A character vector specifying the column name of longitude
#' in 'data'.
#' @param latitude A character vector specifying the column name of latitude
#' in 'data'.
#' @param cell_ids A numeric vector indicating the IDs of cells that serve as
#' geographic centers for buffers. Default is NULL.
#' @param buffer_ngbs Number of neighboring pixels around occurrence points
#' used to build the buffer for sampling.
#' @param raster_mask An object of class SpatRaster used to obtain pixel IDs.
#' @param process_ngbs_by Numeric parameter to improve memory management.
#' It process neighbor cells by a quantity specified by the user.
#' @param n_bg Number of background pixels to sample.
#' @param progress Logical. If \code{TRUE}, show computation progress.
#' @return A numeric vector of cell IDs to be sampled for environmental
#' background data.
#' @examples
#' \donttest{
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
#' # Generete a sample using pixel IDs
#' samp_02 <- tenm::cells2samp(data = abronia,
#'                             longitude = NULL,
#'                             latitude = NULL,
#'                             cell_ids = c(256,290,326),
#'                             buffer_ngbs = 4,
#'                             raster_mask = raster_mask,
#'                             process_ngbs_by = 10,
#'                             n_bg = 50000,
#'                             progress =TRUE)
#' }
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
