#' Function to clean duplicated longitude and latitude data
#' @description Clean duplicated longitude and latitude data using threshold
#' distance. In a process that we call "cleaning of duplicated", this function
#' has the main propose of eliminating data occurrence points spliced, through
#' the information of longitude and latitude,identifying the redundant
#' environmental data in a given distance threshold.That means the possibility
#' of be able to debug large occurrence point clouds that could cause an
#' overestimation of the future model. This function also allows you to
#' eliminate duplicate occurrence points per pixel or in a given pixel
#' neighborhood.
#' @param data A data.frame with longitude and latitude data
#' @param longitude A character vector of the column name of longitude.
#' @param latitude A character vector of the column name of latitude.
#' @param threshold A numeric value representing the euclidean distance between
#' coordinates to be considered as a duplicate.
#' @param by_mask Logical. If TRUE the elimination of duplicates will be done
#' using a raster layer as a mask; If False the elimination of duplicates will
#' be done by the distance threshold.
#' @param raster_mask An object of class RasterLayer that will be used to clean
#' duplicates that are present in the same ID pixel.
#' @param n_ngbs Number of pixel neighbors. Remove duplicates depending on how
#' many pixels range you want, 1 is for eliminate duplicates in the same pixel,
#' 2 is a neighborhood of 3 for 3 pixels, 3 is an 5 for 5 vicinity and so on
#' depending on how much area you want to cover.
#' @return Returns a data.frame with coordinate data from species
#' @examples
#' data(abronia)
#' data(suit_1970_2000)
#' # Clean duplicates without raster mask (just by distance threshold)
#' # First check the number of occurrence records
#' print(nrow(abronia))
#' # Clean duplicated records using a distance of ~ 18 km (0.1666667 grades)
#' ab_1 <- tenm::clean_dup(data =abronia,
#'                         longitude = "decimalLongitude",
#'                         latitude = "decimalLatitude",
#'                         threshold = raster::res(suit_1970_2000),
#'                         by_mask = FALSE,
#'                         raster_mask = NULL)
#' # Check number of records
#' print(nrow(ab_1))
#' # Clean duplicates using a raster mask
#' ab_2 <- tenm::clean_dup(data =abronia,
#'                         longitude = "decimalLongitude",
#'                         latitude = "decimalLatitude",
#'                         threshold = raster::res(suit_1970_2000),
#'                         by_mask = TRUE,
#'                         raster_mask = suit_1970_2000,
#'                         n_ngbs = 0)
#' # Check number of records
#' print(nrow(ab_2))
#' @export
#'
clean_dup <- function(data,longitude,latitude,threshold=0.0, by_mask = FALSE,
                      raster_mask = NULL, n_ngbs = 0){
  data <- data[!is.na(data[,longitude]),]
  dat_sp <- sp::SpatialPointsDataFrame(data[,c(longitude ,latitude)],data)
  if(by_mask == TRUE && methods::is(raster_mask, "RasterLayer")){
    nbase <- 2*n_ngbs+1
    ngMat <- base::matrix(rep(1,nbase*nbase),
                          ncol =nbase,byrow = T )
    ngMat[n_ngbs+1,n_ngbs+1] <- 0
    cellids <- raster::cellFromXY(raster_mask, dat_sp)
    ids_nodup <- which(!duplicated (cellids))
    cellids2 <- cellids[ids_nodup]
    dat2 <- dat_sp@data [ids_nodup,]
    if (n_ngbs == 0){
      return(dat2)
    } else {
      dat2$cellid <- cellids2
      cellids2 <- sort(cellids[ids_nodup])
      adj_cells <- raster::adjacent(x = raster_mask,cells=cellids2,
                                    directions = ngMat)

      adj_cellsL <- split(adj_cells[,2], adj_cells[,1])
      targets <- names(adj_cellsL)
      keep <- rep(NA,length(targets ))
      j <- 1
      for(i in seq_along(targets)){
        focal <- targets[i]
        if(focal %in% names(adj_cellsL)){
          vecinos <- as.character(adj_cellsL [[focal]])
          id_bas <- which(names(adj_cellsL) %in% vecinos)
          if(length(id_bas)>0L) adj_cellsL <- adj_cellsL[- id_bas]
          keep[j] <- as.numeric(focal)
          j <- j+1
        }

      }
      ids_nodup2 <- which(dat2$cellid %in% keep)
      dat3 <- dat2[ids_nodup2,-ncol(dat2)]
      return(dat3)
    }

  }
  else{
    dat_sp1 <- sp::remove.duplicates(dat_sp, zero = threshold)
    return(dat_sp1@data)

  }
}
