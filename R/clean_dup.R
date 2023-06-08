#' Function to clean duplicated longitude and latitude data
#' @description Clean up duplicated or redundant occurrence records whiches
#' present overlapping longitude and latitude geographical coordinates regarding
#' a referent system that the user can choose from fourth possible ways to
#' eliminate: distance threshold, per single pixel grain of resolution, by pixel
#' neighborhood, or both combined distance and pixel at the same time.This
#' function has as its main purpose to eliminate occurrence points
#' geographically spliced, using the information of the spatial position to
#' determine it, in order to depurate large clouds of occurrence points that
#' could cause an environmental overestimation in the future model, we call this
#' process "cleaning of spatial duplicated data". The user has to set a distance
#' in "units" that arose an overlapping area from a random choosen record where
#' any other occurrence within this is consider as a duplicated and then
#' eliminated, the so called threshold is the diameter of a circunference in
#' grades asuming the aproximation of 1° = 111.2km in the Equator. Also you can
#' use a raster mask with a determinate pixel grain resolution as a base sift;
#' see in arguments raster_mask and n_ngbs.
#' @param data A dataframe with longitude and latitude of occurrence records
#' belongings to some specie.
#' @param longitude A character vector of the column name "longitude" within
#' the dataframe.
#' @param latitude A character vector of the column name of "latitude" within
#' the dataframe.
#' @param threshold A numeric value representing the euclidean distance between
#' coordinates to be considered as a duplicate. Also it could be view as a
#' value of radio (r) that covers an area.
#' @param by_mask Logical. If TRUE the elimination of duplicates will be done
#' using a raster layer as a mask; If False the elimination of duplicates will
#' be done by the distance threshold.
#' @param raster_mask An object of class RasterLayer that will be used to clean
#' duplicates that are present in the same ID pixel.
#' @param n_ngbs Number of pixel neighbors. Remove duplicates depending on how
#' many pixels range you want, "0" is for eliminate duplicates in the same pixel,
#' that means just one record per single pixel of resolution,"1" is a neighborhood
#' with one-pixel length of 3 for 3 pixels of area, "2" is an 5 for 5 vicinity
#' with 25 pixel area, and so on depending on how much area you want to cover,
#' following the formule: 2*n_ngbs+1
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
