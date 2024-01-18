#' Function to clean duplicated longitude and latitude data
#' @description Clean up duplicated or redundant occurrence records that
#' present overlapping longitude and latitude geographical coordinates regarding
#' a referent system that the user can choose from fourth possible ways to
#' eliminate: distance threshold, per single pixel grain of resolution, by pixel
#' neighborhood, or both combined distance and pixel at the same time.This
#' function has as its main purpose to eliminate occurrence points
#' geographically spliced, using the information of the spatial position to
#' determine it, in order to debug large clouds of occurrence points that
#' could cause an environmental overestimation in the future model, we call this
#' process "cleaning of spatial duplicated data". The user has to set a distance
#' in "units" that arose an overlapping area from a random chosen record where
#' any other occurrence within this is consider as a duplicated and then
#' eliminated. Also you can use a raster mask with a determinate pixel grain
#' resolution as a base sift; see in arguments raster_mask and n_ngbs.
#' @param data A data.frame with longitude and latitude of occurrence records
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
#' @param raster_mask An object of class SpatRaster that will be used to clean
#' duplicates that are present in the same ID pixel.
#' @param n_ngbs Number of pixel neighbors. Remove duplicates depending on how
#' many pixels range you want, 1 is for eliminate duplicates in the same pixel,
#' that means just one record per single pixel of resolution,2 is a neighborhood
#' with one-pixel length of 3 for 3 pixels, 3 is an 5 for 5 vicinity and so on
#' depending on how much area you want to cover.
#' @return Returns a data.frame with coordinate data from a species
#' @examples
#' data(abronia)
#' tempora_layers_dir <- system.file("extdata/bio",package = "tenm")
#' tenm_mask <- terra::rast(file.path(tempora_layers_dir,"1939/bio_01.tif"))
#' # Clean duplicates without raster mask (just by distance threshold)
#' # First check the number of occurrence records
#' print(nrow(abronia))
#' # Clean duplicated records using a distance of ~ 18 km (0.1666667 grades)
#' ab_1 <- tenm::clean_dup(data =abronia,
#'                         longitude = "decimalLongitude",
#'                         latitude = "decimalLatitude",
#'                         threshold = terra::res(tenm_mask),
#'                         by_mask = FALSE,
#'                         raster_mask = NULL)
#' # Check number of records
#' print(nrow(ab_1))
#' # Clean duplicates using a raster mask
#' ab_2 <- tenm::clean_dup(data =abronia,
#'                         longitude = "decimalLongitude",
#'                         latitude = "decimalLatitude",
#'                         threshold = terra::res(tenm_mask)[1],
#'                         by_mask = TRUE,
#'                         raster_mask = tenm_mask,
#'                         n_ngbs = 1)
#' # Check number of records
#' print(nrow(ab_2))
#' @export
#'
clean_dup <- function(data,longitude,latitude,threshold=0.0, by_mask = FALSE,
                      raster_mask = NULL, n_ngbs = 0){
  data <- data[!is.na(data[,longitude]),]


  #dat_sp <- sp::SpatialPointsDataFrame(data[,c(longitude ,latitude)],data)
  if(by_mask == TRUE && methods::is(raster_mask, "SpatRaster")){
    dat_sp <- sf::st_as_sf(data,coords=c(longitude,latitude),
                           crs=sf::st_crs(raster_mask))
    nbase <- 2*n_ngbs+1
    ngMat <- base::matrix(rep(1,nbase*nbase),
                          ncol =nbase,byrow = T )
    ngMat[n_ngbs+1,n_ngbs+1] <- 0
    cellids <- terra::cellFromXY(raster_mask, sf::st_coordinates(dat_sp))
    ids_nodup <- which(!duplicated (cellids))
    cellids2 <- cellids[ids_nodup]
    #dat2 <- dat_sp@data [ids_nodup,]
    coo <- data.frame(sf::st_coordinates(dat_sp))
    names(coo) <- c(longitude,latitude)
    dat2 <- sf::st_drop_geometry(dat_sp)
    dat2 <- data.frame(coo,dat2)
    if(ncol(dat2)>=4){
      dat2 <- dat2[ids_nodup,c(3,1,2,(4:ncol(dat2)))]
    } else{
      dat2 <- dat2[ids_nodup,c(3,1,2)]
    }
    if (n_ngbs == 0){
      return(dat2)
    } else {
      dat2$cellid <- cellids2
      cellids2 <- sort(cellids2)
      adj_cells <- terra::adjacent(x = raster_mask,cells=cellids2,
                                   directions = ngMat,
                                   pairs = TRUE)
      occ_adj_id <- which(adj_cells[,2] %in% cellids2)
      if(length(occ_adj_id)>0L){
        adj_cells <- adj_cells[occ_adj_id,]
      }

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
    dat_sp <- sf::st_as_sf(data,coords=c(longitude,latitude))
    dat_bf <- sf::st_buffer(dat_sp, dist = threshold[1])
    dup_mat <- sf::st_intersects(dat_sp,dat_bf,sparse = FALSE)
    #dup_mat <- sf::st_is_within_distance(dat_sp,dist = threshold[1],
    #                                     sparse = FALSE)
    diag(dup_mat) <- FALSE
    dup_mat[upper.tri(dup_mat)] <- FALSE
    ids_dups <- which(dup_mat,arr.ind=TRUE)
    ids_dupsL <- split(ids_dups[,2],ids_dups[,1],drop = FALSE)
    ids_dups <- unique(unlist(ids_dupsL))
    no_dups <- seq_len(nrow(data))
    if(length(ids_dups)>0L){
      no_dups <- no_dups[-ids_dups]
    }
    return(data[no_dups,])
  }
}
