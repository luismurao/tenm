#' Helper function to obtain layer name from a
#' @description Returns a character vector with the name of the raster layer
#' @param r An object of class RasterLayer
#' @examples
#' tempora_layers_dir <- system.file("extdata/bio",package = "tenm")
#' p1 <- list.files(tempora_layers_dir,full.names=TRUE,
#'                  pattern=".tif$",recursive=TRUE)[1]
#' r1 <- raster::raster(p1)
#' print(tenm::metaras(r1))
#' @export
metaras <- function(r){
  np <- normalizePath(r@file@name)
  sysinf <- Sys.info()
  os <- sysinf['sysname']
  if(os == "Windows") patt <- "\\\\" else patt <- "/"
  vs <- unlist(strsplit(np,split = patt))
  vs <- unlist(strsplit(vs[length(vs)],split="[.]"))
  lname <- vs[1]
  return(lname)
}
