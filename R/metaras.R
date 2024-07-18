#' Helper function to obtain layer name from a raster layer
#' @description Returns a character vector with the name of the raster layer.
#' @param r An object of class SpatRaster representing the raster layer.
#' @return A character vector with the name of the raster layer.
#' @examples
#' tempora_layers_dir <- system.file("extdata/bio",package = "tenm")
#' p1 <- list.files(tempora_layers_dir,full.names=TRUE,
#'                  pattern=".tif$",recursive=TRUE)[1]
#' r1 <- terra::rast(p1)
#' print(tenm::metaras(r1))
#' @export
metaras <- function(r){
  f1 <- paste0("r@",names(attributes(r))[1])
  f2 <- paste0(f1,"@.xData$names")
  f3 <- eval(parse(text = f2))
  return(f3)
  #np <- normalizePath()
  #sysinf <- Sys.info()
  #os <- sysinf['sysname']
  #if(os == "Windows") patt <- "\\\\" else patt <- "/"
  #vs <- unlist(strsplit(np,split = patt))
  #vs <- unlist(strsplit(vs[length(vs)],split="[.]"))
  #lname <- vs[1]
  #return(lname)
}
