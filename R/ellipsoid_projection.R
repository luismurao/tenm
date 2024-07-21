#' ellipsoid_projection: function to project an ellipsoid model
#' @description Function to project an ellipsoid model using the shape matrix
#' (covariance matrix) of the niche variables.
#' @param envlayers A SpatRaster object of the niche variables.
#' @param centroid A vector with the values of the centers of the ellipsoid
#' (see \code{\link[tenm]{cov_center}}).
#' @param covar The shape matrix (covariance) of the ellipsoid
#' (see \code{\link[tenm]{cov_center}}).
#' @param level The proportion of points  to be included inside the ellipsoid
#' @param output The output distance: two possible values "suitability" or
#' "mahalanobis". By default the function uses "suitability".
#' @param plot Logical If \code{TRUE} a plot of the niche will be shown.
#' @param size The size of the points of the niche plot.
#' @param xlab1 For x label for 2-dimensional histogram
#' @param ylab1 For y label for 2-dimensional histogram
#' @param zlab1 For z label for 2-dimensional histogram
#' @param alpha Control the transparency of the 3-dimensional ellipsoid
#' @param ... Arguments passed to \code{\link[rgl]{plot3d}} function from rgl
#' @return Returns a SpatRaster of suitability values.
#' @export
#'
#' @examples
#' \donttest{
#' library(tenm)
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
#' abtc <- tenm::clean_dup_by_date(abt,threshold = 10/60)
#' future::plan("multisession",workers=2)
#' abex <- tenm::ex_by_date(this_species = abtc,train_prop=0.7)
#' abbg <- tenm::bg_by_date(this_species = abex,
#'                          buffer_ngbs=10,n_bg=50000)
#' future::plan("sequential")
#' mod <- tenm::cov_center(data = abex$env_data,
#'                         mve = TRUE,
#'                         level = 0.975,
#'                         vars = c("bio_05","bio_06","bio_12"))
#' layers_path <-   list.files(file.path(tempora_layers_dir,
#'                                       "2016"),
#'                             pattern = ".tif$",full.names = TRUE)
#' elayers <- terra::rast(layers_path)
#' nmod <- ellipsoid_projection(envlayers = elayers[[names(mod$centroid)]],
#'                              centroid = mod$centroid,
#'                              covar = mod$covariance,
#'                              level = 0.99999,
#'                              output = "suitability",
#'                              size = 3,
#'                              plot = TRUE)
#' }
#'
ellipsoid_projection <- function(envlayers,centroid,covar,level=0.95,
                                 output="suitability",plot=TRUE,size,
                                 xlab1="niche var 1",ylab1= "niche var 2",
                                 zlab1="S", alpha=0.1,...){

  if(methods::is(envlayers, "SpatRaster")){
    suitRaster <- envlayers[[1]]
    names(suitRaster) <- output
    nonaids <- which(!is.na(suitRaster[]))
    env_vars <- 1:terra::nlyr(envlayers) |> purrr::map_dfc(function(x){
      val <- envlayers[[x]][]
      dfv <- data.frame(val[nonaids])
      names(dfv) <- names(envlayers[[x]])
      return(dfv)
    })

  }
  else{
    stop("envlayers should be of class 'SpatRaster'")
  }
  # Calculating distance to the centroid
  mahalanobisD <- stats::mahalanobis(env_vars,
                                     center = centroid,
                                     cov = covar)

  suit <- function( mahalanobisD){
    expo <- exp(-0.5* mahalanobisD)
    return(expo)
  }
  # Computing the suitabilities
  if(output =="suitability"){
    suits <- suit( mahalanobisD)
  } else if(output == "mahalanobis"){
    suits <- mahalanobisD
  }

  rm(list = c("mahalanobisD"))
  suitVals <- rep(NA,terra::ncell(envlayers[[1]]))
  suitVals[nonaids] <- suits
  suitRaster[] <- suitVals
  rm(list=c("suitVals"))

  if(dim(env_vars)[2]==2 && plot==TRUE){

    x <- seq(from = centroid[1]/2,to =centroid[1]*2 ,length=100)
    x <- sort(x)
    y <- seq(from = centroid[2]/2,to =centroid[2]*2 ,length=100)
    y <- sort(y)

    suit1 <- function(x,y) {
      maha1 <- stats::mahalanobis(cbind(x,y),
                                  center = centroid,
                                  cov = covar)
      expo <- exp(-0.5* maha1)
      return(expo)
    }
    #z <- x %o% y
    z <- outer(x,y,FUN = suit1)
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(oldpar))
    p1 <- graphics::persp(x,y,z, box=T,xlab=xlab1,
                          ylab=ylab1,zlab=zlab1, col="blue",
                          theta = 55, phi = 30,r = 40,
                          d = 0.1, expand = 0.5,
                          ticktype = "detailed", nticks=5,
                          cex.lab=1.5, cex.axis=1.3,
                          cex.main=1.5, cex.sub=1.5)


    ranges <- t(sapply(list(x,y,z),range))
    means <- rowMeans(ranges)

    ## label offset distance, as a fraction of the plot width
    labelspace <- 0.12  ## tweak this until you like the result

    xpos <- min(x)-(diff(range(x)))*labelspace
    ypos <- min(y)-(diff(range(y)))*labelspace
    labelbot3d <- c(xpos,ypos,min(z))
    labeltop3d <- c(xpos,ypos,max(z))
    labelmid3d <- c(xpos,ypos,mean(range(z)))

    trans3dfun <- function(v) { grDevices::trans3d(v[1],v[2],v[3],p1) }
    labelbot2d <- trans3dfun(labelbot3d)
    labelmid2d <- trans3dfun(labelmid3d)
    labeltop2d <- trans3dfun(labeltop3d)
    labelang <- 180/pi*atan2(labeltop2d$y-labelbot2d$y,labeltop2d$x-labelbot2d$x)
    graphics::par(xpd=NA,srt=labelang)  ## disable clipping and set string rotation
    graphics::text(labelmid2d[1]$x,labelmid2d$y,zlab1,cex=1.5)


  }

  if(dim(env_vars)[2]==3 && plot==TRUE){

    data1 <- env_vars
    dfd <- dim(data1)[1] - 1
    dfn <- dim(data1)[2] - 1
    # Ellipsoid radius
    #ell.radius_E <- sqrt(dfn * qf(level, dfn, dfd))


    ellips_E  <- rgl::ellipse3d(covar,centre = centroid,level = 0.99)


    if(dfd > 50000)
      np <- 50000
    else
      np <- dim(data1)[1]

    toSam <- sample(1:length(data1[,1]),np)
    data1 <- data1[toSam,]
    if(output == "suitability"){
      suits2 <- suits[toSam]
    } else if(output == "mahalanobis"){

      suits2 <- suit(suits[toSam])
    }

    rgl::plot3d(data1,size = size,col=grDevices::hsv(suits2*.71,.95,.9),
                xlab = xlab1, ylab = ylab1, zlab = zlab1,...)
    rgl::wire3d(ellips_E, col=4, lit=FALSE,alpha=alpha,...)
  }
  rm(list=c("env_vars","suits"))

  return(suitRaster)


}
