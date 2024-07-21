#' Function to plot ellipsoid models in E-space
#' @description The function plots 2D and 3D ellipsoids using environmental
#' information as coordinates.
#' @param x Numeric vector representing the x coordinate of the ellipsoid.
#' @param y Numeric vector representing the y coordinate of the ellipsoid.
#' @param z Numeric vector representing the z coordinate of the ellipsoid.
#' Defaults to NULL.
#' @param col Plot color
#' @param xlab Character vector with the name of the x-axis label.
#' @param ylab Character vector with the name of the y-axis label.
#' @param zlab Character vector with the name of the z-axis label
#' (if plotting in 3D).
#' @param mve Logical. If \code{TRUE}, fits a minimum volume ellipsoid model.
#' @param level Numeric value indicating the proportion of points to be
#' included inside the ellipsoid model.
#' @param semiaxes Logical. If \code{TRUE}, shows semi-axes of the ellipsoid.
#' @param lwd_axes Line width for ellipsoid semi-axes.
#' @param lty_axes Line type for ellipsoid semi-axes.
#' @param add Logical. If \code{TRUE}, add plot to existing plot (for 2D plots only).
#' @param ... Additional arguments to pass to base::plot,
#' rgl::plot3d, rgl::wire3d, or other plotting functions
#' @return A 2-dimensional or 3-dimensional plot depending on the input
#' coordinates.
#'
#' @export
#' @examples
#' x <- rnorm(100)
#' y <- rnorm(100)
#' z <- rnorm(100)
#' # 2 dimensional plot
#' plot_ellipsoid(x, y, col = "darkgreen", xlab = "X-axis", ylab = "Y-axis",
#'                mve = TRUE, level = 0.95)
#' # 3 dimensional plot
#' plot_ellipsoid(x, y, z, col = "blue", xlab = "X-axis", ylab = "Y-axis",
#'                zlab = "Z-axis", mve = TRUE, level = 0.95)
#' \donttest{
#' # Examples using functions of the package
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
#' abex <- tenm::ex_by_date(abtc,train_prop=0.7)
#' future::plan("sequential")
#' x <- abex$temporal_df$bio_05
#' y <- abex$temporal_df$bio_06
#' z <- abex$temporal_df$bio_12
#' # 2D ellipsoid
#' tenm::plot_ellipsoid(x = x, y=y, semiaxes= TRUE,xlim=c(140,390))
#' tenm::plot_ellipsoid(x = x+100, y=y, semiaxes= TRUE,add=TRUE)
#' # 3D ellipsoid
#' tenm::plot_ellipsoid(x = x, y=y, z=z ,semiaxes= FALSE)
#' tenm::plot_ellipsoid(x = x+100, y=y, z=z ,semiaxes= FALSE,add=TRUE)
#' }

plot_ellipsoid <- function(x,y,z=NULL,xlab="x",ylab="y",zlab="x",mve=TRUE,
                           level=0.975,col=NULL,
                           lwd_axes=2,lty_axes=2,
                           semiaxes=FALSE,add=FALSE,...){
  # Color palette taken from the RColorBrewer package
  dark2 <- c("#1B9E77","#D95F02","#7570B3",
             "#E7298A","#66A61E","#E6AB02",
             "#A6761D","#666666","#ff7f00")

  e_data <- data.frame(x,y)
  if(is.numeric(z)) e_data <- data.frame(e_data,z)
  e_metadata <- tenm::cov_center(data = e_data,mve = mve,
                                 level = level,vars = names(e_data))
  ndim <- length(e_metadata$centroid)
  r1 <- stats::qchisq(level, df = ndim)
  radius <- sqrt(r1)
  sigma <- e_metadata$covariance
  centroid <- e_metadata$centroid
  axes_coordinates <- e_metadata$axis_coordinates
  oldpar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(oldpar))
  if(is.null(col)) col <- sample(dark2,1)
  if(ndim==2){

    r <- sigma[1, 2]
    theta <- seq(0,2*pi,by=0.005)
    scale <- sqrt(c(sigma[1, 1],
                    sigma[2,2]))
    if (scale[1] > 0) r <- r/scale[1]
    if (scale[2] > 0) r <- r/scale[2]
    r <- min(max(r,-1),1)
    d <- acos(r)
    xx <- radius * scale[1] * cos(theta + d/2) + centroid[1]
    yy <- radius * scale[2] * cos(theta - d/2) + centroid[2]
    if(add){
      graphics::lines(xx,yy,lwd=3,col=col,...)
    }
    else{
      plot(xx,yy,type="l",xlab=xlab,ylab=ylab,lwd=3,col=col,...)
    }


    if(semiaxes){

      graphics::segments(x0 = axes_coordinates[[1]][1,1],
                         y0 = axes_coordinates[[1]][1,2],
                         x1 = axes_coordinates[[1]][2,1],
                         y1 = axes_coordinates[[1]][2,2],
                         col="gray70",lwd=lwd_axes,lty_axes,...)
      graphics::segments(x0 = axes_coordinates[[2]][1,1],
                         y0 = axes_coordinates[[2]][1,2],
                         x1 = axes_coordinates[[2]][2,1],
                         y1 = axes_coordinates[[2]][2,2],
                         col="gray70",lwd=lwd_axes,
                         lty_axes,...)
    }
  }
  if(ndim ==3){

    rgl::plot3d(x,y,z=z,xlab=xlab,ylab=ylab,zlab=zlab,type="n",add=add,...)
    ellips_E <- rgl::ellipse3d(sigma, centre = centroid,
                               level = level)
    rgl::wire3d(ellips_E, col=col,lit = FALSE,...)

    if(semiaxes){

      rgl::segments3d(x = c(axes_coordinates[[1]][1,1],
                            axes_coordinates[[1]][2,1]),
                      y = c(axes_coordinates[[1]][1,2],
                            axes_coordinates[[1]][2,2]),
                      z = c(axes_coordinates[[1]][1,3],
                            axes_coordinates[[1]][2,3]),
                      lwd=3,...)

      rgl::segments3d(x = c(axes_coordinates[[2]][1,1],
                            axes_coordinates[[2]][2,1]),
                      y = c(axes_coordinates[[2]][1,2],
                            axes_coordinates[[2]][2,2]),
                      z = c(axes_coordinates[[2]][1,3],
                            axes_coordinates[[2]][2,3]),
                      lwd=3,...)


      rgl::segments3d(x = c(axes_coordinates[[3]][1,1],
                            axes_coordinates[[3]][2,1]),
                      y = c(axes_coordinates[[3]][1,2],
                            axes_coordinates[[3]][2,2]),
                      z = c(axes_coordinates[[3]][1,3],
                            axes_coordinates[[3]][2,3]),
                      lwd=3,...)
    }

  }

}
