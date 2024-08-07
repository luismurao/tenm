% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ellipsoid_projection.R
\name{ellipsoid_projection}
\alias{ellipsoid_projection}
\title{ellipsoid_projection: function to project an ellipsoid model}
\usage{
ellipsoid_projection(
  envlayers,
  centroid,
  covar,
  level = 0.95,
  output = "suitability",
  plot = TRUE,
  size,
  xlab1 = "niche var 1",
  ylab1 = "niche var 2",
  zlab1 = "S",
  alpha = 0.1,
  ...
)
}
\arguments{
\item{envlayers}{A SpatRaster object of the niche variables.}

\item{centroid}{A vector with the values of the centers of the ellipsoid
(see \code{\link[tenm]{cov_center}}).}

\item{covar}{The shape matrix (covariance) of the ellipsoid
(see \code{\link[tenm]{cov_center}}).}

\item{level}{The proportion of points  to be included inside the ellipsoid}

\item{output}{The output distance: two possible values "suitability" or
"mahalanobis". By default the function uses "suitability".}

\item{plot}{Logical If \code{TRUE} a plot of the niche will be shown.}

\item{size}{The size of the points of the niche plot.}

\item{xlab1}{For x label for 2-dimensional histogram}

\item{ylab1}{For y label for 2-dimensional histogram}

\item{zlab1}{For z label for 2-dimensional histogram}

\item{alpha}{Control the transparency of the 3-dimensional ellipsoid}

\item{...}{Arguments passed to \code{\link[rgl]{plot3d}} function from rgl}
}
\value{
Returns a SpatRaster of suitability values.
}
\description{
Function to project an ellipsoid model using the shape matrix
(covariance matrix) of the niche variables.
}
\examples{
\donttest{
library(tenm)
data("abronia")
tempora_layers_dir <- system.file("extdata/bio",package = "tenm")
abt <- tenm::sp_temporal_data(occs = abronia,
                              longitude = "decimalLongitude",
                              latitude = "decimalLatitude",
                              sp_date_var = "year",
                              occ_date_format="y",
                              layers_date_format= "y",
                              layers_by_date_dir = tempora_layers_dir,
                              layers_ext="*.tif$")
abtc <- tenm::clean_dup_by_date(abt,threshold = 10/60)
future::plan("multisession",workers=2)
abex <- tenm::ex_by_date(this_species = abtc,train_prop=0.7)
abbg <- tenm::bg_by_date(this_species = abex,
                         buffer_ngbs=10,n_bg=50000)
future::plan("sequential")
mod <- tenm::cov_center(data = abex$env_data,
                        mve = TRUE,
                        level = 0.975,
                        vars = c("bio_05","bio_06","bio_12"))
layers_path <-   list.files(file.path(tempora_layers_dir,
                                      "2016"),
                            pattern = ".tif$",full.names = TRUE)
elayers <- terra::rast(layers_path)
nmod <- ellipsoid_projection(envlayers = elayers[[names(mod$centroid)]],
                             centroid = mod$centroid,
                             covar = mod$covariance,
                             level = 0.99999,
                             output = "suitability",
                             size = 3,
                             plot = TRUE)
}

}
