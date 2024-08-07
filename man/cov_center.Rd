% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cov_center.R
\name{cov_center}
\alias{cov_center}
\title{Function to compute the covariance matrix of an ellipsoid niche model.}
\usage{
cov_center(data, mve = TRUE, level, vars = NULL)
}
\arguments{
\item{data}{A data.frame or matrix containing numeric values of variables
used to model the niche.}

\item{mve}{Logical. If \code{TRUE}, computes a minimum volume ellipsoid using
the \code{\link[MASS]{cov.mve}} function from the MASS package. If
\code{FALSE}, uses the covariance matrix of the input data.}

\item{level}{Proportion of data to be used for computing the ellipsoid,
applicable when mve is \code{TRUE}.}

\item{vars}{Vector specifying column indexes or names of variables in
the input data used to fit the ellipsoid model.}
}
\value{
A list containing the following components:
\itemize{
\item \code{centroid}: Centroid (mean vector) of the ellipsoid.
\item \code{covariance_matrix}: Covariance matrix based on the input data.
\item \code{volume}: Volume of the ellipsoid.
\item \code{semi_axes_lengths}: Lengths of semi-axes of the ellipsoid.
\item \code{axis_coordinates}: Coordinates of ellipsoid axes.
}
}
\description{
Computes the covariance matrix, niche centroid, volume, and other
ellipsoid parameter based on the values of niche variables from
occurrence points.
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
abex <- tenm::ex_by_date(abtc,train_prop=0.7)
future::plan("sequential")
mod <- tenm::cov_center(data = abex$env_data,
                        mve = TRUE,
                        level = 0.975,
                        vars = c("bio_05","bio_06","bio_12"))
# Print model parameters
print(mod)
}
}
