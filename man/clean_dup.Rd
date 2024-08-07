% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clean_dup.R
\name{clean_dup}
\alias{clean_dup}
\title{Function to thin longitude and latitude data}
\usage{
clean_dup(
  data,
  longitude,
  latitude,
  threshold = 0,
  by_mask = FALSE,
  raster_mask = NULL,
  n_ngbs = 0
)
}
\arguments{
\item{data}{A data.frame with longitude and latitude of occurrence records.}

\item{longitude}{A character vector indicating the column name of the
"longitude" variable.}

\item{latitude}{A character vector indicating the column name of the
"latitude" variable.}

\item{threshold}{A numeric value representing the distance threshold between
coordinates to be considered duplicates. Units depend on whether
\code{by_mask} is \code{T} or \code{F}. If \code{T}, the user needs to specify the number
of pixels that define the neighborhood of duplicates (see n_ngbs parameter).}

\item{by_mask}{Logical. If \code{T}, the thinning process will use a raster layer
as a mask for defining distance in pixel units.}

\item{raster_mask}{An object of class SpatRaster that serves as a reference
to thin the occurrence data. Required if \code{by_mask} is \code{T}.}

\item{n_ngbs}{Number of pixels used to define the neighborhood matrix that
helps determine which occurrences are duplicates:
\itemize{
\item 0 removes occurrences within the same pixel, keeping one.
\item 1 considers duplicates all occurrences within a distance of one pixel.
\item n considers duplicates all occurrences within a distance of n pixels.
}}
}
\value{
Returns a data.frame with cleaned occurrence records, excluding
duplicates based on the specified criteria.
}
\description{
Cleans up duplicated or redundant occurrence records that present overlapping
longitude and latitude coordinates. Thinning can be performed using either a
geographical distance threshold or a pixel neighborhood approach.
}
\details{
This function cleans up duplicated occurrences based on the specified
distance threshold. If \code{by_mask} is \code{T}, the distance is interpreted as
pixel distance using the provided raster_mask; otherwise, it is interpreted
as geographic distance.
}
\examples{
data(abronia)
tempora_layers_dir <- system.file("extdata/bio",package = "tenm")
tenm_mask <- terra::rast(file.path(tempora_layers_dir,"1939/bio_01.tif"))
# Clean duplicates without raster mask (just by distance threshold)
# First check the number of occurrence records
print(nrow(abronia))
# Clean duplicated records using a distance of ~ 18 km (0.1666667 grades)
ab_1 <- tenm::clean_dup(data =abronia,
                        longitude = "decimalLongitude",
                        latitude = "decimalLatitude",
                        threshold = terra::res(tenm_mask),
                        by_mask = FALSE,
                        raster_mask = NULL)
# Check number of records
print(nrow(ab_1))
# Clean duplicates using a raster mask
ab_2 <- tenm::clean_dup(data =abronia,
                        longitude = "decimalLongitude",
                        latitude = "decimalLatitude",
                        threshold = terra::res(tenm_mask)[1],
                        by_mask = TRUE,
                        raster_mask = tenm_mask,
                        n_ngbs = 1)
# Check number of records
print(nrow(ab_2))
}
