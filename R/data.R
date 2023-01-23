#' Occurrence records of Abronia graminea
#'
#' A dataset containing occurrence records for Abronia graminea.
#' The data was downloaded from GBIF (GBIF, 2022).
#'
#' @format A data frame with 106 rows and 5 variables:
#' \describe{
#'   \item{species}{Scientific name of the species}
#'   \item{decimalLongitude}{Longitude}
#'   \item{decimalLatitude}{Latitude}
#'   \item{year}{Observation year}
#'   \item{gbif_doi}{DOI id for citing the dataset}
#'   ...
#' }
#' @source GBIF.org (22 February 2022) GBIF Occurrence Download \url{https://doi.org/10.15468/dl.teyjm9}
"abronia"


#' Potential distribution of abronia graminea
#'
#' A dataset containing a niche model for Abronia graminea.
#' The model was generated using the time-specific modeling approach implemented
#' in the \emph{tenm} package.
"suit_1970_2000"

#' Colors for plotting
#'
#' A string vector of colors for plotting the vignette example.

"colors"
