
# R classes for tenmp package
# February 2022
# Version 0.3.2
# Licence GPL v3

#' S3 classes to organize data and results of \code{tenmp} objects
#' @importFrom methods new
#' @author Luis Osorio-Olvera
#' @exportClass sp.temporal.modeling
#' @export
#'

methods::setClass("sp.temporal.modeling",
                   representation(
                     temporal_df = "data.frame",
                     sp_date_var = "character",
                     lon_lat_vars = "character",
                     layers_ext = "character"
                   ),
                   prototype (
                   ),
                   validity = function(object)	{
                     return(TRUE)
                   })


#' S3 classes to organize data and results of \code{tenmp} objects
#' @importFrom methods new
#' @author Luis Osorio-Olvera
#' @exportClass sp.temporal.env
#' @export
#'

methods::setClass("sp.temporal.env",contains = c("sp.temporal.modeling"),
                  representation(
                    env_data = "data.frame"
                  ),
                  prototype (
                  ),
                  validity = function(object)	{
                    return(TRUE)
                  })



#' S3 classes to organize data and results of \code{tenmp} objects
#' @importFrom methods new
#' @author Luis Osorio-Olvera
#' @exportClass sp.temporal.bg
#' @export
#'

methods::setClass("sp.temporal.bg",contains = c("sp.temporal.modeling","sp.temporal.env"),
                  representation(
                    env_bg = "data.frame"
                  ),
                  prototype (
                  ),
                  validity = function(object)	{
                    return(TRUE)
                  })


#' S3 classes to organize data and results of \code{tenmp} objects
#' @importFrom methods new
#' @author Luis Osorio-Olvera
#' @exportClass sp.temporal.selection
#' @export
#'

methods::setClass("sp.temporal.selection",contains = c("sp.temporal.modeling",
                                                       "sp.temporal.env",
                                                       "sp.temporal.bg"),
                  representation(
                    mods_table = "data.frame"
                  ),
                  prototype (
                  ),
                  validity = function(object)	{
                    return(TRUE)
                  })
