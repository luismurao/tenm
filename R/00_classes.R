
# R classes for tenm package
# February 2022
# Version 0.3.2
# Licence GPL v3

#' S3 classes to organize data and results of \code{tenm} objects
#' @importFrom methods new
#' @author Luis Osorio-Olvera
#' @return
#' This object is a list comprising four elements: a) A data.frame
#' containing occurrence records and layer information. b) A character vector
#' specifying variable names. c) A character vector indicating the names of
#' longitude and latitude variables. d) A character vector denoting the
#' layers extension.
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


#' S3 classes to organize data and results of \code{tenm} objects
#' @importFrom methods new
#' @author Luis Osorio-Olvera
#' @exportClass sp.temporal.env
#' @return An object of class 'sp.temporal.env' inheriting information from
#' 'sp.temporal.env'. This object adds a data.frame of environmental values
#' associated with occurrence data.
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



#' S3 classes to organize data and results of \code{tenm} objects
#' @importFrom methods new
#' @author Luis Osorio-Olvera
#' @return An object of class 'sp.temporal.bg'. The object inherits information
#' from objects of classes 'sp.temporal.modeling' and 'sp.temporal.env'. This
#' class adds environmental background information.
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


#' S3 classes to organize data and results of \code{tenm} objects
#' @importFrom methods new
#' @author Luis Osorio-Olvera
#' @exportClass sp.temporal.selection
#' @return An object of class 'sp.temporal.selection'. This object inherits
#' information from objects of classes 'sp.temporal.modeling', 'sp.temporal.env'
#' and 'sp.temporal.bg'. The object stores the results of
#' the model calibration and selection process in a data.frame.
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
