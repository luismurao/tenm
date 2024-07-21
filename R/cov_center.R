
#' Function to compute the covariance matrix of an ellipsoid niche model.
#' @description
#' Computes the covariance matrix, niche centroid, volume, and other
#' ellipsoid parameter based on the values of niche variables from
#' occurrence points.
#' @param data A data.frame or matrix containing numeric values of variables
#'   used to model the niche.
#' @param mve Logical. If \code{TRUE}, computes a minimum volume ellipsoid using
#'   the \code{\link[MASS]{cov.mve}} function from the MASS package. If
#'   \code{FALSE}, uses the covariance matrix of the input data.
#' @param level Proportion of data to be used for computing the ellipsoid,
#'   applicable when mve is \code{TRUE}.
#' @param vars Vector specifying column indexes or names of variables in
#' the input data used to fit the ellipsoid model.
#' @return A list containing the following components:
#'   - `centroid`: Centroid (mean vector) of the ellipsoid.
#'   - `covariance_matrix`: Covariance matrix based on the input data.
#'   - `volume`: Volume of the ellipsoid.
#'   - `semi_axes_lengths`: Lengths of semi-axes of the ellipsoid.
#'   - `axis_coordinates`: Coordinates of ellipsoid axes.

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
#' abex <- tenm::ex_by_date(abtc,train_prop=0.7)
#' future::plan("sequential")
#' mod <- tenm::cov_center(data = abex$env_data,
#'                         mve = TRUE,
#'                         level = 0.975,
#'                         vars = c("bio_05","bio_06","bio_12"))
#' # Print model parameters
#' print(mod)
#' }
#' @export


cov_center <- function (data, mve = TRUE, level, vars = NULL)
{

  data <- data[, vars]
  if (mve) {

    NDquntil <- function(nD, level) {
      n <- floor(nD * level)
      if (n > nD)
        n <- nD
      return(n)
    }
    n <- NDquntil(dim(data)[1], level)
    cent_var <- MASS::cov.mve(data, quantile.used = n)
    centroid <- cent_var$center
    vari <- cent_var$cov
  }
  else {
    centroid <- colMeans(data)
    vari <- stats::cov(data)
  }
  # Ellipsoid volume
  ellip_vol <- function(n, axis_length) {
    term1 <- 2 * pi^(n/2)
    term2 <- n * gamma(n/2)
    term3 <- prod(axis_length)
    term4 <- (term1/term2) * term3
    return(term4)
  }

  sigmaI <- solve(vari)/stats::qchisq(level, df = dim(data)[2])
  sigIEigenSys <- eigen(sigmaI)
  sigIEval <- sigIEigenSys$values
  sigIEvec <- sigIEigenSys$vectors
  stds <- 1/sqrt(sigIEval)
  axis_length <- NULL
  for (i in 1:dim(sigmaI)[1]) {
    axis_length[i] <- stds[i] * 2
  }
  names(axis_length) <- letters[1:dim(vari)[1]]
  n <- dim(vari)[1]

  vol2 <- ellip_vol(n, axis_length/2)
  axis_coordinates <- list()
  for (i in 1:dim(vari)[1]) {
    assign(paste0("l", i, "_inf"),
           centroid - sigIEvec[,i] * stds[i])
    assign(paste0("l", i, "_sup"),
           centroid + sigIEvec[,i] * stds[i])
    coord_matrix <- matrix(c(eval(parse(text = paste0("l",
                                                      i, "_sup"))),
                             eval(parse(text = paste0("l", i, "_inf")))),
                           byrow = T, nrow = 2)
    colnames(coord_matrix) <- names(centroid)
    rownames(coord_matrix) <- paste0("vec_", 1:2)
    axis_coordinates[[i]] <- coord_matrix
  }
  return(list(centroid = centroid,
              covariance = vari,
              niche_volume = vol2,
              SemiAxis_length = axis_length/2,
              axis_coordinates = axis_coordinates))
}
