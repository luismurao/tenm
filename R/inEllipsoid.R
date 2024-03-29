
#' inEllipsoid: Determine if a point is inside or outside an ellipsoid
#'
#' @description Determine if a point is inside or outside an ellipsoid.
#' @param centroid A numeric vector of centroids for each environmental variable
#' @param eShape Shape matrix of the ellipsoid (can be a covariance matrix or a
#' minimum volume ellipsoid).
#' @param env_data A data frame with the environmental training data.
#' @param level Proportion of points to be included in the ellipsoids. This
#' parameter is equivalent to the error (E) proposed by Peterson et al. (2008).
#' @return A data.frame with 2 columns. The first "in_Ellipsoid" binary response
#' with values 1 (inside the ellipsoid) and zeros (outside the ellipsoid); the
#' second "mh_dist" Mahalanobis distance to centroid.
#' @examples
#' \dontrun{
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
#' future::plan("multisession",workers=10)
#' abex <- tenm::ex_by_date(abtc,train_prop=0.7)
#' varcorrs <- tenm::correlation_finder(environmental_data = abex$env_data[,-ncol(abex$env_data)],
#'                                      method = "spearman",
#'                                      threshold = 0.8,
#'                                      verbose = FALSE)
#' future::plan("sequential")
#' mod <- tenm::cov_center(data = abex$env_data,
#'                         mve = TRUE,
#'                         level = 0.975,
#'                         vars = c("bio_05","bio_06","bio_12"))
#' in_elip <- tenm::inEllipsoid(centroid = mod$centroid,
#'                        eShape = mod$covariance,
#'                        env_data = abex$env_data[,c("bio_05","bio_06","bio_12")],
#'                        level = 0.975)
#' # 1 = Inside the ellipsoid; 0 = Outside the ellipsoid
#' print(in_elip)
#' }
#' @export

inEllipsoid <- function(centroid,eShape,env_data,level){

  mh_dist <- stats::mahalanobis(env_data,
                                center = centroid,
                                cov =eShape)
  in_Ellipsoid <- mh_dist <= stats::qchisq(level,
                                           length(centroid))
  in_Ellipsoid <- in_Ellipsoid*1
  in_Ellipsoid_mh <- data.frame(in_Ellipsoid,mh_dist )

  return(in_Ellipsoid_mh)
}
