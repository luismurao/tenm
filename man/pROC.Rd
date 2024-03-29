% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PartialROC.R
\name{pROC}
\alias{pROC}
\title{Partial ROC calculation for Niche Models}
\usage{
pROC(
  continuous_mod,
  test_data,
  n_iter = 1000,
  E_percent = 5,
  boost_percent = 50,
  rseed = FALSE,
  sub_sample = TRUE,
  sub_sample_size = 1000
)
}
\arguments{
\item{continuous_mod}{a SpatRaster or a numeric vector of the ecological
niche model to be evaluated. If a numeric vector is provided it should
contain the values of the predicted suitability.}

\item{test_data}{A numerical matrix, data.frame, or a numeric vector. If it
is data.frame or matrix it should contain coordinates of the occurrences used
to test the ecological niche model to be evaluated; columns must be:
longitude and latitude. If numeric vector it should contain the values of the
predicted suitability.}

\item{n_iter}{(numeric) number of bootstrap iterations to be performed;
default = 1000.}

\item{E_percent}{(numeric) value from 0 to 100 that will be used as a
threshold (E); default = 5.}

\item{boost_percent}{(numeric) value from 0 to 100 representing the percent
of testing data to be used for performing the bootstrap process for
calculating the partial ROC; default = 50.}

\item{rseed}{Logical. Whether or not to set a random seed. Default FALSE.}

\item{sub_sample}{Logical. Indicates whether the test should run using a
subsample of size sub_sample_size. It is recommended for big rasters.}

\item{sub_sample_size}{Numeric. Size of the sample to be used for computing
pROC values.}
}
\value{
A data.frame containing the AUC values and AUC ratios calculated for
each iteration.
}
\description{
pROC applies partial ROC tests to continuous niche models.
}
\details{
Partial ROC is calculated following Peterson et al.
(2008; \doi{10.1016/j.ecolmodel.2007.11.008}).
This function is a modification of the PartialROC function, available
at \url{https://github.com/narayanibarve/ENMGadgets}.
}
\examples{
data(abronia)
suit_1970_2000 <- terra::rast(system.file("extdata/suit_1970_2000.tif",
                                          package = "tenm"))
print(suit_1970_2000)
proc_test <- tenm::pROC(continuous_mod = suit_1970_2000,
                        test_data = abronia[,c("decimalLongitude",
                                               "decimalLatitude")],
                        n_iter = 500, E_percent=5,
                        boost_percent=50)
print(proc_test$pROC_summary)

}
\references{
Peterson, A.T. et al. (2008) Rethinking receiver operating
characteristic analysis applications in ecological niche modeling.
Ecol. Modell., 213, 63–72. \doi{10.1016/j.ecolmodel.2007.11.008}
}
