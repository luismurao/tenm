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
\item{continuous_mod}{A SpatRaster or numeric vector of the ecological
niche model to be evaluated. If a numeric vector is provided, it should
contain the values of the predicted suitability.}

\item{test_data}{A numerical matrix, data.frame, or numeric vector:
\itemize{
\item If data.frame or matrix, it should contain coordinates of the
occurrences used to test the ecological niche model.
Columns must be: longitude and latitude.
\item If numeric vector, it should contain the values of the predicted
suitability.
}}

\item{n_iter}{Number of bootstrap iterations to perform for partial ROC
calculations. Default is 1000.}

\item{E_percent}{Numeric value from 0 to 100 used as the threshold (E)
for partial ROC calculations. Default is 5.}

\item{boost_percent}{Numeric value from 0 to 100 representing the
percentage of testing data to use for bootstrap iterations in partial ROC.
Default is 50.}

\item{rseed}{Logical. Whether or not to set a random seed for
reproducibility. Default is \code{FALSE}.}

\item{sub_sample}{Logical. Indicates whether to use a subsample of
the test data. Recommended for large datasets.}

\item{sub_sample_size}{Size of the subsample to use for computing pROC
values when sub_sample is \code{TRUE}.}
}
\value{
A list of two elements:
\itemize{
\item "pROC_summary": a data.frame containing the mean
AUC value, AUC ratio calculated for each iteration and the p-value of the
test.
\item "pROC_results": a data.frame with four columns containing the AUC
(auc_model), partial AUC (auc_pmodel), partial AUC of the random model
(auc_prand) and the AUC ratio (auc_ratio) for each iteration.
}
}
\description{
Apply partial ROC tests to continuous niche models.
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
