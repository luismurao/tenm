% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ex_by_date.R
\name{ex_by_date}
\alias{ex_by_date}
\title{Extract environmental data by date}
\usage{
ex_by_date(this_species, train_prop = 0.7)
}
\arguments{
\item{this_species}{Species Temporal Data object see
\code{\link[tenm]{sp_temporal_data}}.}

\item{train_prop}{Numeric. Training proportion. This is the proportion at
which each record (occurrences points) will be used as training or testing.
A 0.7 train_prop indicates that 70 of the data will be used to train the
model, and the other 30 to model itself.}
}
\value{
An object of class sp.temporal.env that contains a temporal
data.frame (temporal_df) with the following columns: latitude, longitude,
year, layer_dates, layers_path, cell_ids_year, and environmental data.
}
\description{
Function to extract environmental data by date. This function
generates training and testing data sets using a random partition with a
definite proportion.
}
\examples{
\dontrun{
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
abex <- tenm::ex_by_date(this_species = abtc,
                         train_prop=0.7)
future::plan("sequential")
}
}
