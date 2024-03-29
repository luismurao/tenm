% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sp_temporal_data.R
\name{sp_temporal_data}
\alias{sp_temporal_data}
\title{Function to create a Species Temporal Data object (STD object).}
\usage{
sp_temporal_data(
  occs,
  longitude,
  latitude,
  sp_date_var,
  occ_date_format = "y",
  layers_date_format = "y",
  layers_by_date_dir,
  layers_ext = "*.tif$"
)
}
\arguments{
\item{occs}{A data.frame with information about
the occurrence records of the specie that is being modeled. It is fundamental
to count with exact geographical coordinates of longitude and latitude where
the specie was detected or at least the nearest, also a temporal column
indicating the time at which was the record.}

\item{longitude}{If occs is a data.frame the user must indicate the variable
name of longitude data.}

\item{latitude}{If occs is a data.frame the user must indicate the variable
name of latitude data.}

\item{sp_date_var}{A date variable indicating the date of each observation.
The name of the variable where is stored the date of each observation.}

\item{occ_date_format}{Occurrences date format. It is the format with which
the dates of the occurrence points are organized. The possible options are
"y" for years; "ym" for years and months; "ymd" for year, month and day;
"mdy" for month, day and year; "my" for month and year; "dmy" for day, month
and year.}

\item{layers_date_format}{Raster layers of environmental information data
format. The possible options are "y" for years; "ym" for years and months;
"ymd" for year, month and day; "mdy" for month, day and year; "my" for month
and year; "dmy" for day, month and year.}

\item{layers_by_date_dir}{A directory which has contain inside other folders
organized by date with the raster layers of environmental information.}

\item{layers_ext}{This is the extension or path of each raster layer archive.
In other words, this is the object where is stored each one of the paths
that leads to the location of the environmental raster layers into the inner
memory.}
}
\value{
Returns a  sp.temporal.modeling object (list) with the coordinates
of each occurrences points, the years of observation and the path to the
temporal layers.
}
\description{
This function creates an object of class sp.temporal.modeling
that contains a list with four attributes. Inside the list, there is a data
frame called temporal_df containing five columns: the first two columns are
coordinates of longitude and latitude that came from the occurrences
records, the third is the date variable at which the species were
registered, the fourth is the format of dates for each layer, organized in
a particular pattern, for example year/month/day, year/month, just year or
some other arrangement like month/year, month/year/day, day/month/year and
the fifth one is the path where is stored the bioclimatic layer
corresponding to each year. The other three variables present in the list
are character class objects: sp_date_var[1] is the name of the date
variable column available in the occurrences records, lon_lat_vars[2] are
both of the columns with the coordinates of longitude and latitude
correspondingly, layers_ext[1] is the final extension format of the
environmental information (“.tif$”).
}
\examples{
library(tenm)
#A data.frame with occurrences points information of Abronia graminea.
# See help(abronia)
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
}
