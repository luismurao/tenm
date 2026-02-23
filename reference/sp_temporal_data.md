# Function to create a Species Temporal Data object (STD object).

Creates an object of class sp.temporal.modeling that contains a list
with four attributes:

- temporal_df: A data frame with the following columns:

  - Longitude: Longitude coordinates of occurrence records.

  - Latitude: Latitude coordinates of occurrence records.

  - Date: Date variable indicating when the species were observed.

  - Layer Dates: Format of dates for each layer of environmental data.

  - Layers Path: Path to the bioclimatic layer corresponding to each
    year.

- sp_date_var: Name of the date variable column in the occurrence
  records.

- lon_lat_vars: Names of the longitude and latitude columns.

- layers_ext: Final extension format of the environmental information
  (e.g., ".tif").

## Usage

``` r
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
```

## Arguments

- occs:

  A data.frame with information about occurrence records of the species
  being modeled.

- longitude:

  Name of the variable in 'occs' containing longitude data.

- latitude:

  Name of the variable in 'occs' containing latitude data.

- sp_date_var:

  Name of the date variable.

- occ_date_format:

  Format of dates in occurrence records. Options: "y", "ym", "ymd",
  "mdy", "my", "dmy".

- layers_date_format:

  Format of dates in raster layers. Options: "y", "ym", "ymd", "mdy",
  "my", "dmy".

- layers_by_date_dir:

  Directory containing folders organized by date with raster layers of
  environmental information.

- layers_ext:

  Extension or path of each raster layer archive (e.g., ".tif").

## Value

Returns a sp.temporal.modeling object (list) with the coordinates of
each occurrences points, the years of observation and the path to the
temporal layers.

## Details

The format of dates for each layer can be organized in a particular
pattern, for example year/month/day ("ymd"), year/month ("ym"), just
year ("y") or some other arrangement like month/year ("my"),
month/year/day ("myd"), day/month/year ("dmy").

## Examples

``` r
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
```
