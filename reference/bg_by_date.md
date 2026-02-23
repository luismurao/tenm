# Function to obtain environmental background organized by date

Function to retrieve background data from occurrence records. The
background data is organized as a function of the dated environmental
data.

## Usage

``` r
bg_by_date(
  this_species,
  buffer_ngbs = NULL,
  buffer_distance = 1000,
  n_bg = 50000,
  process_ngbs_by = 100
)
```

## Arguments

- this_species:

  An object of class sp.temporal.env representing species occurrence
  data organized by date. See [`ex_by_date`](ex_by_date.md).

- buffer_ngbs:

  Number of pixel neighbors used to build the buffer around each
  occurrence point.

- buffer_distance:

  Distance (in the same units as raster layers) used to create a buffer
  around occurrence points to sample background data.

- n_bg:

  Number of background points to sample.

- process_ngbs_by:

  Numeric parameter to improve memory management. It process neighbor
  cells by a quantity specified by the user.

## Value

An object of class sp.temporal.bg containing background data organized
by date. The object is a list with the following components:

- "bg_df": A data.frame with columns for longitude, latitude, year,
  layer_date, layer_path, cell_ids_year, and environmental information.

- Other metadata relevant to background sampling.

## Details

This function retrieves background data around species occurrence
points, sampled based on the dated environmental data provided in
`this_species`. Background points are sampled within a buffer around
each occurrence point. The function returns an object of class
sp.temporal.bg, which contains background data organized by date. This
object is the input of the function
[`tenm_selection`](tenm_selection.md).

## Examples

``` r
# \donttest{
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
#This code is for running in parallel
future::plan("multisession",workers=2)
abex <- tenm::ex_by_date(this_species = abtc,train_prop=0.7)
abbg <- tenm::bg_by_date(this_species = abex,
                         buffer_ngbs=10,n_bg=50000)
future::plan("sequential")
# }
```
