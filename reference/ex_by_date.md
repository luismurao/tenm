# Extract environmental data by date

Function to extract environmental data by date. It generates training
and testing datasets using a random partition with a specified
proportion.

## Usage

``` r
ex_by_date(this_species, train_prop = 0.7)
```

## Arguments

- this_species:

  Species Temporal Data object. See
  [`sp_temporal_data`](sp_temporal_data.md) for details.

- train_prop:

  Numeric. Proportion of data to use for training. For example, a
  train_prop of 0.7 means 70% of the data will be used for training and
  30% for testing.

## Value

An object of class sp.temporal.env that consists in a list of five
elements:

1.  "temporal_df": a temporal data.frame (temporal_df) with the
    following columns: latitude, longitude, year, layer_dates,
    layers_path, cell_ids_year, and environmental data.

2.  "sp_date_var": Name of date variable.

3.  "lon_lat_vars": Names of the longitude and latitude variables.

4.  "layers_ext": Environmental layers extension.

5.  "env_data": Environmental data of occurrences.

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
future::plan("multisession",workers=2)
abex <- tenm::ex_by_date(this_species = abtc,
                         train_prop=0.7)
future::plan("sequential")
# }
```
