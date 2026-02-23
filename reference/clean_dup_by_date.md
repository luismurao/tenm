# Function to thin occurrence data Cleans up duplicated longitude and latitude data by year using a specified distance threshold. The distance can be specified as a geographic distance or, if a raster_mask is provided, as a pixel distance.

Function to thin occurrence data Cleans up duplicated longitude and
latitude data by year using a specified distance threshold. The distance
can be specified as a geographic distance or, if a raster_mask is
provided, as a pixel distance.

## Usage

``` r
clean_dup_by_date(
  this_species,
  threshold,
  by_mask = FALSE,
  raster_mask = NULL,
  n_ngbs = 0,
  process_ngbs_by = 100,
  progress = TRUE
)
```

## Arguments

- this_species:

  An object of class sp.temporal.modeling representing species
  occurrence data organized by date. See
  [`sp_temporal_data`](sp_temporal_data.md).

- threshold:

  A numeric value representing the distance threshold between
  coordinates to be considered duplicates. Units depend on whether
  `by_mask` is `TRUE` or `FALSE`. If `TRUE`, the user needs to specify
  the number of pixels that define the neighborhood of duplicates (see
  n_ngbs parameter).

- by_mask:

  Logical. If `TRUE`, the thinning process will use a raster layer as a
  mask for defining distance in pixel units.

- raster_mask:

  An object of class SpatRaster that serves as a reference to thin the
  occurrence data. Required if `by_mask` is `TRUE`.

- n_ngbs:

  Number of pixels used to define the neighborhood matrix that helps
  determine which occurrences are duplicates:

  - 0 removes occurrences within the same pixel, keeping one.

  - 1 considers duplicates all occurrences within a distance of one
    pixel.

  - n considers duplicates all occurrences within a distance of n
    pixels.

- process_ngbs_by:

  Numeric parameter to improve memory management. It process neighbor
  cells by a quantity specified by the user.

- progress:

  Logical. If `TRUE`, show computation progress.

## Value

Returns a data.frame with cleaned occurrence records, excluding

An object of class sp.temporal.modeling containing a temporal data.frame
with cleaned occurrence data, including columns for longitude, latitude,
date variable, layers_dates, and layers_path.

## Details

This function is based on [`clean_dup`](clean_dup.md). It cleans up
duplicated occurrences based on the specified threshold. If `by_mask` is
`TRUE`, the distance is interpreted as pixel distance using the provided
raster_mask; otherwise, it is interpreted as geographic distance.

## Examples

``` r
library(tenm)
data("abronia")
tempora_layers_dir <- system.file("extdata/bio",package = "tenm")
tenm_mask <- terra::rast(file.path(tempora_layers_dir,"1939/bio_01.tif"))
# Clean duplicates without raster mask (just by distance threshold)
abt <- tenm::sp_temporal_data(occs = abronia,
                              longitude = "decimalLongitude",
                              latitude = "decimalLatitude",
                              sp_date_var = "year",
                              occ_date_format="y",
                              layers_date_format= "y",
                              layers_by_date_dir = tempora_layers_dir,
                              layers_ext="*.tif$")
abtc1 <- tenm::clean_dup_by_date(abt,threshold = terra::res(tenm_mask)[1])
# Check number of records
print(nrow(abtc1$temporal_df))
#> [1] 40
# Clean duplicates using a raster mask
abtc2 <- tenm::clean_dup_by_date(this_species = abt,
                                by_mask = TRUE,
                                threshold = terra::res(tenm_mask)[1],
                                raster_mask = tenm_mask,
                                n_ngbs = 0)
# Check number of records
print(nrow(abtc2$temporal_df))
#> [1] 50

abtc3 <- tenm::clean_dup_by_date(this_species = abt,
                                by_mask = TRUE,
                                threshold = terra::res(tenm_mask)[1],
                                raster_mask = tenm_mask,
                                n_ngbs = 2)
# Check number of records
print(nrow(abtc3$temporal_df))
#> [1] 38
```
