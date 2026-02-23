# Function to thin longitude and latitude data

Cleans up duplicated or redundant occurrence records that present
overlapping longitude and latitude coordinates. Thinning can be
performed using either a geographical distance threshold or a pixel
neighborhood approach.

## Usage

``` r
clean_dup(
  data,
  longitude,
  latitude,
  threshold = 0,
  by_mask = FALSE,
  raster_mask = NULL,
  n_ngbs = 0,
  process_ngbs_by = 100,
  progress = TRUE
)
```

## Arguments

- data:

  A data.frame with longitude and latitude of occurrence records.

- longitude:

  A character vector indicating the column name of the "longitude"
  variable.

- latitude:

  A character vector indicating the column name of the "latitude"
  variable.

- threshold:

  A numeric value representing the distance threshold between
  coordinates to be considered duplicates. Units depend on whether
  `by_mask` is `T` or `F`. If `T`, the user needs to specify the number
  of pixels that define the neighborhood of duplicates (see n_ngbs
  parameter).

- by_mask:

  Logical. If `T`, the thinning process will use a raster layer as a
  mask for defining distance in pixel units.

- raster_mask:

  An object of class SpatRaster that serves as a reference to thin the
  occurrence data. Required if `by_mask` is `T`.

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
duplicates based on the specified criteria.

## Details

This function cleans up duplicated occurrences based on the specified
distance threshold. If `by_mask` is `T`, the distance is interpreted as
pixel distance using the provided raster_mask; otherwise, it is
interpreted as geographic distance.

## Examples

``` r
data(abronia)
tempora_layers_dir <- system.file("extdata/bio",package = "tenm")
tenm_mask <- terra::rast(file.path(tempora_layers_dir,"1939/bio_01.tif"))
# Clean duplicates without raster mask (just by distance threshold)
# First check the number of occurrence records
print(nrow(abronia))
#> [1] 106
# Clean duplicated records using a distance of ~ 18 km (0.1666667 grades)
ab_1 <- tenm::clean_dup(data =abronia,
                        longitude = "decimalLongitude",
                        latitude = "decimalLatitude",
                        threshold = terra::res(tenm_mask),
                        by_mask = FALSE,
                        raster_mask = NULL)
# Check number of records
print(nrow(ab_1))
#> [1] 10
# Clean duplicates using a raster mask
ab_2 <- tenm::clean_dup(data =abronia,
                        longitude = "decimalLongitude",
                        latitude = "decimalLatitude",
                        threshold = terra::res(tenm_mask)[1],
                        by_mask = TRUE,
                        raster_mask = tenm_mask,
                        n_ngbs = 1)
# Check number of records
print(nrow(ab_2))
#> [1] 9
```
