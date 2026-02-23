# Helper function to randomly select cell IDs for generating environmental background data.

This function returns pixel IDs to be sampled for generating
environmental background data around species occurrence points.

## Usage

``` r
cells2samp(
  data,
  longitude,
  latitude,
  cell_ids = NULL,
  buffer_ngbs = 2,
  raster_mask,
  process_ngbs_by = 10,
  n_bg = 50000,
  progress = TRUE
)
```

## Arguments

- data:

  A data.frame containing longitude and latitude data of occurrence
  points.

- longitude:

  A character vector specifying the column name of longitude in 'data'.

- latitude:

  A character vector specifying the column name of latitude in 'data'.

- cell_ids:

  A numeric vector indicating the IDs of cells that serve as geographic
  centers for buffers. Default is NULL.

- buffer_ngbs:

  Number of neighboring pixels around occurrence points used to build
  the buffer for sampling.

- raster_mask:

  An object of class SpatRaster used to obtain pixel IDs.

- process_ngbs_by:

  Numeric parameter to improve memory management. It process neighbor
  cells by a quantity specified by the user.

- n_bg:

  Number of background pixels to sample.

- progress:

  Logical. If `TRUE`, show computation progress.

## Value

A numeric vector of cell IDs to be sampled for environmental background
data.

## Examples

``` r
# \donttest{
# cells to sample
data(abronia)
temporal_layer <- system.file("extdata/bio/2016/bio_01.tif",package = "tenm")
raster_mask <- terra::rast(temporal_layer)
set.seed(123)
samp_01 <- tenm::cells2samp(data = abronia,
                            longitude = "decimalLongitude",
                            latitude = "decimalLatitude",
                            cell_ids = NULL,
                            buffer_ngbs = 4,
                            raster_mask = raster_mask,
                            process_ngbs_by = 10,
                            n_bg = 50000,
                            progress =TRUE)
# Generete a sample using pixel IDs
samp_02 <- tenm::cells2samp(data = abronia,
                            longitude = NULL,
                            latitude = NULL,
                            cell_ids = c(256,290,326),
                            buffer_ngbs = 4,
                            raster_mask = raster_mask,
                            process_ngbs_by = 10,
                            n_bg = 50000,
                            progress =TRUE)
# }
```
