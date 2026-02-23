# Helper function to obtain layer name from a raster layer

Returns a character vector with the name of the raster layer.

## Usage

``` r
metaras(r)
```

## Arguments

- r:

  An object of class SpatRaster representing the raster layer.

## Value

A character vector with the name of the raster layer.

## Examples

``` r
tempora_layers_dir <- system.file("extdata/bio",package = "tenm")
p1 <- list.files(tempora_layers_dir,full.names=TRUE,
                 pattern=".tif$",recursive=TRUE)[1]
r1 <- terra::rast(p1)
print(tenm::metaras(r1))
#> [1] "bio_01"
```
