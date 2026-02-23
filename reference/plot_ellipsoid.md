# Function to plot ellipsoid models in E-space

The function plots 2D and 3D ellipsoids using environmental information
as coordinates.

## Usage

``` r
plot_ellipsoid(
  x,
  y,
  z = NULL,
  xlab = "x",
  ylab = "y",
  zlab = "x",
  mve = TRUE,
  level = 0.975,
  col = NULL,
  lwd_axes = 2,
  lty_axes = 2,
  semiaxes = FALSE,
  add = FALSE,
  ...
)
```

## Arguments

- x:

  Numeric vector representing the x coordinate of the ellipsoid.

- y:

  Numeric vector representing the y coordinate of the ellipsoid.

- z:

  Numeric vector representing the z coordinate of the ellipsoid.
  Defaults to NULL.

- xlab:

  Character vector with the name of the x-axis label.

- ylab:

  Character vector with the name of the y-axis label.

- zlab:

  Character vector with the name of the z-axis label (if plotting in
  3D).

- mve:

  Logical. If `TRUE`, fits a minimum volume ellipsoid model.

- level:

  Numeric value indicating the proportion of points to be included
  inside the ellipsoid model.

- col:

  Plot color

- lwd_axes:

  Line width for ellipsoid semi-axes.

- lty_axes:

  Line type for ellipsoid semi-axes.

- semiaxes:

  Logical. If `TRUE`, shows semi-axes of the ellipsoid.

- add:

  Logical. If `TRUE`, add plot to existing plot (for 2D plots only).

- ...:

  Additional arguments to pass to base::plot, rgl::plot3d, rgl::wire3d,
  or other plotting functions

## Value

A 2-dimensional or 3-dimensional plot depending on the input
coordinates.

## Examples

``` r
x <- rnorm(100)
y <- rnorm(100)
z <- rnorm(100)
# 2 dimensional plot
plot_ellipsoid(x, y, col = "darkgreen", xlab = "X-axis", ylab = "Y-axis",
               mve = TRUE, level = 0.95)

# 3 dimensional plot
plot_ellipsoid(x, y, z, col = "blue", xlab = "X-axis", ylab = "Y-axis",
               zlab = "Z-axis", mve = TRUE, level = 0.95)
# \donttest{
# Examples using functions of the package
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
abex <- tenm::ex_by_date(abtc,train_prop=0.7)
future::plan("sequential")
x <- abex$temporal_df$bio_05
y <- abex$temporal_df$bio_06
z <- abex$temporal_df$bio_12
# 2D ellipsoid
tenm::plot_ellipsoid(x = x, y=y, semiaxes= TRUE,xlim=c(140,390))
tenm::plot_ellipsoid(x = x+100, y=y, semiaxes= TRUE,add=TRUE)

# 3D ellipsoid
tenm::plot_ellipsoid(x = x, y=y, z=z ,semiaxes= FALSE)
tenm::plot_ellipsoid(x = x+100, y=y, z=z ,semiaxes= FALSE,add=TRUE)
# }
```
