# inEllipsoid: Determine if a point is inside or outside an ellipsoid

Determine if a point is inside or outside an ellipsoid based on a
confidence level.

## Usage

``` r
inEllipsoid(centroid, eShape, env_data, level)
```

## Arguments

- centroid:

  Numeric vector of centroids for each environmental variable.

- eShape:

  Shape matrix of the ellipsoid (can be a covariance matrix or a minimum
  volume ellipsoid).

- env_data:

  Data frame with the environmental data.

- level:

  Proportion of points to be included in the ellipsoids, equivalent to
  the error (E) proposed by Peterson et al. (2008).

## Value

A data.frame with 2 columns:

- "in_Ellipsoid": Binary response indicating if each point is inside (1)
  or outside (0) the ellipsoid.

- "mh_dist": Mahalanobis distance from each point to the centroid.

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
abex <- tenm::ex_by_date(abtc,train_prop=0.7)
varcorrs <- tenm::correlation_finder(environmental_data = abex$env_data[,-ncol(abex$env_data)],
                                     method = "spearman",
                                     threshold = 0.8,
                                     verbose = FALSE)
#> Warning: the standard deviation is zero
future::plan("sequential")
mod <- tenm::cov_center(data = abex$env_data,
                        mve = TRUE,
                        level = 0.975,
                        vars = c("bio_05","bio_06","bio_12"))
in_elip <- tenm::inEllipsoid(centroid = mod$centroid,
                       eShape = mod$covariance,
                       env_data = abex$env_data[,c("bio_05","bio_06","bio_12")],
                       level = 0.975)
# 1 = Inside the ellipsoid; 0 = Outside the ellipsoid
print(in_elip)
#>    in_Ellipsoid    mh_dist
#> 1             1  1.5208724
#> 2             1  1.1441787
#> 3             1  4.5220005
#> 4             1  0.6103583
#> 5             1  0.8606538
#> 6             1  2.2926090
#> 7             1  2.7047455
#> 8             1  0.7001683
#> 9             1  1.4891237
#> 10            1  0.6497573
#> 11            1  0.9349719
#> 12            1  5.7577780
#> 13            1  0.6612829
#> 14            1  2.8157858
#> 15            1  5.0221566
#> 16            1  0.4349905
#> 17            1  0.9724952
#> 18            1  1.0877318
#> 19            1  0.1253176
#> 20            1  4.8374752
#> 21            1  1.9092557
#> 22            1  7.8798423
#> 23            1  1.8857259
#> 24            1  3.2271328
#> 25            1  6.6748831
#> 26            1  2.1050325
#> 27            1  3.4891962
#> 28            1  3.0973598
#> 29            1  0.8824323
#> 30            0  9.8155714
#> 31            1  1.7366118
#> 32            1  2.0110087
#> 33            0 13.4444583
#> 34            1  0.6095582
#> 35            1  4.3824177
#> 36            1  6.4840384
#> 37            1  5.6598968
#> 38            1  7.5065473
#> 39            1  3.8067799
#> 40            1  1.6922563
# }
```
