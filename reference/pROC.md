# Partial ROC calculation for Niche Models

Apply partial ROC tests to continuous niche models.

## Usage

``` r
pROC(
  continuous_mod,
  test_data,
  n_iter = 1000,
  E_percent = 5,
  boost_percent = 50,
  rseed = FALSE,
  sub_sample = TRUE,
  sub_sample_size = 1000
)
```

## Arguments

- continuous_mod:

  A SpatRaster or numeric vector of the ecological niche model to be
  evaluated. If a numeric vector is provided, it should contain the
  values of the predicted suitability.

- test_data:

  A numerical matrix, data.frame, or numeric vector:

  - If data.frame or matrix, it should contain coordinates of the
    occurrences used to test the ecological niche model. Columns must
    be: longitude and latitude.

  - If numeric vector, it should contain the values of the predicted
    suitability.

- n_iter:

  Number of bootstrap iterations to perform for partial ROC
  calculations. Default is 1000.

- E_percent:

  Numeric value from 0 to 100 used as the threshold (E) for partial ROC
  calculations. Default is 5.

- boost_percent:

  Numeric value from 0 to 100 representing the percentage of testing
  data to use for bootstrap iterations in partial ROC. Default is 50.

- rseed:

  Logical. Whether or not to set a random seed for reproducibility.
  Default is `FALSE`.

- sub_sample:

  Logical. Indicates whether to use a subsample of the test data.
  Recommended for large datasets.

- sub_sample_size:

  Size of the subsample to use for computing pROC values when sub_sample
  is `TRUE`.

## Value

A list of two elements:

- "pROC_summary": a data.frame containing the mean AUC value, AUC ratio
  calculated for each iteration and the p-value of the test.

- "pROC_results": a data.frame with four columns containing the AUC
  (auc_model), partial AUC (auc_pmodel), partial AUC of the random model
  (auc_prand) and the AUC ratio (auc_ratio) for each iteration.

## Details

Partial ROC is calculated following Peterson et al. (2008;
[doi:10.1016/j.ecolmodel.2007.11.008](https://doi.org/10.1016/j.ecolmodel.2007.11.008)
). This function is a modification of the PartialROC function, available
at <https://github.com/narayanibarve/ENMGadgets>.

## References

Peterson, A.T. et al. (2008) Rethinking receiver operating
characteristic analysis applications in ecological niche modeling. Ecol.
Modell., 213, 63–72.
[doi:10.1016/j.ecolmodel.2007.11.008](https://doi.org/10.1016/j.ecolmodel.2007.11.008)

## Examples

``` r
data(abronia)
suit_1970_2000 <- terra::rast(system.file("extdata/suit_1970_2000.tif",
                                          package = "tenm"))
print(suit_1970_2000)
#> class       : SpatRaster 
#> size        : 23, 20, 1  (nrow, ncol, nlyr)
#> resolution  : 0.1666667, 0.1666667  (x, y)
#> extent      : -99.16667, -95.83333, 17.16667, 21  (xmin, xmax, ymin, ymax)
#> coord. ref. : lon/lat WGS 84 (EPSG:4326) 
#> source      : suit_1970_2000.tif 
#> name        : suit_1970_2000 
#> min value   :   4.614139e-37 
#> max value   :   8.332142e-01 
proc_test <- tenm::pROC(continuous_mod = suit_1970_2000,
                        test_data = abronia[,c("decimalLongitude",
                                               "decimalLatitude")],
                        n_iter = 500, E_percent=5,
                        boost_percent=50)
print(proc_test$pROC_summary)
#>              Mean_AUC Mean_pAUC_ratio_at_5%               P_value 
#>             0.9024062             1.5066030             0.0000000 
```
