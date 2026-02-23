# Compute omission rate and statistical metrics for ellipsoid models.

Computes omission rate and statistical metrics for ellipsoid models
using environmental data.

## Usage

``` r
ellipsoid_omr(
  env_data,
  env_test = NULL,
  env_bg,
  cf_level,
  mve = TRUE,
  proc = FALSE,
  proc_iter = 100,
  rseed = TRUE
)
```

## Arguments

- env_data:

  A data frame containing the environmental data used for modeling.

- env_test:

  A data frame with environmental testing data. Default is NULL. If
  provided, the selection process includes p-values from a binomial
  test.

- env_bg:

  Environmental data sampled from the calibration area to compute the
  approximated prevalence of the model.

- cf_level:

  Proportion of points to be included in the ellipsoids. Equivalent to
  the error (E) proposed by Peterson et al. (2008).
  [doi:10.1016/j.ecolmodel.2007.11.008](https://doi.org/10.1016/j.ecolmodel.2007.11.008)
  .

- mve:

  Logical. If `TRUE`, computes a minimum volume ellipsoid using
  [`cov.rob`](https://rdrr.io/pkg/MASS/man/cov.rob.html) from the MASS
  package. If `FALSE`, uses the covariance matrix of the input data.

- proc:

  Logical. If `TRUE`, performs a partial ROC test.

- proc_iter:

  Numeric. Total number of iterations for the partial ROC bootstrap.

- rseed:

  Logical. If `TRUE`, sets a random seed for the partial ROC bootstrap.
  Default is `TRUE`.

## Value

A data.frame with the following columns:

- "fitted_vars": Names of variables that were fitted.

- "nvars": Number of fitted variables

- "om_rate_train": Omission rate of the training data.

- "non_pred_train_ids": Row IDs of non-predicted training data.

- "om_rate_test"': Omission rate of the testing data.

- "non_pred_test_ids": Row IDs of non-predicted testing data.

- "bg_prevalence": Approximated prevalence of the model (see details).

- "pval_bin": p-value of the binomial test.

- "pval_proc": p-value of the partial ROC test.

- "env_bg_paucratio": Environmental partial AUC ratio value.

- "env_bg_auc": Environmental AUC value.

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
edata <- abex$env_data
etrain <- edata[edata$trian_test=="Train",c("bio_05","bio_06","bio_12")]
etest <- edata[edata$trian_test=="Test",c("bio_05","bio_06","bio_12")]
bg <- abbg$env_bg[,c("bio_05","bio_06","bio_12")]
eor <- ellipsoid_omr(env_data=etrain,env_test=etest,env_bg=bg,
                     cf_level=0.975,proc=TRUE)
eor
#>                  fitted_vars nvars om_rate_train non_pred_train_ids
#> P_value bio_05,bio_06,bio_12     3        0.0625              18,28
#>         om_rate_test non_pred_test_ids bg_prevalence pval_bin pval_proc
#> P_value            0                        0.505534        0         0
#>         env_bg_paucratio env_bg_auc
#> P_value         1.406124    0.73792
# }
```
