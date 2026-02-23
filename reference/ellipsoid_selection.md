# ellipsoid_selection: Performs models selection for ellipsoid models

The function performs model selection for ellipsoid models using three
criteria: a) the omission rate, b) the significance of partial ROC and
binomial tests and c) the AUC value.

## Usage

``` r
ellipsoid_selection(
  env_train,
  env_test = NULL,
  env_vars,
  nvarstest,
  level = 0.95,
  mve = TRUE,
  env_bg = NULL,
  omr_criteria,
  parallel = FALSE,
  ncores = NULL,
  comp_each = 100,
  proc = FALSE,
  proc_iter = 100,
  rseed = TRUE
)
```

## Arguments

- env_train:

  A data frame with the environmental training data.

- env_test:

  A data frame with the environmental testing data. Default is NULL.

- env_vars:

  A vector with the names of environmental variables used in the
  selection process. To help choosing which variables to use see
  [`correlation_finder`](correlation_finder.md).

- nvarstest:

  A vector indicating the number of variables to fit the ellipsoids
  during model selection.

- level:

  Proportion of points to be included in the ellipsoids, equivalent to
  the error (E) proposed by Peterson et al. (2008).

- mve:

  Logical. If `TRUE`, a minimum volume ellipsoid will be computed. using
  [`cov.rob`](https://rdrr.io/pkg/MASS/man/cov.rob.html) from MASS. If
  `FALSE`, the covariance matrix of the input data will be used.

- env_bg:

  Environmental data to compute the approximated prevalence of the
  model, should be a sample of the environmental layers of the
  calibration area.

- omr_criteria:

  Omission rate criteria: the allowable omission rate for the selection
  process. Default is NULL (see details).

- parallel:

  Logical. If `TRUE`, computations will run in parallel. Default is `F`.

- ncores:

  Number of cores to use for parallel processing. Default uses all
  available cores minus one.

- comp_each:

  Number of models to run in each job in parallel computation. Default
  is 100.

- proc:

  Logical. If `TRUE`, a partial ROC test will be run.

- proc_iter:

  Numeric. Total iterations for the partial ROC bootstrap.

- rseed:

  Logical. If `TRUE`, set a random seed for partial ROC bootstrap.
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

- "mean_omr_train_test": Mean value of omission rates (train and test).

- "rank_by_omr_train_test": Rank value of importance in model selection
  by omission rate.

- "rank_omr_aucratio": Rank value by AUC ratio.

## Details

Model selection occurs in environmental space (E-space). For each
variable combination specified in nvarstest, the omission rate (omr) in
E-space is computed using [`inEllipsoid`](inEllipsoid.md) function.
Results are ordered by omr of the testing data. If env_bg is provided,
an estimated prevalence is computed and results are additionally ordered
by partial AUC. Model selection can be run in parallel. For more details
and examples go to [`ellipsoid_omr`](ellipsoid_omr.md) help.

## References

Peterson, A.T. et al. (2008) Rethinking receiver operating
characteristic analysis applications in ecological niche modeling. Ecol.
Modell. 213, 63–72.
[doi:10.1016/j.ecolmodel.2007.11.008](https://doi.org/10.1016/j.ecolmodel.2007.11.008)

## Author

Luis Osorio-Olvera <luismurao@gmail.com>

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
abex <- tenm::ex_by_date(this_species = abtc,train_prop=0.7)
abbg <- tenm::bg_by_date(this_species = abex,
                         buffer_ngbs=10,n_bg=50000)
future::plan("sequential")
varcorrs <- tenm::correlation_finder(environmental_data =
                                     abex$env_data[,-ncol(abex$env_data)],
                                     method = "spearman",
                                     threshold = 0.8,
                                     verbose = FALSE)
#> Warning: the standard deviation is zero
edata <- abex$env_data
etrain <- edata[edata$trian_test=="Train",] |> data.frame()
etest <- edata[edata$trian_test=="Test",] |> data.frame()
bg <- abbg$env_bg
res1 <- tenm::ellipsoid_selection(env_train = etrain,
                                  env_test = etest,
                                  env_vars = varcorrs$descriptors,
                                  nvarstest = 3,
                                  level = 0.975,
                                  mve = TRUE,
                                  env_bg = bg,
                                  omr_criteria = 0.1,
                                  parallel = FALSE,proc = TRUE)
#> -------------------------------------------------------------------
#>      **** Starting model selection process ****
#> -------------------------------------------------------------------
#> 
#> A total number of 84 models will be created for combinations of 9 variables taken by 3 
#> 
#> -------------------------------------------------------------------
#>   **A total number of 84 models will be tested **
#> 
#> -------------------------------------------------------------------
#>   54 models passed omr_criteria for train data
#>   9 models passed omr_criteria for test data
#>   9 models passed omr_criteria for train and test data
head(res1)
#>            fitted_vars nvars om_rate_train non_pred_train_ids om_rate_test
#> 1 bio_01,bio_04,bio_07     3       0.06250              18,31            0
#> 2 bio_01,bio_03,bio_04     3       0.06250               3,18            0
#> 3 bio_01,bio_02,bio_04     3       0.09375            3,18,31            0
#> 4 bio_01,bio_03,bio_12     3       0.06250               3,18            0
#> 5 bio_01,bio_07,bio_12     3       0.06250              18,31            0
#> 6 bio_04,bio_07,bio_12     3       0.06250              21,28            0
#>   non_pred_test_ids bg_prevalence pval_bin pval_proc env_bg_paucratio
#> 1                       0.4661502        0         0         1.486787
#> 2                       0.4692861        0         0         1.447409
#> 3                       0.4809998        0         0         1.442396
#> 4                       0.4851503        0         0         1.388009
#> 5                       0.4973252        0         0         1.380808
#> 6                       0.7312304        0         0         1.327986
#>   env_bg_auc mean_omr_train_test rank_by_omr_train_test rank_omr_aucratio
#> 1  0.8009887            0.031250                      1                 1
#> 2  0.7605788            0.031250                      2                 2
#> 3  0.7605862            0.046875                      6                 3
#> 4  0.6881975            0.031250                      3                 4
#> 5  0.7192737            0.031250                      4                 5
#> 6  0.6983763            0.031250                      5                 6
# }
```
