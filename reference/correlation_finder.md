# Function to find strong correlations within environmental predictors

This function identifies variables with strong correlations based on a
specified threshold.

## Usage

``` r
correlation_finder(
  environmental_data,
  method = "spearman",
  threshold,
  verbose = TRUE
)
```

## Arguments

- environmental_data:

  A matrix or data.frame containing environmental data.

- method:

  Method used to estimate the correlation matrix. Possible options
  include "spearman" (Spearman's rank correlation), "pearson" (Pearson's
  correlation), or "kendall" (Kendall's tau correlation).

- threshold:

  Correlation threshold value. Variables with absolute correlation
  values greater than or equal to this threshold are considered strongly
  correlated.

- verbose:

  Logical. If `TRUE`, prints verbose output detailing correlations.

## Value

A list with two elements:

- `not_correlated_vars`: A vector containing names of variables that are
  not strongly correlated.

- `correlation_values`: A list with correlation values for all pairs of
  variables.

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
future::plan("sequential")
envdata <- abex$env_data[,-ncol(abex$env_data)]
ecors <- tenm::correlation_finder(environmental_data =envdata,
                                  method="spearman",
                                  threshold = 0.7 )
#> Warning: the standard deviation is zero
#> *****************************************************************
#> 
#>  Here is a list of variables that can summarize your niche
#>  information, according to the threshold of 0.7 :
#> 
#>   bio_01 bio_02 bio_03 bio_04 bio_12 bio_14 bio_15 
#> 
#> *****************************************************************
#> 
#> ----------------------------------------------------------------
#> 
#> Correlation list:
#> 
#> Variable bio_01 is strongly correlated with:
#> 
#>    bio_01    bio_05    bio_06    bio_08    bio_09    bio_10    bio_11 
#> 1.0000000 0.9720096 0.8622147 0.9762799 0.9389854 0.9884480 0.9708369 
#> ----------------------------------------------------------------
#> 
#> Variable bio_02 is strongly correlated with:
#> 
#>   bio_02   bio_07 
#> 1.000000 0.750235 
#> ----------------------------------------------------------------
#> 
#> Variable bio_03 is strongly correlated with:
#> 
#> bio_03 
#>      1 
#> ----------------------------------------------------------------
#> 
#> Variable bio_04 is strongly correlated with:
#> 
#> bio_04 
#>      1 
#> ----------------------------------------------------------------
#> 
#> Variable bio_05 is strongly correlated with:
#> 
#>    bio_01    bio_05    bio_06    bio_08    bio_09    bio_10    bio_11 
#> 0.9720096 1.0000000 0.8287244 0.9468420 0.9218126 0.9528169 0.9538006 
#> ----------------------------------------------------------------
#> 
#> Variable bio_06 is strongly correlated with:
#> 
#>    bio_01    bio_05    bio_06    bio_08    bio_09    bio_10    bio_11 
#> 0.8622147 0.8287244 1.0000000 0.8314238 0.8774418 0.8339515 0.8896244 
#> ----------------------------------------------------------------
#> 
#> Variable bio_07 is strongly correlated with:
#> 
#>   bio_02   bio_07 
#> 0.750235 1.000000 
#> ----------------------------------------------------------------
#> 
#> Variable bio_08 is strongly correlated with:
#> 
#>    bio_01    bio_05    bio_06    bio_08    bio_09    bio_10    bio_11 
#> 0.9762799 0.9468420 0.8314238 1.0000000 0.9187488 0.9767573 0.9329921 
#> ----------------------------------------------------------------
#> 
#> Variable bio_09 is strongly correlated with:
#> 
#>    bio_01    bio_05    bio_06    bio_08    bio_09    bio_10    bio_11 
#> 0.9389854 0.9218126 0.8774418 0.9187488 1.0000000 0.9263277 0.9321939 
#> ----------------------------------------------------------------
#> 
#> Variable bio_10 is strongly correlated with:
#> 
#>    bio_01    bio_05    bio_06    bio_08    bio_09    bio_10    bio_11 
#> 0.9884480 0.9528169 0.8339515 0.9767573 0.9263277 1.0000000 0.9473264 
#> ----------------------------------------------------------------
#> 
#> Variable bio_11 is strongly correlated with:
#> 
#>    bio_01    bio_05    bio_06    bio_08    bio_09    bio_10    bio_11 
#> 0.9708369 0.9538006 0.8896244 0.9329921 0.9321939 0.9473264 1.0000000 
#> ----------------------------------------------------------------
#> 
#> Variable bio_12 is strongly correlated with:
#> 
#>    bio_12    bio_13    bio_16    bio_18    bio_19 
#> 1.0000000 0.8849691 0.9095685 0.8344122 0.7131398 
#> ----------------------------------------------------------------
#> 
#> Variable bio_13 is strongly correlated with:
#> 
#>    bio_12    bio_13    bio_16    bio_18 
#> 0.8849691 1.0000000 0.9462376 0.7578118 
#> ----------------------------------------------------------------
#> 
#> Variable bio_14 is strongly correlated with:
#> 
#>    bio_14    bio_17 
#> 1.0000000 0.7848596 
#> ----------------------------------------------------------------
#> 
#> Variable bio_15 is strongly correlated with:
#> 
#> bio_15 
#>      1 
#> ----------------------------------------------------------------
#> 
#> Variable bio_16 is strongly correlated with:
#> 
#>    bio_12    bio_13    bio_16    bio_18 
#> 0.9095685 0.9462376 1.0000000 0.7283047 
#> ----------------------------------------------------------------
#> 
#> Variable bio_17 is strongly correlated with:
#> 
#>    bio_14    bio_17    bio_19 
#> 0.7848596 1.0000000 0.8555739 
#> ----------------------------------------------------------------
#> 
#> Variable bio_18 is strongly correlated with:
#> 
#>    bio_12    bio_13    bio_16    bio_18 
#> 0.8344122 0.7578118 0.7283047 1.0000000 
#> ----------------------------------------------------------------
#> 
#> Variable bio_19 is strongly correlated with:
#> 
#>    bio_12    bio_17    bio_19 
#> 0.7131398 0.8555739 1.0000000 
#> ----------------------------------------------------------------
#> 
# }
```
