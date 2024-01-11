  <!-- badges: start -->
[![R-CMD-check](https://github.com/luismurao/tenm/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/luismurao/tenm/actions/workflows/R-CMD-check.yaml)
[![Lifecycle: stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
<!-- badges: end -->

# tenm <a href="https://luismurao.github.io/tenm/"><img src="man/figures/logo.png" align="right" height="139" /></a>

An R package with a set of functions to calibrate time-specific ecological 
niche models. 

# Installation

```r
if (!require('devtools')) install.packages('devtools')
devtools::install_github('luismurao/tenm')
# If you want to build vignette, install pandoc before and then
devtools::install_github('luismurao/tenm',build_vignettes=TRUE)
```

## Acknowledgments

CONACYT Ciencia de Frontera CF-2023-I-1156.
