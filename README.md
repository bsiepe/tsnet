
<!-- README.md is generated from README.Rmd. Please edit that file -->

# tsnet

<!-- badges: start -->
<!-- badges: end -->

The goal of tsnet is to include helpful functions for dynamic network
modelling in psychology and surrounding fields. Currently, the core of
the package is the comparison of Bayesian gVAR models as estimated in
the `BGGM` package.

## Installation

You can install the development version of `tsnet` from
[GitHub](https://github.com/bsiepe/tsnet) with:

``` r
# install.packages("devtools")
devtools::install_github("bsiepe/tsnet")
```

## Getting Started

This is a basic example which shows you how to solve a common problem:

``` r
library(tsnet)
#> Registered S3 methods overwritten by 'BFpack':
#>   method               from
#>   get_estimates.lm     bain
#>   get_estimates.t_test bain
#> Registered S3 method overwritten by 'GGally':
#>   method from   
#>   +.gg   ggplot2
library(BGGM)

# Load data
data <- BGGM::ifit
```

## References

If you use the package, please cite the preprint that introduces the
package and the test:

Siepe, B.S. & Heck, D.W. (2023). Bayesian Estimation and Comparison of
Idiographic Network Models.
