statr
================

[![Build
Status](https://travis-ci.org/MGallow/statr.svg?branch=master)](https://travis-ci.org/MGallow/statr)
[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/statr)](https://cran.r-project.org/package=statr)

## Overview

`statr` is a personal R package that I have created for
organizational/convenience purposes. This project is purely
experimental\! A (possibly incomplete) list of functions contained in
the package can be found below:

  - `tidy()` tidy’s R package code and updates documentation
  - `CVsplit()` splits data objects into training and testing sets
  - `data_gen()` generates data for linear regression settings
  - `dense()` generates dense matrices
  - `denseQR()` generates dense matrices via spectral decomposition
  - `tridiag()` generates tri-diagonal matrices
  - `derivative()` approximates the derivative for a given function
  - `diagnostic()` creates diagnostic plots using ggplot (residual and
    QQ)
  - `dsearch()` is a dichotomous search algorithm for minimizing a
    univariate function
  - `bsearch()` is a bi-section search algorithm for minimizing a
    univariate function
  - `LASSO()` calculates lasso regression coefficient with optimal
    tuning
  - `RIDGE()` calculates ridge regression coefficient with optimal
    tuning
  - `scatter()` creates a scatterplot using ggplot

See [vignette](https://mgallow.github.io/statr/) or
[manual](https://github.com/MGallow/statr/blob/master/statr.pdf).

## Installation

The easiest way to install is from the development version from Github:

``` r
# install.packages("devtools")
devtools::install_github("MGallow/statr")
```

If there are any issues/bugs, please let me know:
[github](https://github.com/MGallow/statr/issues). You can also contact
me via my [website](https://mgallow.github.io/). Contributions are
welcome\!

## Usage

``` r
library(statr)
library(magrittr)

# we will use the iris data set
X = dplyr::select(iris, -c(Species, Sepal.Length)) %>% as.matrix
y = dplyr::select(iris, Sepal.Length) %>% as.matrix
y_class = ifelse(dplyr::select(iris, Species) == "setosa", 1, 0)

# let us split the data for testing and training
CV = CVsplit(X, y)

# we can do some exploratory analysis
# plot Sepal.Length v Sepal.Width
iris %>% scatter(Sepal.Length, Sepal.Width)
```

![](README_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

``` r
# plot diagnostic plots
iris %>% diagnostic(Sepal.Length, Sepal.Width)
```

![](README_files/figure-gfm/unnamed-chunk-2-2.png)<!-- -->

``` r
# use the training data to fit ridge regression
RIDGE(CV$X.train, CV$Y.train)
```

    ## $betas
    ##                     [,1]
    ## Sepal.Width  0.006182588
    ## Petal.Length 0.008078388
    ## Petal.Width  0.002630500
    ## 
    ## $lam
    ## [1] 2233.731

``` r
# or lasso regression
statr::LASSO(CV$X.train, CV$Y.train)
```

    ## $betas
    ##                   [,1]
    ## Sepal.Width  0.0000000
    ## Petal.Length 0.3305527
    ## Petal.Width  0.0000000
    ## 
    ## $lam
    ## [1] 17.29497

``` r
# we can also generate our own data
data = data_gen(p = 10, r = 5, n = 100)
CV = CVsplit(data$X, data$Y)

# and again fit a ridge regression
statr::RIDGE(CV$X.train, CV$Y.train)
```

    ## $betas
    ##              s0           s0           s0           s0           s0
    ## V1   0.22361780 -0.115761718  0.104229847 -0.072600428 -0.154270994
    ## V2   0.02305304 -0.016964647  0.006759119 -0.041432524 -0.046074158
    ## V3   0.15874949  0.099007474  0.212525724  0.006138727  0.017029184
    ## V4   0.09539731 -0.002189673  0.064124189  0.092066074  0.026877249
    ## V5  -0.04285759 -0.041315878  0.076181904  0.104015778 -0.026374616
    ## V6   0.00579039  0.031399701  0.164809901  0.149226048 -0.031209118
    ## V7   0.04179391  0.033475039 -0.017651981  0.044214675 -0.135049107
    ## V8   0.24074144 -0.105872848 -0.214755207 -0.048318457 -0.076911694
    ## V9  -0.06354331 -0.082274759 -0.411792495 -0.054637518  0.055296303
    ## V10  0.11876787 -0.077535702 -0.048367468  0.057315345 -0.009928922
    ## 
    ## $lam
    ## [1] 0.7241429

``` r
# we can also generate random matrices with may be useful
# for other applications
# tridiagonal matrices
tridiag(p = 5)$Omega %>% round(5)
```

    ##          [,1]     [,2]     [,3]     [,4]     [,5]
    ## [1,]  1.96078 -1.37255  0.00000  0.00000  0.00000
    ## [2,] -1.37255  2.92157 -1.37255  0.00000  0.00000
    ## [3,]  0.00000 -1.37255  2.92157 -1.37255  0.00000
    ## [4,]  0.00000  0.00000 -1.37255  2.92157 -1.37255
    ## [5,]  0.00000  0.00000  0.00000 -1.37255  1.96078

``` r
# dense matrices
dense(p = 5)$Omega %>% round(5)
```

    ##          [,1]     [,2]     [,3]     [,4]     [,5]
    ## [1,]  8.04348 -1.95652 -1.95652 -1.95652 -1.95652
    ## [2,] -1.95652  8.04348 -1.95652 -1.95652 -1.95652
    ## [3,] -1.95652 -1.95652  8.04348 -1.95652 -1.95652
    ## [4,] -1.95652 -1.95652 -1.95652  8.04348 -1.95652
    ## [5,] -1.95652 -1.95652 -1.95652 -1.95652  8.04348

``` r
# compound symmetric matrices
compound(p = 5)$Omega %>% round(5)
```

    ##      [,1] [,2] [,3] [,4] [,5]
    ## [1,]  1.0  0.9  0.9  0.9  0.9
    ## [2,]  0.9  1.0  0.9  0.9  0.9
    ## [3,]  0.9  0.9  1.0  0.9  0.9
    ## [4,]  0.9  0.9  0.9  1.0  0.9
    ## [5,]  0.9  0.9  0.9  0.9  1.0
