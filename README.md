statr
================

Alternatively, check out the [vignette](https://htmlpreview.github.io/?https://github.com/MGallow/statr/blob/master/Vignette.html) or [manual](https://github.com/MGallow/statr/blob/master/statr.pdf).

Overview
--------

`statr` is a personal R package that I have created for organizational/convenience purposes. A (possibly incomplete) list of functions contained in the package can be found below:

-   `tidy()` tidy's R package code and updates documentation
-   `timeit()` prints the computation time of a function
-   `scatter()` creates a scatterplot using ggplot
-   `dsearch()` is a dichotomous search algorithm for minimizing a univariate function
-   `bsearch()` is a bi-section search algorithm for minimizing a univariate function
-   `linearr()` computes the linear regression coefficient estimates (ridge regularization and weights optional)
-   `predict_linearr()` generates predictions and loss metric for linear regression
-   `logisticr()` computes the coefficient estimates for logistic regression (ridge and bridge regularization optional)
-   `predict_logisticr()` generates predictions and loss metrics for logistic

Installation
------------

The easiest way to install is from the development version from Github:

``` r
# install.packages("devtools")
devtools::install_github("MGallow/statr")
```

If there are any issues/bugs, please let me know: [github](https://github.com/MGallow/statr/issues). You can also contact me via my [website](http://users.stat.umn.edu/~gall0441/). Contributions are welcome!

Usage
-----

``` r
library(statr)

#we will use the iris data set
X = dplyr::select(iris, -c(Species, Sepal.Length))
y = dplyr::select(iris, Sepal.Length)
y_class = ifelse(dplyr::select(iris, Species) == "setosa", 1, 0)

#plot Sepal.Length v Sepal.Width
scatter(iris, Sepal.Length, Sepal.Width)
```

![](README_files/figure-markdown_github/unnamed-chunk-2-1.png)

``` r
#ridge regression
linearr(X, y, lam = 0.1)
```

    ## $coefficients
    ##                 [,1]
    ## intercept  1.8778524
    ##            0.6462400
    ##            0.7023063
    ##           -0.5415988

``` r
#ridge logistic regression (IRLS)
logisticr(X, y_class, lam = 0.1, penalty = "ridge")
```

    ## $coefficients
    ##                [,1]
    ## intercept  6.276283
    ##            1.540809
    ##           -3.641782
    ##           -1.630507
    ## 
    ## $MSE
    ## [1] 5.13471e-05
    ## 
    ## $log.loss
    ## [1] 0.3956525
    ## 
    ## $misclassification
    ## [1] 0
    ## 
    ## $total.iterations
    ## [1] 11
    ## 
    ## $gradient
    ##                   [,1]
    ## intercept 4.536561e-11
    ##           1.175446e-10
    ##           1.671454e-10
    ##           5.304043e-11

``` r
#ridge logistic regression (MM)
logisticr(X, y_class, lam = 0.1, penalty = "ridge", method = "MM")
```

    ## $coefficients
    ##                [,1]
    ## intercept  6.276226
    ##            1.540818
    ##           -3.641769
    ##           -1.630502
    ## 
    ## $MSE
    ## [1] 5.134801e-05
    ## 
    ## $log.loss
    ## [1] 0.3956563
    ## 
    ## $misclassification
    ## [1] 0
    ## 
    ## $total.iterations
    ## [1] 5459
    ## 
    ## $gradient
    ##                   [,1]
    ## intercept 1.831808e-06
    ##           5.390560e-06
    ##           9.993386e-06
    ##           3.516598e-06

``` r
#bridge logistic regression (MM)
fit = logisticr(X, y_class, lam = 0.1, alpha = 1.2, penalty = "bridge")
```

    ## [1] "using MM algorithm..."

``` r
fit
```

    ## $coefficients
    ##                 [,1]
    ## intercept 13.0803079
    ##            0.6170122
    ##           -5.7168400
    ##           -0.2300421
    ## 
    ## $MSE
    ## [1] 2.749702e-05
    ## 
    ## $log.loss
    ## [1] 0.1878654
    ## 
    ## $misclassification
    ## [1] 0
    ## 
    ## $total.iterations
    ## [1] 26021
    ## 
    ## $gradient
    ##                   [,1]
    ## intercept 1.790042e-06
    ##           5.268013e-06
    ##           9.997899e-06
    ##           3.518574e-06

``` r
#time bridge logistic regression
timeit(logisticr(X, y_class, lam = 0.1, alpha = 1.2, penalty = "bridge"))
```

    ## [1] "using MM algorithm..."

    ##    user  system elapsed 
    ##   0.235   0.000   0.235

``` r
#predict using bridge logistic regression estimates
predict_logisticr(fit, X[1:3,], y_class[1:3])
```

    ## $fitted.values
    ##           [,1]
    ## [1,] 0.9992467
    ## [2,] 0.9989747
    ## [3,] 0.9994881
    ## 
    ## $class
    ##      [,1]
    ## [1,]    1
    ## [2,]    1
    ## [3,]    1
    ## 
    ## $MSE
    ## [1] 6.269174e-07
    ## 
    ## $log.loss
    ## [1] 0.002291457
    ## 
    ## $misclassification
    ## [1] 0
