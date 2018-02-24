statr
================

Alternatively, check out the [vignette](https://htmlpreview.github.io/?https://github.com/MGallow/statr/blob/master/Vignette.html) or [manual](https://github.com/MGallow/statr/blob/master/statr.pdf).

Overview
--------

`statr` is a personal R package that I have created for organizational/convenience purposes. A (possibly incomplete) list of functions contained in the package can be found below:

-   `tidy()` tidy's R package code and updates documentation
-   `timeit()` prints the computation time of a function
-   `scatter()` creates a scatterplot using ggplot
-   `diagnostic()` creates diagnostic plots using ggplot (residual and QQ)
-   `dsearch()` is a dichotomous search algorithm for minimizing a univariate function
-   `bsearch()` is a bi-section search algorithm for minimizing a univariate function

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
#plot diagnostic plots
diagnostic(iris, Sepal.Length, Sepal.Width)
```

![](README_files/figure-markdown_github/unnamed-chunk-2-2.png)
