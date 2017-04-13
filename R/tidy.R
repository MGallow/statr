## Matt Galloway


#' @title Tidy
#' @description Tidy's R code and updates package documentation. Directly uses Yihui Xie's 'formatR' package
#'
#' @examples
#' tidy()

# no argument function to tidy code in package
tidy = function() {

    formatR::tidy_dir("R")
    devtools::document()

}
