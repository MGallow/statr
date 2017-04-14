## Matt Galloway


#' @title Tidy
#' @description It tidys package R code and updates package documentation. Directly uses Yihui Xie's 'formatR' package.
#' @export
#' @examples
#' tidy()

# notice that there is no argument
tidy = function() {
    
    formatR::tidy_dir("R")
    devtools::document()
    
}
