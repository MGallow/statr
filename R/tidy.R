## Matt Galloway


#' @title Tidy
#' @description tidys package R code and updates package documentation. Directly uses Yihui Xie's 'formatR' package.
#' @export

# notice that there is no argument
tidy = function() {
    
    formatR::tidy_dir("R")
    devtools::document()
    
    pkg = as.package(".")
    system(paste0("R CMD Rd2pdf --force ", shQuote(pkg$path)))
    
    devtools::build(vignettes = FALSE)
    devtools::reload()
    
}
