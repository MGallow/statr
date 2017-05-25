# Matt Galloway

#' @title Time-It
#' @description Simple function that prints the computation time of a function
#'
#' @param f the function to time
#' @return returns the elapsed time
#' @export
#' @examples
#' timeit(lm(dist ~ speed, cars))


timeit = function(f) {
    
    # time function
    system.time(f)
    
}
