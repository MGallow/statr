# Matt Galloway

#' @title Time function
#' @description Simple function that prints the computation time of a function
#'
#' @param f the function to time
#' @return returns the elapsed time
#' @export
#' @examples
#' timeit(lm(dist ~ speed, cars))
#'
#'


# we define the timeit function this will print the computation time of function
timeit = function(f) {

    # set start time
    start = proc.time()

    # function to time
    fun = f

    # print elapsed time
    print(proc.time() - start)

}
