## Matt Galloway Augmented from Adam Rothman's STAT 8054
## code




#' @title Derivative
#' @description Takes the approximate derivative for a given function
#'
#' @param g the derivative of the function to minimize, where dg(u, ...) is the function evaluated at u.
#' @param x value to evaluate the derivative at
#' @param delta defaults to 10e-8
#' @export

derivative = function(g, x, delta = 1e-07) {
    
    # define derivate
    deriv = (g(x + delta) - g(x))/delta
    
    return(deriv)
    
}



##--------------------------------------------




#' @title Bisection search
#' @description Minimizes a univariate strictly pseudoconvex function over the interval [a, b]. This is augmented code from Adam Rothman's STAT 8054 course (2017).
#'
#' @param dg the derivative of the function to minimize, where dg(u, ...) is the function evaluated at u.
#' @param a left endpoint of the initial interval of uncertainty.
#' @param b right endpoint of the initial interval of uncertainty.
#' @param L the maximum length of the final interval of uncertainty.
#' @param quiet should the function stay quiet?
#' @param ... additional argument specifications for dg
#' @return returns the midpoint of the final interval of uncertainty.
#' @export

bsearch = function(dg, a, b, L = 1e-07, quiet = FALSE) {
    
    # initial estimate
    est = mean(c(a, b))
    
    # continue until interval smaller than L
    while (b - a > L) {
        
        # compute gradient at midpoint
        dgm = dg(est)
        
        # if gradient less than 0...
        if (dgm < 0) {
            
            # function is decreasing at est ## new interval is [est,
            # b]
            a = est
            
            # if gradient great than 0...
        } else if (dgm > 0) {
            
            # function is increasing at est ## new interval is [a,
            # mm]
            b = est
            
            # if gradient equal to 0...
        } else {
            
            # est is a stationary point
            b = est
            a = est
        }
        
        # create new interval and calculate estimate
        if (!quiet) 
            cat("new interval is", a, b, "\n")
        est = mean(c(a, b))
        
    }
    
    return(est)
    
}





##--------------------------------------------

#' @title Dichotomous search
#' @description Minimizes a univariate strictly quasiconvex function over the interval [a, b]. This is augmented code from Adam Rothman's STAT 8054 course (2017).
#'
#' @param g the function to minimize, where g(u, ...) is the function evaluated at u.
#' @param a left endpoint of the initial interval of uncertainty.
#' @param b right endpoint of the initial interval of uncertainty.
#' @param L the maximum length of the final interval of uncertainty.
#' @param eps search parameter, must be less than L/2
#' @param quiet should the function stay quiet?
#' @param ... additional argument specifications for g
#' @return returns the midpoint of the final interval of uncertainty.
#' @export

dsearch = function(g, a, b, L = 1e-07, eps = (L/2.1), quiet = FALSE) {
    
    # initial estimate
    est = mean(c(a, b))
    
    # continue until interval smaller than L
    while (b - a > L) {
        
        # calulcate direction of decreasing function
        lam = est - eps
        mu = est + eps
        g.at.lam = g(lam)
        g.at.mu = g(mu)
        
        # if increasing...
        if (g.at.lam < g.at.mu) {
            
            # cut off new interval at mu
            b = mu
            
            # if decreasing...
        } else if (g.at.lam > g.at.mu) {
            
            # cut off new interval at lam
            a = lam
            
            # if flat...
        } else {
            
            # create new interval
            b = mu
            a = lam
            
        }
        
        # create new interval and calculate new estimate
        if (!quiet) 
            cat("new interval is", a, b, "\n")
        est = mean(c(a, b))
        
    }
    
    return(est)
    
}
