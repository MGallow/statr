## Matt Galloway


#' @title Linear
#' @description Computes the linear regression coefficient estimates (ridge-penalization and weights, optional)
#'
#' @param X matrix or data frame
#' @param y matrix or data frame of response values
#' @param lam optional tuning parameter for ridge regularization term. Defaults to 'lam = 0'
#' @param weights optional vector of weights for weighted least squares
#' @param intercept add column of ones if not already present. Defaults to TRUE
#' @param kernel use linear kernel to compute ridge regression coefficeients. Defaults to TRUE when p >> n
#' @return returns the coefficient estimates
#' @export
#' @examples
#'
#' Weighted ridge regression
#' library(dplyr)
#' X = dplyr::select(iris, -c(Species, Sepal.Length))
#' y = dplyr::select(iris, Sepal.Length)
#' linearr(X, y, lam = 0.1, weights = rep(1:150))
#'
#' Kernelized ridge regression
#' linearr(X, y, lam = 0.1, kernel = T)

linearr = function(X, y, lam = 0, weights = NULL, intercept = TRUE, 
    kernel = FALSE) {
    
    # checks
    n = dim(X)[1]
    if (is.null(weights)) {
        weights = rep(1, n)
    }
    if (length(weights) != n) 
        stop("weights must be length ", n)
    if (all(lam >= 0) == FALSE) 
        stop("lam must be nonnegative!")
    if ((kernel == TRUE) & (lam == 0)) 
        stop("must specify lam to use kernel!")
    
    
    # initialization
    X = as.matrix(X)
    y = as.matrix(y)
    
    # execute linearc
    linear = linearc(X, y, lam, weights, intercept, kernel)
    
    # add intercept name, if needed
    betas = linear$coefficients
    grads = gradient_linearc(betas, X, y, lam, weights, intercept)
    if (intercept) {
        b1 = as.matrix(betas[1])
        rownames(b1) = "intercept"
        betas = rbind(b1, as.matrix(betas[-1, ]))
        g1 = as.matrix(grads[1])
        rownames(g1) = "intercept"
        grads = rbind(g1, as.matrix(grads[-1, ]))
    }
    
    # generate MSE
    fit = predict_linearr(linear, as.matrix(X), y)
    
    returns = list(coefficients = betas, MSE = fit$MSE, gradient = grads)
    return(returns)
}


