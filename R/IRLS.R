## Matt Galloway


#' @title Logit
#' @description Computes the logit for u
#'
#' @param u some number. Ex: X %*% beta
#' @return returns the logit of u
#' @examples
#' logit(X %*% beta)


logitr = function(u) {
    
    # calculate logit probabilities
    exp(u)/(1 + exp(u))
    
}



##------------------------------------------------------------------------------------


#' @title Gradient of Logistic Regression (IRLS)
#' @description Computes the gradient of logistic regression (optional ridge regularization term). We use this to determine if the KKT conditions are satisfied. This function is to be used with the 'IRLS' function.
#'
#' @param betas beta estimates (includes intercept)
#' @param X matrix or data frame
#' @param y response vector of 0,1
#' @param lam tuning parameter for ridge regularization term
#' @param vec vector to specify which coefficients will be penalized
#' @return returns the gradient
#' @examples
#' gradient_IRLS_logistic(betas, X, y, lam = 0.1, penalty = 'ridge')

gradient_IRLS_logistic = function(betas, X, y, lam = 0, 
    vec) {
    
    # gradient for beta
    t(X) %*% (logitr(X %*% betas) - y) + lam * as.matrix(vec) * 
        betas
    
}



##--------------------------------------------------------------------------------------------


#' @title Iterative Re-Weighted Least Squares
#' @description Computes the logistic regression coefficient estimates using the iterative re-weighted least squares (IRLS) algorithm. This function is to be used with the 'logisticr' function.
#'
#' @param betas beta estimates (includes intercept)
#' @param X matrix or data frame
#' @param y matrix or vector of response 0,1
#' @param lam tuning parameter for regularization term
#' @param vec optional vector to specify which coefficients will be penalized
#' @param intercept Defaults to TRUE
#' @param tol tolerance - used to determine algorithm convergence
#' @param maxit maximum iterations
#' @return returns beta estimates (includes intercept), total iterations, and gradients.
#' @examples
#' IRLS(X, y, n.list = c(rep(1, n)), lam = 0.1, alpha = 1.5)


# calculates the coefficient estimates for logistic
# regression (IRLS)
IRLS = function(X, y, lam = 0, intercept = TRUE, tol = 10^(-5), 
    maxit = 1e+05, vec) {
    
    # initialize
    n = dim(X)[1]
    p = dim(X)[2]
    X = as.matrix(X)
    y = as.matrix(y)
    betas = as.matrix(rep(0.1, p))/n
    weights = rep(1, n)
    iteration = 1
    grads = gradient_IRLS_logistic(betas, X, y, lam, vec)
    
    # IRLS algorithm
    while ((iteration < maxit) & (max(abs(grads)) > tol)) {
        
        # update working data
        Xb = X %*% betas
        P = logitr(Xb)
        weights = as.numeric(P * (1 - P))
        z = (y - P)/weights + Xb
        
        # calculate new betas
        betas = linearr(X = X, y = z, lam = 0.1, weights = weights, 
            intercept = intercept, kernel = FALSE)$coefficients
        
        # calculate updated gradients
        grads = gradient_IRLS_logistic(betas, X, y, lam, 
            vec)
        iteration = iteration + 1
    }
    
    returns = list(coefficients = betas, total.iterations = iteration, 
        gradient = grads)
    return(returns)
}


