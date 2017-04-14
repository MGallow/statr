## Matt Galloway

# VERIFY ACCURACY - add formulas



##-------------------------------------------------

#' @title Gradient of Bridge-Penalized Regression
#' @description Computes the gradient of bridge regression. We use this to determine if the KKT conditions are satisfied.
#'
#' @param XX crossprod(X, X) - X is centered, no intercept
#' @param Xy crossprod(X, y) - y is centered
#' @param b beta estimates (intercept not included)
#' @param lam tuning parameter for regularization term
#' @param alpha tuning parameter. Alpha is contained in (1, 2).
#' @return returns the gradient.
#' @examples
#' bridge_gradient(XX, Xy, b = betas, lam = 0.1, alpha = 1.5)


bridge_gradient = function(XX, Xy, b, lam, alpha) {
    
    # gradient for beta
    grad_b = -Xy + XX %*% b + lam * (abs(b)^(alpha - 1)) * 
        sign(b)
    
    returns = list(grad = c(grad_b))
    return(returns)
    
}




#' @title Bridge Regression
#' @description Computes the bridge-penalized linear regression estimates. We assume alpha in (1, 2). Uses the MM algorithm to compute.
#'
#' @param X matrix or data frame
#' @param y matrix or data frame of response values
#' @param lam tuning parameter for regularization term
#' @param alpha tuning parameter. alpha in (1, 2).
#' @param tol tolerance - used to determine algorithm convergence
#' @param maxit maximum iterations
#' @return returns beta estimates, total iterations, and gradients.
#' @export
#' @examples
#' library(dplyr)
#' X = dplyr::select(iris, -c(Species, Sepal.Length))
#' y = dplyr::select(iris, Sepal.Length)
#' Ridge(X, y, lam = 0.1, alpha = 1.5)


# Bridge Regression
bridge = function(X, y, lam, alpha, tol = 10^(-5), maxit = 1000) {
    
    # if first column of ones, then remove
    n = dim(X)[1]
    if (all(X[, 1] == rep(1, n))) {
        X = X[, -1]
    }
    
    # dimensions of data
    X = as.matrix(X)
    y = as.matrix(y)
    p = dim(X)[2]
    
    # center the data
    X_bar = (1/n) * as.matrix(t(rep(1, n))) %*% X
    y_bar = mean(y)
    
    X_center = X - rep(1, n) %*% X_bar
    y_center = y - y_bar
    
    # initialize
    iteration = 1
    b = as.matrix(rep(1, p))
    
    # calculate outside loop (saves computation time)
    XX = t(X_center) %*% X_center
    Xy = t(X_center) %*% y_center
    grads = bridge_gradient(XX, Xy, b, lam, alpha)
    
    
    # MM algorithm
    while ((iteration < maxit) & (max(grads$grad) > tol)) {
        
        # adjust v vector
        v = as.numeric((b^2)^(alpha/2 - 1))
        
        # calculate new beta estimate with qr.solve
        b = qr.solve(XX + lam * diag(v), Xy)
        
        # calculate updated gradients
        grads = bridge_gradient(XX, Xy, b, lam, alpha)
        iteration = iteration + 1
        
    }
    
    
    # calculate intercept
    b1 = y_bar - X_bar %*% b
    
    
    returns = list(b1 = b1, b = b, total.iterations = iteration, 
        grads = grads)
    return(returns)
}

