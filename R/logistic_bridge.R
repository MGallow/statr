## Matt Galloway

# VERIFY ACCURACY - add formulas (Q function) to bridge
# logistic



#' @title Softmax
#' @description Computes the probability for an individual observation (x) being in class '1'.
#'
#' @param x data vector
#' @param b beta estimates (includes intercept)
#' @return returns the estimated probability
#' @examples
#' softmax(x, b)


softmax = function(x, b) {
    
    # calculate softmax probabilities
    exp(x %*% b)/(1 + exp(x %*% b))
    
}



##------------------------------------------------------------------------------------


#' @title Gradient of Bridge-Penalized Logistic Regression
#' @description Computes the gradient of bridge-penalized logistic regression. We use this to determine if the KKT conditions are satisfied.
#'
#' @param X matrix or data frame
#' @param y response vector of 0,1
#' @param N diag(n.list) - where n.list is a vector of the number of replicates for each observation
#' @param b beta estimates (includes intercept)
#' @param lam tuning parameter for regularization term
#' @param alpha tuning parameter. Alpha is contained in (1, 2).
#' @return returns the gradient.
#' @examples
#' logistic_bridge_gradient(X, y, N = diag(n.list), b = betas, lam = 0.1, alpha = 1.5)


gradient_logistic_bridge = function(X, y, N, b, lam, alpha) {
    
    # if first column not vector of ones, add it
    n = dim(X)[1]
    if (all(X[, 1] != rep(1, n))) {
        X = cbind(1, X)
    }
    
    # do not penalize the intercept set first element equal
    # to zero
    b_reg = b
    b_reg[1] = 0
    
    # gradient for beta
    grad_b = t(X) %*% N %*% (softmax(X, b) - y) + lam * 
        (abs(b_reg)^(alpha - 1)) * sign(b_reg)
    
    returns = list(grad = c(grad_b))
    return(returns)
    
}



##--------------------------------------------------------------------------------------------


#' @title Bridge-Penalized Logistic Regression
#' @description Computes the bridge-penalized logistic regression estimates. We assume alpha is contained in (1, 2). Uses the MM algorithm to compute.
#'
#' @param X Matrix or data frame
#' @param y response vector of 0,1
#' @param n.list vector containing the number of replicates for each observation
#' @param lam tuning parameter for regularization term
#' @param alpha tuning parameter. Alpha is contained in (1, 2).
#' @param tol tolerance - used to determine algorithm convergence
#' @param maxit maximum iterations
#' @return returns beta estimates (includes intercept), total iterations, and gradients.
#' @export
#' @examples
#' bridge_logistic(X, y, n.list = c(rep(1, n)), lam = 0.1, alpha = 1.5)


# Bridge penalized logistic regression estimators

# calculates the beta estimates for bridge logistic
# regression
bridge_logistic = function(X, y, n.list, lam, alpha, tol = 10^(-5), 
    maxit = 1e+05) {
    
    # if first column not vector of ones, add it
    n = dim(X)[1]
    if (all(X[, 1] != rep(1, n))) {
        X = cbind(1, X)
    }
    
    # dimensions of data
    p = dim(X)[2]
    X = as.matrix(X)
    y = as.matrix(y)
    
    # initialize
    iteration = 1
    delta = 10^(-5)
    b = as.matrix(rep(1, p))
    N = diag(n.list)
    grads = gradient_logistic_bridge(X, y, N, b, lam, alpha)
    
    # part of Q maximization function
    gamma = t(X) %*% ((0.25 + delta) * N) %*% X
    
    
    # MM algorithm
    while ((iteration < maxit) & (max(grads$grad) > tol)) {
        
        # adjust v vector
        v = as.numeric((b^2)^(alpha/2 - 1))
        v[1] = 0
        
        # qrsolve
        b = qr.solve(gamma + lam * diag(v), t(X) %*% N %*% 
            (y - softmax(X, b)) + gamma %*% b)
        
        # calculate updated gradients
        grads = gradient_logistic_bridge(X, y, N, b, lam, 
            alpha)
        iteration = iteration + 1
        
    }
    
    
    returns = list(b = b, total.iterations = iteration, 
        grads = grads)
    return(returns)
}


