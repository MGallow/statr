## Matt Galloway


#' @title Gradient of Logistic Regression (MM)
#' @description Computes the gradient of logistic regression (optional ridge regularization term). We use this to determine if the KKT conditions are satisfied. This function is to be used with the 'MM' function.
#'
#' @param betas beta estimates (includes intercept)
#' @param X matrix or data frame
#' @param y response vector of 0,1
#' @param lam tuning parameter for ridge regularization term
#' @param alpha optional tuning parameter for bridge regularization term. Defaults to 'alpha = 1.5'
#' @param gamma indicator function. 'gamma = 1' for ridge, 'gamma = 0' for bridge. Defaults to 'gamma = 1'
#' @param vec vector to specify which coefficients will be penalized
#' @return returns the gradient
#' @examples
#' gradient_MM_logistic(betas, X, y, lam = 0.1, alpha = 1.5, penalty = 'bridge')

gradient_MM_logistic = function(betas, X, y, lam = 0, alpha = 1.5, 
    gamma = 1, vec) {
    
    # gradient for beta
    t(X) %*% (logitr(X %*% betas) - y) + lam * (gamma * 
        vec * betas + (1 - gamma) * abs(vec * betas)^(alpha - 
        1) * sign(betas))
    
}



##------------------------------------------------------------------------------------



#' @title Majorize-Minimization function
#' @description This function utilizes the MM algorithm. It will be used to compute the logistic regression coefficient estimates. This function is to be used with the 'logisticr' function.
#'
#' @param X matrix or data frame
#' @param y matrix or vector of response 0,1
#' @param lam optional tuning parameter for ridge regularization term. Defaults to 'lam = 0'
#' @param alpha optional tuning parameter for bridge regularization term. Defaults to 'alpha = 1.5'
#' @param vec optional vector to specify which coefficients will be penalized
#' @param gamma gamma indicator function. 'gamma = 1' for ridge, 'gamma = 0' for bridge. Defaults to 'gamma = 1'
#' @param intercept defaults to TRUE
#' @param tol tolerance - used to determine algorithm convergence
#' @param maxit maximum iterations
#' @return returns beta estimates (includes intercept), total iterations, and gradients.
#' @examples
#' MM(X, y)


# calculates the coefficient estimates for logistic
# regression (MM)
MM = function(X, y, lam = 0, alpha = 1.5, gamma = 1, intercept = TRUE, 
    tol = 10^(-5), maxit = 1e+05, vec = NULL) {
    
    # initialize
    n = dim(X)[1]
    p = dim(X)[2]
    X = as.matrix(X)
    y = as.matrix(y)
    delta = 10^(-5)
    betas = as.matrix(rep(0.1, p))/n
    iteration = 1
    grads = gradient_MM_logistic(betas, X, y, lam, alpha, 
        gamma, vec)
    Z = t(X) %*% X * (0.25 + delta)
    
    # MM algorithm
    while ((iteration < maxit) & (max(abs(grads)) > tol)) {
        
        # update d vector
        d = as.numeric((betas^2)^(alpha/2 - 1))
        d[1] = 0
        
        # qrsolve
        betas = qr.solve(Z + lam * diag(gamma * (vec - d) + 
            d), t(X) %*% (y - logitr(X %*% betas)) + Z %*% 
            betas)
        rownames(betas) = NULL
        
        # calculate updated gradients
        grads = gradient_MM_logistic(betas, X, y, lam, alpha, 
            gamma, vec)
        iteration = iteration + 1
        
    }
    
    # add intercept name, if needed
    rownames(betas) = NULL
    if (intercept == TRUE) {
        b1 = as.matrix(betas[1])
        rownames(b1) = "intercept"
        betas = rbind(b1, as.matrix(betas[-1, ]))
    }
    
    returns = list(coefficients = betas, total.iterations = iteration, 
        gradient = grads)
    return(returns)
}

