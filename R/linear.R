## Matt Galloway


##------------------------------------------------------------------------------------


#' @title Gradient of Linear Regression
#' @description Computes the gradient of linear regression (optional ridge regularization term). This function is to be used with the 'Linearr' function.
#'
#' @param betas beta estimates (includes intercept)
#' @param X matrix or data frame
#' @param y response vector of 0,1
#' @param lam tuning parameter for ridge regularization term
#' @param weights option vector of weights for weighted least squares
#' @param vec vector to specify which coefficients will be penalized
#' @return returns the gradient
#' @examples
#'
#' gradient_linear(betas, X, y, lam = 0.1)

gradient_linear = function(betas, X, y, lam = 0, weights = NULL, 
    vec) {
    
    # add intercept, if needed
    if (ncol(X) != nrow(betas)) {
        X = cbind(1, X)
    }
    
    # gradient for beta
    -t(X) %*% diag(weights) %*% y + t(X) %*% diag(weights) %*% 
        X %*% betas + lam * vec * betas
    
}



##--------------------------------------------------------------------------------------------



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
    p = dim(X)[2]
    if (is.null(weights)) {
        weights = rep(1, n)
    }
    if (length(weights) != n) 
        stop("weights must be length ", n)
    if (length(lam) > 1) 
        stop("lam must be a scalar!")
    if (lam < 0) 
        stop("lam must be nonnegative!")
    if ((kernel == TRUE) & (lam == 0)) 
        stop("must specify lam to use kernel!")
    
    
    # initialization
    X. = as.matrix(X)
    y. = as.matrix(y)
    if (intercept == TRUE) {
        # if first column of ones, then remove it
        if (all(X.[, 1] == rep(1, n))) {
            X. = X.[, -1]
        }
        
        # center the data
        X_bar = (as.matrix(t(weights)) %*% X.)/sum(weights)
        y_bar = sum(weights * y.)/sum(weights)
        
        X. = X. - rep(1, n) %*% X_bar
        y. = y. - y_bar
        
    }
    
    root = sqrt(weights)
    X. = root * X.
    y. = root * y.
    
    
    # if p > n, linear kernel ridge regression
    if ((p > n) | (kernel == TRUE)) {
        
        # SVD
        svd = svd(X.)
        
        # adjust d vector for regularization and diagonalize
        d_adj = diag(1/(svd$d^2 + lam))
        
        # calculate beta estimates
        betas = t(X.) %*% svd$u %*% d_adj %*% t(svd$u) %*% 
            y.
        rownames(betas) = NULL
        
    } else {
        # compute normal ridge regression
        
        # SVD
        svd = svd(X.)
        
        # adjust d vector for regularization and diagonalize
        d_adj = diag(svd$d/(svd$d^2 + lam))
        
        # calculate beta estimates
        betas = svd$v %*% d_adj %*% t(svd$u) %*% y.
        rownames(betas) = NULL
        
    }
    
    # add intercept, if needed
    if (intercept == TRUE) {
        
        # calculate intercept
        b1 = y_bar - X_bar %*% betas
        rownames(b1) = c("intercept")
        betas = rbind(b1, betas)
        
    }
    
    returns = list(coefficients = betas)
    return(returns)
}


