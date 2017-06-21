## Matt Galloway


#' @title Linear
#' @description Computes the linear regression coefficient estimates (ridge-penalization and weights, optional)
#'
#' @param X matrix or data frame
#' @param y matrix or data frame of response values
#' @param lam optional tuning parameter for ridge regularization term. If passing a list of values, the function will choose the optimal value based on K-fold cross validation. Defaults to 'lam = seq(0, 2, 0.1)'
#' @param alpha optional tuning parameter for bridge regularization term. If passing a list of values, the function will choose the optimal value based on K-fold cross validation. Defaults to 'alpha = 1.5'
#' @param penalty choose from c('none', 'ridge', 'bridge'). Defaults to 'none'
#' @param weights optional vector of weights for weighted least squares
#' @param intercept add column of ones if not already present. Defaults to TRUE
#' @param kernel use linear kernel to compute ridge regression coefficeients. Defaults to TRUE when p >> n (for 'SVD')
#' @param method optimization algorithm. Choose from 'SVD' or 'MM'. Defaults to 'SVD'
#' @param tol tolerance - used to determine algorithm convergence for 'MM'. Defaults to 10^-5
#' @param maxit maximum iterations for 'MM'. Defaults to 10^5
#' @param vec optional vector to specify which coefficients will be penalized
#' @param init optional initialization for MM algorithm
#' @param K specify number of folds for cross validation, if necessary
#' @return returns the selected tuning parameters, coefficient estimates, MSE, and gradients
#' @export
#' @examples
#'
#' Weighted ridge regression
#' library(dplyr)
#' X = dplyr::select(iris, -c(Species, Sepal.Length))
#' y = dplyr::select(iris, Sepal.Length)
#' linearr(X, y, lam = 0.1, penalty = 'ridge', weights = rep(1:150))
#'
#' Kernelized ridge regression
#' linearr(X, y, lam = 0.1, penalty = 'ridge', kernel = T)

linearr = function(X, y, lam = seq(0, 2, 0.1), alpha = 1.5, 
    penalty = "none", weights = NULL, intercept = TRUE, kernel = FALSE, 
    method = "SVD", tol = 1e-05, maxit = 1e+05, vec = NULL, 
    init = 1, K = 5) {
    
    # checks
    n = dim(X)[1]
    p = dim(X)[2]
    X = as.matrix(X)
    y = as.matrix(y)
    if (penalty == "none") {
        lam = 0
        alpha = 1.5
    }
    if (is.null(weights)) {
        weights = rep(1, n)
    }
    vec_ = vec
    if (is.null(vec)) {
        vec_ = rep(1, p)
    }
    if (length(weights) != n) 
        stop("weights must be length ", n)
    if (all(alpha >= 2 | alpha <= 1)) 
        stop("alpha must be between 1 and 2!")
    if (all(lam >= 0) == FALSE) 
        stop("lam must be nonnegative!")
    if (kernel & all(lam == 0)) 
        stop("must specify lam to use kernel!")
    if (kernel & (penalty == "bridge")) 
        stop("cannot use kernel with bridge penalty!")
    if (all(lam == 0) & (penalty != "none")) {
        print("No penalty used: lam = 0")
        penalty = "none"
    }
    
    if ((penalty == "bridge") & (method != "MM")) {
        print("using MM algorithm...")
        method = "MM"
    }
    if (penalty %in% c("none", "ridge", "bridge") == FALSE) 
        stop("incorrect penalty!")
    if (method %in% c("SVD", "MM") == FALSE) 
        stop("incorrect method!")
    if (intercept) {
        # if no first column of ones, then add it
        if (all(X[, 1] != rep(1, n))) {
            X = cbind(1, X)
            p = dim(X)[2]
        }
        # do not penalize intercept, if not specified
        if (is.null(vec)) {
            vec_ = c(0, rep(1, p - 1))
        }
    }
    if (length(init) > 1) {
        if (p != length(init)) 
            stop("initialization wrong dimension!")
    }
    
    
    # CV needed?
    if ((length(lam) > 1 | length(alpha) > 1) & (penalty != 
        "none")) {
        
        # execute CV_logisticc
        CV = CV_linearc(X, y, lam, alpha, penalty, weights, 
            intercept, kernel, method, tol, maxit, vec_, init, 
            K)
        lam = CV$best.lam
        alpha = CV$best.alpha
    }
    
    # execute linearc
    linear = linearc(X, y, lam, alpha, penalty, weights, intercept, 
        kernel, method, tol, maxit, vec_, init)
    
    
    # add intercept name, if needed
    betas = linear$coefficients
    grads = linear$gradient
    if (intercept) {
        b1 = as.matrix(betas[1])
        rownames(b1) = "intercept"
        betas = rbind(b1, as.matrix(betas[-1, ]))
        g1 = as.matrix(grads[1])
        rownames(g1) = "intercept"
        grads = rbind(g1, as.matrix(grads[-1, ]))
    }
    
    # generate fitted values
    fit = predict_linearc(linear$coefficients, as.matrix(X), 
        y)
    
    # misc
    if (penalty == "none") {
        lam = NaN
    }
    if (penalty != "bridge") {
        alpha = NaN
    }
    parameters = matrix(c(lam, alpha), ncol = 2)
    colnames(parameters) = c("lam", "alpha")
    
    returns = list(parameters = parameters, coefficients = betas, 
        MSE = fit$MSE, gradient = grads)
    return(returns)
}


