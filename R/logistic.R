## Matt Galloway



#' @title Logistic Regression
#' @description Computes the coefficient estimates for logistic regression. ridge regularization and bridge regularization optional.
#'
#' @param X matrix or data frame
#' @param y matrix or vector of response values 0,1
#' @param lam optional tuning parameter(s) for ridge regularization term. If passing a list of values, the function will choose optimal value based on K-fold cross validation. Defaults to `lam = seq(0, 2, 0.1)`
#' @param alpha optional tuning parameter for bridge regularization term. If passing a list of values, the function will choose the optimal value based on K-fold cross validation. Defaults to 'alpha = 1.5'
#' @param penalty choose from c('none', 'ridge', 'bridge'). Defaults to 'none'
#' @param intercept Defaults to TRUE
#' @param method optimization algorithm. Choose from 'IRLS' or 'MM'. Defaults to 'IRLS'
#' @param tol tolerance - used to determine algorithm convergence. Defaults to 10^-5
#' @param maxit maximum iterations. Defaults to 10^5
#' @param vec optional vector to specify which coefficients will be penalized
#' @param init optional initialization for MM algorithm
#' @param criteria specify the criteria for cross validation. Choose from c('mse', 'logloss', 'misclass'). Defauls to 'logloss'
#' @param K specify number of folds for cross validation, if necessary
#' @return returns selected tuning parameters, beta estimates (includes intercept), MSE, log loss, misclassification rate, total iterations, and gradients.
#' @export
#' @examples
#' Logistic Regression
#' library(dplyr)
#' X = dplyr::select(iris, -Species)
#' y = dplyr::select(iris, Species)
#' y$Species = ifelse(y$Species == 'setosa', 1, 0)
#' logisticr(X, y)
#'
#' ridge Logistic Regression with IRLS
#' logistir(X, y, lam = 0.1, penalty = 'ridge')
#'
#' ridge Logistic Regression with MM
#' logisticr(X, y, lam = 0.1, penalty = 'ridge', method = 'MM')
#'
#' bridge Logistic Regression
#' (Defaults to MM -- IRLS will return error)
#' logisticr(X, y, lam = 0.1, alpha = 1.5, penalty = 'bridge')


logisticr = function(X, y, lam = seq(0, 2, 0.1), alpha = 1.5, 
    penalty = "none", intercept = TRUE, method = "IRLS", tol = 1e-05, 
    maxit = 1e+05, vec = NULL, init = 1, criteria = "logloss", 
    K = 5) {
    
    # checks
    n = dim(X)[1]
    p = dim(X)[2]
    X = as.matrix(X)
    y = as.matrix(y)
    if (penalty == "none") {
        lam = 0
        alpha = 1.5
    }
    if (all(y == 1 | y == 0) == FALSE) 
        stop("y must be binary!")
    vec_ = vec
    if (is.null(vec)) {
        vec_ = rep(1, p)
    }
    if (all(alpha >= 2 | alpha <= 1)) 
        stop("alpha must be between 1 and 2!")
    if (all(lam >= 0) == FALSE) 
        stop("lam must be nonnegative!")
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
    if (criteria %in% c("mse", "logloss", "misclass") == FALSE) 
        stop("incorrect criteria!")
    if (method %in% c("IRLS", "MM") == FALSE) 
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
        CV = CV_logisticc(X, y, lam, alpha, penalty, intercept, 
            method, tol, maxit, vec_, init, criteria, K)
        lam = CV$best.lam
        alpha = CV$best.alpha
    }
    
    # execute logisticc
    logistic = logisticc(X, y, lam, alpha, penalty, intercept, 
        method, tol, maxit, vec_, init)
    
    
    # add intercept name, if needed
    betas = logistic$coefficients
    grads = logistic$gradient
    if (intercept == TRUE) {
        b1 = as.matrix(betas[1])
        rownames(b1) = "intercept"
        betas = rbind(b1, as.matrix(betas[-1, ]))
        g1 = as.matrix(grads[1])
        rownames(g1) = "intercept"
        grads = rbind(g1, as.matrix(grads[-1, ]))
    }
    
    # generate fitted values
    fit = predict_logisticc(logistic$coefficients, as.matrix(X), 
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
        MSE = fit$MSE, log.loss = fit$log.loss, misclassification = fit$misclassification, 
        total.iterations = logistic$total.iterations, gradient = grads)
    return(returns)
    
}
