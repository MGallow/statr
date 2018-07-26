## Matt Galloway


#' @title Ridge regression
#' @description calculate ridge regression coefficients using the optimal tuning parameter from the glmnet package.
#'
#' @param X nxp data matrix. Each row corresponds to a single observation and each column contains n observations of a single feature/variable.
#' @param Y nxr response matrix. Each row corresponds to a single response and each column contains n response of a single feature/response.
#' @param lam tuning parameter
#' @param intercept option to include intercept, defaults to FALSE
#' @param standardize option to standardize the data, defaults to FALSE
#' @param ... other options to pass to glmnet
#' @return betas, lam
#' @export

# we define the ridge regression function
RIDGE = function(X, Y, lam = NULL, intercept = FALSE, standardize = FALSE, 
    ...) {
    
    # which family?
    if (ncol(Y) > 1) {
        family = "mgaussian"
    } else {
        family = "gaussian"
    }
    
    # is lam specified?
    if (is.null(lam)) {
        lam = glmnet::cv.glmnet(x = X, y = Y, alpha = 0, 
            intercept = intercept, family = family, standardize = standardize, 
            ...)$lambda.min
    }
    
    # calculate coefficients
    betas = glmnet::glmnet(x = X, y = Y, standardize = standardize, 
        intercept = intercept, family = family, alpha = 0, 
        lambda = lam, ...)
    
    if (ncol(Y) > 1) {
        betas = do.call(cbind, coef(betas))
    } else {
        betas = coef(betas)
    }
    if (intercept) {
        betas = as.matrix(betas)
    } else {
        betas = as.matrix(betas[-1, ])
    }
    
    returns = list(betas = betas, lam = lam)
    return(returns)
    
}


##-----------------------------------------------------




#' @title Lasso regression
#' @description calculate lasso regression coefficients using the optimal tuning parameter from the glmnet package.
#'
#' @param X nxp data matrix. Each row corresponds to a single observation and each column contains n observations of a single feature/variable.
#' @param Y nxr response matrix. Each row corresponds to a single response and each column contains n response of a single feature/response.
#' @param lam tuning parameter
#' @param intercept option to include intercept, defaults to FALSE
#' @param standardize option to standardize the data, defaults to FALSE
#' @param ... other options to pass to glmnet
#' @return betas, lam
#' @export

# we define the ridge regression function
LASSO = function(X, Y, lam = NULL, intercept = FALSE, standardize = FALSE, 
    ...) {
    
    # which family?
    if (ncol(Y) > 1) {
        family = "mgaussian"
    } else {
        family = "gaussian"
    }
    
    # is lam specified?
    if (is.null(lam)) {
        lam = glmnet::cv.glmnet(x = X, y = Y, alpha = 1, 
            intercept = intercept, family = family, standardize = standardize, 
            ...)$lambda.min
    }
    
    # calculate coefficients
    betas = glmnet::glmnet(x = X, y = Y, standardize = standardize, 
        intercept = intercept, family = family, alpha = 1, 
        lambda = lam, ...)
    
    if (ncol(Y) > 1) {
        betas = do.call(cbind, coef(betas))
    } else {
        betas = coef(betas)
    }
    if (intercept) {
        betas = as.matrix(betas)
    } else {
        betas = as.matrix(betas[-1, ])
    }
    
    returns = list(betas = betas, lam = lam)
    return(returns)
    
}


##-----------------------------------------------------



#' @title Mean Frobenius norm
#' @description calculates the average squared value of an object. That is, all elements are squared, summed, and divided by the total number of elements
#'
#' @param X object
#' @return norm
#' @export

# we define the frobenius norm function
fro = function(X) {
    
    mean(X^2)
    
}


##-----------------------------------------------------



#' @title CV split
#' @description splits data objects into training and testing sets
#'
#' @param X nxp data matrix. Each row corresponds to a single observation and each column contains n observations of a single feature/variable.
#' @param Y nxr response matrix. Each row corresponds to a single response and each column contains n response of a single feature/response.
#' @param split fraction of objects devoted to training set
#' @param N option to provide number of objects devoted to training set
#' @return X.train, Y.train, X.test, Y.test
#' @export

# we define the CVsplit function
CVsplit = function(X, Y, split = 0.5, N = NULL) {
    
    # checks
    if ((split <= 0) || (split >= 1)) {
        stop("split must be contained in c(0, 1)!")
    }
    
    # specify leave.out
    if (!is.null(N)) {
        if (N > nrow(X)) {
            stop("N cannot exceed number of observations in X!")
        }
        leave.out = sample(nrow(X), N)
    } else {
        leave.out = sample(nrow(X), floor(nrow(X) * split))
    }
    
    # training sets
    X.train = X[leave.out, , drop = FALSE]
    Y.train = Y[leave.out, , drop = FALSE]
    
    # testing sets
    X.test = X[-leave.out, , drop = FALSE]
    Y.test = Y[-leave.out, , drop = FALSE]
    
    returns = list(X.train = X.train, Y.train = Y.train, 
        X.test = X.test, Y.test = Y.test)
    return(returns)
    
}


##-----------------------------------------------------
