## Matt Galloway Augmented from Adam Rothman's
## STAT 8931 code


#' @title Linear Discriminant Analysis
#' @description this function fit the LDA model
#'
#' @param X n x p matrix where the ith row is the values of the predictor for the ith case
#' @param y n entry response vector where the ith entry is the response category in {1, ..., C} for the ith case
#' @return returns a list with the parameter estimates
#' @export
#' @examples LDA(X, y, method = 'ridge', lam = seq(0.1, 2, 0.1))

# we define the LDA function
LDA = function(X, y, method = c("MLE", "diagonal", 
    "ridge"), lam = NULL) {
    
    # which method to use?
    method = match.arg(method)
    
    # y has values in {1,...,C}
    C = max(y)
    n = length(y)
    
    # allocate memory
    pi.hats = numeric(C)
    mu.hats = list()
    Sigma.inv.hats = list()
    
    # center the X values
    X_center = NULL
    
    # loop over all response categories
    for (k in 1:C) {
        
        # indices for response category k
        indices = which(y == k)
        
        # sample size for response category k
        nk = length(indices)
        
        # pi estimate for category k
        pi.hats[k] = nk/n
        
        # center X's in category k
        X_k = X[indices, , drop = FALSE]
        mu.hats[[k]] = apply(X_k, 2, mean)
        X_center = rbind(X_center, scale(X_k, 
            center = mu.hats[[k]], scale = FALSE))
    }
    
    # if method is MLE
    if (method == "MLE") {
        
        # calculate sample covariance and precision
        S = crossprod(X_center)/n
        Sigma.inv.hat = qr.solve(S)
        
        # if method is diagonal
    } else if (method == "diagonal") {
        
        # calcualte sample covariance and precision
        S = crossprod(X_center)/n
        Sigma.inv.hat = diag(1/diag(S))
        
        # if method is ridge
    } else {
        
        # calculate sample covariance and precision
        # from sigma_ridge function
        fit = CV_sigma_ridge(X = X_center, lam = lam)
        picked.ridge = fit$best.lam
        Sigma.inv.hat = fit$omega.hat
    }
    
    for (k in 1:C) {
        
        Sigma.inv.hats[[k]] = Sigma.inv.hat
    }
    
    return(list(pi.hats = pi.hats, mu.hats = mu.hats, 
        Sigma.inv.hats = Sigma.inv.hats, picked.ridge = picked.ridge))
}



##-------------------------------------------------------------------------------



#' @title Quadratic Discriminant Analysis
#' @description this function fit the QDA model
#'
#' @param X n x p matrix where the ith row is the values of the predictor for the ith case
#' @param y n entry response vector where the ith entry is the response category in {1, ..., C} for the ith case
#' @return returns a list with the parameter estimates
#' @export
#' @examples QDA(X, y, method = 'ridge', lam = seq(0.1, 2, 0.1))

# we define the QDA function
QDA = function(X, y, method = c("MLE", "diagonal", 
    "ridge"), lam = NULL) {
    
    # which method to use?
    method = match.arg(method)
    
    
    ## y has values in {1,...,C}
    C = max(y)
    n = length(y)
    
    
    # allocate memory
    pi.hats = numeric(C)
    mu.hats = list()
    Sigma.inv.hats = list()
    picked.ridge = numeric(C)
    
    # loop over all response categories
    for (k in 1:C) {
        
        # indices for response category k indices for
        # response category k
        indices = which(y == k)
        
        # sample size for response category k
        nk = length(indices)
        
        # pi estimate for category k
        pi.hats[k] = nk/n
        
        # center X's in category k
        X_k = X[indices, , drop = FALSE]
        mu.hats[[k]] = apply(X_k, 2, mean)
        
        # if method is MLE
        if (method == "MLE") {
            
            # calculate sample covariance and precision
            S = cov(X_k) * ((nk - 1)/nk)
            Sigma.inv.hats[[k]] = qr.solve(S)
            
            # if method is diagonal
        } else if (method == "diagonal") {
            
            # calculate sample covariance and precision
            S = cov(X_k) * ((nk - 1)/nk)
            Sigma.inv.hats[[k]] = diag(1/diag(S))
            
            # if method is ridge
        } else {
            
            # calculate sample covariance and precision
            # from sigma_ridge function
            fit = CV_sigma_ridge(X = X_k, lam = lam)
            picked.ridge[k] = fit$best.lam
            Sigma.inv.hats[[k]] = fit$omega.hat
            
        }
    }
    
    
    return(list(pi.hats = pi.hats, mu.hats = mu.hats, 
        Sigma.inv.hats = Sigma.inv.hats, picked.ridge = picked.ridge))
}



##-----------------------------------------------------------------------------------



#' @title Predict QDA
#' @description this function classifies test data using a fitted QDA model
#'
#' @param fit this is a list with elements pi.hats, mu.hats, and Sigma.hats where pi.hats is a list of C response category sample proportions, mu.hats is a list of C p-dimensional sample mean proportions, Sigma.hats is a list of C p by p Sample covariance matrices
#' @param Xtest this is a matrix with ntest rows and p column, each row is a test case
#' @return returns a vector of ntest entries, where the ith entry is the estimated response category (some value in {1, ..., C}) for the ith test case.
#' @export
#' @examples predict_QDA(model, Xtest)

# we define the predict_QDA function
predict_QDA = function(fit, Xtest) {
    
    # allocate memory
    ntest = nrow(Xtest)
    C = length(fit$pi.hats)
    score.mat = matrix(NA, nrow = ntest, ncol = C)
    
    # loop over all response categories
    for (k in 1:C) {
        
        ## compute all ntest discriminant scores for
        ## category k
        tec = scale(Xtest, center = fit$mu.hats[[k]], 
            scale = FALSE)
        eout = eigen(fit$Sigma.inv.hats[[k]], 
            symmetric = TRUE, only.values = TRUE)
        ld = sum(log(eout$values))
        
        # score matrix
        score.mat[, k] = 0.5 * ld - 0.5 * diag(tec %*% 
            fit$Sigma.inv.hats[[k]] %*% t(tec)) + 
            log(fit$pi.hats[k])
        
    }
    
    # determine the best category for each of the
    # ntest cases
    pred.classes = apply(score.mat, 1, which.max)
    
    return(pred.classes)
}
