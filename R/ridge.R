## Matt Galloway


#' @title Ridge Regression
#' @description Computes the ridge (L2)-penalized linear regression estimates
#'
#' @param X matrix or data frame
#' @param y matrix or data frame of response values
#' @param lam tuning parameter for regularization term
#' @param K blah
#' @return returns the coefficient estimates
#' @examples
#' library(dplyr)
#' X = dplyr::select(iris, -c(Species, Sepal.Length))
#' y = dplyr::select(iris, Sepal.Length)
#' Ridge(X, y, 0.1)



ridge = function(X, y, lam, K = NULL) {

    # dimensions of data
    n = dim(X)[1]
    p = dim(X)[2]
    X = as.matrix(X)
    y = as.matrix(y)

    # if first column of ones, then remove
    if (all(X[, 1] == rep(1, n))) {
        X = X[, -1]
    }


    if (n > p) {

        # center the data
        X_bar = (1/n) * as.matrix(t(rep(1, n))) %*% X
        y_bar = mean(y)

        X_center = X - rep(1, n) %*% X_bar
        y_center = y - y_bar

        # SVD
        svd = svd(X_center %*% t(X_center))

        # adjust d vector for regularization and diagonalize
        d_adj = diag(1/(svd$d + lam))


        # calculate beta estimates
        b = t(X_center) %*% svd$v %*% d_adj %*% t(svd$u) %*% y_center
        b1 = as.numeric(y_bar - X_bar %*% b)

    }

    returns = list(Intercept = b1, Coefficients = b)
    return(returns)
}
