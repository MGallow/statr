## Matt Galloway

# VERIFY ACCURACY


#' @title Ridge Regression
#' @description Computes the ridge (L2)-penalized linear regression estimates
#'
#' @param X matrix or data frame
#' @param y matrix or data frame of response values
#' @param lam tuning parameter for regularization term
#' @param K blah
#' @return returns the coefficient estimates
#' @export
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





##-----------------------------------------------------------


#' @title Ridge_CV
#' @description Computes the ridge (L2)-penalized linear regression estimates
#'
#' @param X_train training data (matrix or data frame)
#' @param X_test validation data (matrix or data frame)
#' @param y_train training response
#' @param y_test validation reponse
#' @param lam.vec vector of lambda values
#' @return returns the coefficient estimates
#' @export
#' @examples
#' library(dplyr)
#' X = dplyr::select(iris, -c(Species, Sepal.Length))
#' y = dplyr::select(iris, Sepal.Length)
#' Ridge_CV(X, y, 0.1)


# K-fold cross validation

#tweaked version of the function from question two
#this function returns CV error for each lambda
Ridge_CV = function(X_train, X_test, y_train, y_test, lam.vec){

  #dimensions of data
  n = dim(X_train)[1]
  p = dim(X_train)[2]

  #center the data
  X_bar = (1/n)*as.matrix(t(rep(1, n))) %*% X_train
  y_bar = mean(y_train)

  X_center = X_train - rep(1, n) %*% X_bar
  y_center = y_train - y_bar

  #SVD
  svd = svd(X_center)

  #define matrix g := t(u) %*% y
  g = t(svd$u) %*% y_center

  #loop over all lambda values
  CV_error = NULL
  for(i in lam.vec){

    #adjust d vector for regularization and diagonalize
    d_adj = diag(svd$d/(svd$d^2 + i))

    #calculate beta estimates
    b = svd$v %*% d_adj %*% g
    b1 = y_bar - X_bar %*% b

    #calculate fitted values
    fitted = X_test %*% b + as.numeric(b1)

    #calculate error
    error = sum((y_test - fitted)^2)

    #append error to vector
    CV_error = c(CV_error, error)

  }

  return(CV_error)

}



##--------------------------------------------------------------------------

#' @title CV
#' @description Cross validation for ridge regression
#'
#' @param X matrix or data frame
#' @param y matrix or data frame
#' @param lam.vec vector of lambda values
#' @param K number of folds in CV
#' @param FUN Ridge_CV object
#' @return returns the coefficient estimates
#' @export
#' @examples
#' library(dplyr)
#' X = dplyr::select(iris, -c(Species, Sepal.Length))
#' y = dplyr::select(iris, Sepal.Length)
#' CV(X, y, lam.vec, K = 5)


#this next function performs CV
#ans selects optimal lambda (over all folds)
#it will then return beta estimates
CV = function(X, y, lam.vec, K = NULL, FUN = Ridge_CV){

  #dimensions of data
  X = as.matrix(X)
  y = as.matrix(y)
  n = dim(X)[1]
  p = dim(X)[2]
  lam = lam.vec
  CV_errors = NULL

  #if first column of ones, then remove
  if (all(X[, 1] == rep(1, n))){
    X = X[, -1]
  }

  #if K = NULL (no folds) calculate ridge regr.
  if (!is.null(K)){

    #designate folds and shuffle -- ensures randomized folds
    kfold = cut(seq(1, n), breaks = K, labels = FALSE)
    kfold = kfold[sample(n)]

    #parse data into folds and perform CV
    CV_errors = rep(0, length(lam.vec))
    for(j in 1:K){

      #separate into training and testing data
      index = which(kfold == j, arr.ind = TRUE)
      X_train = X[-index,]
      X_test = X[index,]
      y_train = y[-index,]
      y_test = y[index,]

      #append CV errors
      CV_errors = CV_errors + FUN(X_train, X_test, y_train, y_test, lam.vec)

    }

    #determine optimal lambda based on CV_errors
    lam.optimal = which.min(CV_errors)

    #run ridge regression with optimal lambda using function from question 2
    lam = lam.vec[lam.optimal]
  }

  ridge = Ridge(X, y, lam)


  returns = list(b1 = ridge$b1, b = ridge$b, best.lam = lam, cv.error = CV_errors)
  return(returns)
}



## -------------------------------------------------------------------------



#' @title CV
#' @description Cross validation for ridge regression
#'
#' @param X matrix or data frame
#' @param y matrix or data frame
#' @param lam lambda 1
#' @param lam2 lambda 2
#' @param A A matrix
#' @param K number of folds in CV
#' @return returns the coefficient estimates
#' @export
#' @examples
#' A = matrix(c(1, -1, 0, 0, -1, 2, -1, 0, 0, -1, 2, -1, 0, 0, -1, 1), ncol = 4)
#' New_Ridge(X, y, lam, lam2, A)

# altered ridge regression estimators

#calculates the beta estimates for the altered ridge regression
#this accounts for the second penalization term
New_Ridge = function(X, y, lam, lam2, A, K = NULL){

  #dimensions of data
  n = dim(X)[1]
  p = dim(X)[2]
  X = as.matrix(X)
  y = as.matrix(y)

  #if first column of ones, then remove
  if (all(X[, 1] == rep(1, n))){
    X = X[, -1]
  }

  #center the data
  X_bar = (1/n)*as.matrix(t(rep(1, n))) %*% X
  y_bar = mean(y)

  X_center = X - rep(1, n) %*% X_bar
  y_center = y - y_bar
  new_lam = lam*diag(p) + lam2*diag(p) %*% A

  #we will not use SVD this time for simplicity
  b = ginv(t(X_center) %*% X_center + new_lam) %*% t(X_center) %*% y_center
  b1 = y_bar - X_bar %*% b

  returns = list(b1 = b1, b = b)
  return(returns)
}