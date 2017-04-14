## Matt Galloway

# VERIFY ACCURACY - Add formulas, gradients, other outputs, combine functions? Can we put into one function if length(lam) > 1 then.... if K != NULL....


#' @title Ridge Regression
#' @description Computes the ridge (L2)-penalized linear regression estimates
#'
#' @param X matrix or data frame
#' @param y matrix or data frame of response values
#' @param lam tuning parameter for regularization term
#' @param K optional K-fold cross validation
#' @return returns the coefficient estimates
#' @export
#' @examples
#' library(dplyr)
#' X = dplyr::select(iris, -c(Species, Sepal.Length))
#' y = dplyr::select(iris, Sepal.Length)
#' Ridge(X, y, 0.1)



ridge = function(X, y, lam, K = NULL) {

  # if first column of ones, then remove
  n = dim(X)[1]
  if (all(X[, 1] == rep(1, n))) {
    X = X[, -1]
  }

  #dimensions of data
  p = dim(X)[2]
  X = as.matrix(X)
  y = as.matrix(y)

  # center the data
  X_bar = (1/n) * as.matrix(t(rep(1, n))) %*% X
  y_bar = mean(y)

  X_center = X - rep(1, n) %*% X_bar
  y_center = y - y_bar

  #if p > n, linear kernel ridge regression
  if(p > n){

    # SVD
    svd = svd(X_center %*% t(X_center))

    # adjust d vector for regularization and diagonalize
    d_adj = diag(1/(svd$d + lam))

    # calculate beta estimates
    b = t(X_center) %*% svd$v %*% d_adj %*% t(svd$u) %*% y_center
    b1 = as.numeric(y_bar - X_bar %*% b)

  }

  #else compute normal ridge regression
  else {

    # SVD
    svd = svd(X_center)

    # adjust d vector for regularization and diagonalize
    d_adj = diag(svd$d/(svd$d^2 + lam))

    # calculate beta estimates
    b = svd$v %*% d_adj %*% t(svd$u) %*% y_center
    b1 = as.numeric(y_bar - X_bar %*% b)

  }


  returns = list(Intercept = b1, Coefficients = b)
  return(returns)
}




##-----------------------------------------------------------


#' @title Ridge2
#' @description This is a special formulation of ridge regression used for cross validation. This function should not be used outside the ridge cross validation function. See: `ridge`.
#'
#' @param X_train training data (matrix or data frame)
#' @param X_test validation data (matrix or data frame)
#' @param y_train training response
#' @param y_test validation reponse
#' @param lam.vec vector of lambda values
#' @return returns RSS for each lambda value in lam.vec


#tweaked version of original ridge function
#this function returns CV error for each lambda
ridge2 = function(X_train, X_test, y_train, y_test, lam.vec){

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

  #define matrix uy := t(u) %*% y
  uy = t(svd$u) %*% y_center

  #loop over all lambda values
  CV_error = NULL
  for(lam in lam.vec){

    #adjust d vector for regularization and diagonalize
    d_adj = diag(svd$d/(svd$d^2 + lam))

    #calculate beta estimates
    b = svd$v %*% d_adj %*% uy
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

#' @title Ridge Cross Validation
#' @description Cross validation for ridge regression. This function will find the optimal lambda tuning parameter given a vector of optional values.
#'
#' @param X matrix or data frame
#' @param y matrix or data frame
#' @param lam.vec vector of lambda values
#' @param K number of folds in CV (2 or more) - defaults to 3
#' @return returns the coefficient estimates, optimal lambda, and cross validation errors
#' @export
#' @examples
#' library(dplyr)
#' X = dplyr::select(iris, -c(Species, Sepal.Length))
#' y = dplyr::select(iris, Sepal.Length)
#' CV(X, y, lam.vec, K = 5)


#ridge cross validation
#will select the optimal lambda over range of values
#returns beta estimates
ridge_CV = function(X, y, lam.vec, K = 3){

  #if first column of ones, then remove
  n = dim(X)[1]
  if (all(X[, 1] == rep(1, n))){
    X = X[, -1]
  }

  #dimensions of data
  X = as.matrix(X)
  y = as.matrix(y)
  p = dim(X)[2]
  lam = lam.vec
  CV_errors = NULL


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
    CV_errors = CV_errors + ridge2(X_train, X_test, y_train, y_test, lam.vec)

  }

  #determine optimal lambda based on CV_errors
  lam.optimal = which.min(CV_errors)

  #run ridge regression with optimal lambda using function from question 2
  lam = lam.vec[lam.optimal]

  #once we have optimal lambda, perform normal ridge
  ridge = ridge(X, y, lam)


  returns = list(b1 = ridge$b1, b = ridge$b, best.lam = lam, cv.error = CV_errors)
  return(returns)
}