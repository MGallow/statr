## Matt Galloway

# VERIFY ACCURACY



#' @title Softmax
#' @description Computes the probability for an individual observation (x) being in class "1".
#'
#' @param x data vector
#' @param b beta estimates
#' @return returns the estimated probability
#' @examples
#' softmax(x, b)


softmax = function(x, b){

  #calculate softmax probabilities
  exp(x %*% b)/(1 + exp(x %*% b))

}



##------------------------------------------------------------------------------------


#' @title Logistic Bridge Gradient
#' @description Computes the gradient of bridge-penalized logistic regression. Use to determine if KKT conditions are satisfied
#'
#' @param X matrix or data frame
#' @param y matrix or data frame
#' @param N diag(n.list) - where n.list is vector of the number of replicates for each observation
#' @param lam tuning parameter for regularization term
#' @param alpha tuning parameter. alpha in (1, 2).
#' @return returns the gradient.
#' @examples
#' logistic_bridge_gradient(X, y, N = diag(n.list), b = betas, lam = 0.1, alpha = 1.5)


logistic_bridge_gradient = function(X, y, N, b, lam, alpha){

  #gradient for beta
  b_reg = b
  b_reg[1] = 0
  grad_b = t(X) %*% N %*% (softmax(X, b) - y) + lam*(abs(b_reg)^(alpha - 1))*sign(b_reg)

  returns = list(grad = c(grad_b))
  return(returns)

}



##--------------------------------------------------------------------------------------------


#' @title Bridge-Penalized Logistic Regression
#' @description Computes the bridge-penalized logistic regression estimates. We assume alpha in (1, 2). Uses the MM algorithm to compute.
#'
#' @param X matrix or data frame
#' @param y matrix or data frame of response values
#' @param n.list vector of replicates for each observation
#' @param lam tuning parameter for regularization term
#' @param alpha tuning parameter. alpha in (1, 2).
#' @param tol tolerance - used to determine algorithm convergence
#' @param maxit maximum iterations
#' @return returns beta estimates, total iterations, and gradients.
#' @export
#' @examples
#' bridge_logistic(X, y, n.list = c(rep(1, n)), lam = 0.1, alpha = 1.5)


# Bridge penalized logistic regression estimators

#calculates the beta estimates for bridge logistic regression
bridge_logistic = function(X, y, n.list, lam, alpha, tol = 10^(-5), maxit = 100000){

  #dimensions of data
  n = dim(X)[1]
  p = dim(X)[2]
  X = as.matrix(X)
  y = as.matrix(y)

  #initialize
  iteration = 1
  delta = 10^(-5)
  b = as.matrix(rep(1, p))
  N = diag(n.list)
  gamma = t(X) %*% ((0.25 + delta)*N) %*% X
  grads = logistic_bridge_gradient(X, y, N, b, lam, alpha)


  #MM algorithm
  while ((iteration < maxit) & (max(grads$grad) > tol)){

    #adjust v vector
    v = as.numeric((b^2)^(alpha/2 - 1))
    v[1] = 0

    #qrsolve
    b = qr.solve(gamma + lam*diag(v), t(X) %*% N %*% (y - softmax(X, b)) + gamma %*% b)

    #calculate updated gradients
    grads = logistic_bridge_gradient(X, y, N, b, lam, alpha)
    iteration = iteration + 1

  }


  returns = list(b = b, total.iterations = iteration, grads = grads)
  return(returns)
}


