//// Matt Galloway

#include <RcppArmadillo.h>
#include <Rcpp.h>
// #include "linear.cpp"
// #include "IRLS.cpp"

using namespace Rcpp;


//' @title Gradient of Logistic Regression (MM) (c++)
//' @description Computes the gradient of logistic regression (optional ridge regularization term). We use this to determine if the KKT conditions are satisfied. This function is to be used with the 'MM' function.
//'
//' @param betas beta estimates (includes intercept)
//' @param X matrix or data frame
//' @param y response vector of 0,1
//' @param lam tuning parameter for ridge regularization term
//' @param alpha optional tuning parameter for bridge regularization term. Defaults to 'alpha = 1.5'
//' @param gamma indicator function. 'gamma = 1' for ridge, 'gamma = 0' for bridge. Defaults to 'gamma = 1'
//' @param vec vector to specify which coefficients will be penalized
//' @return returns the gradient
//' @examples
//' gradient_MM_logistic(betas, X, y, lam = 0.1, alpha = 1.5, penalty = 'bridge')
//'
// [[Rcpp::export]]
arma::colvec gradient_MM_logisticc(const arma::colvec& betas, const arma::mat& X, const arma::colvec& y, double lam = 0, double alpha = 1.5, double gamma = 1, const arma::colvec& vec = 0) {

  // gradient for beta
  arma::colvec logitc(const arma::colvec& u);
  return arma::trans(X) * (logitc(X * betas) - y) + lam * (gamma * (vec % betas) + (1 - gamma) * arma::pow(arma::abs(vec % betas), alpha - 1) % arma::sign(betas));

}



//------------------------------------------------------------------------------------



//' @title Majorize-Minimization function (c++)
//' @description This function utilizes the MM algorithm. It will be used to compute the logistic regression coefficient estimates. This function is to be used with the 'logisticr' function.
//'
//' @param X matrix or data frame
//' @param y matrix or vector of response 0,1
//' @param lam optional tuning parameter for ridge regularization term. Defaults to 'lam = 0'
//' @param alpha optional tuning parameter for bridge regularization term. Defaults to 'alpha = 1.5'
//' @param vec optional vector to specify which coefficients will be penalized
//' @param gamma gamma indicator function. 'gamma = 1' for ridge, 'gamma = 0' for bridge. Defaults to 'gamma = 1'
//' @param intercept defaults to TRUE
//' @param tol tolerance - used to determine algorithm convergence
//' @param maxit maximum iterations
//' @return returns beta estimates (includes intercept), total iterations, and gradients.
//' @examples
//' MMc(X, y)
//'
// [[Rcpp::export]]
List MMc(const arma::mat& X, const arma::colvec& y, double lam = 0, double alpha = 1.5, double gamma = 1, bool intercept = true, double tol = 1e-5, double maxit = 1e5, const arma::colvec& vec = 0) {

  // initialize
  int n = X.n_rows, p = X.n_cols;
  double delta = 1e-5;
  arma::colvec betas = 0.1*arma::ones<arma::colvec>(p)/n;
  int iteration = 1;
  arma::colvec grads = gradient_MM_logisticc(betas, X, y, lam, alpha, gamma, vec);
  arma::mat Z = arma::trans(X) * X * (0.25 + delta);

  // MM algorithm
  while ((iteration < maxit) & (arma::max(arma::abs(grads)) > tol)) {

    // update d vector
    arma::colvec d = arma::pow((betas % betas), alpha/2 - 1);
    if (intercept) {
      d[0] = 0;
    }

    // qrsolve
    arma::colvec logitc(const arma::colvec& u);
    betas = arma::solve(Z + lam * arma::diagmat(gamma * (vec - d) + d), arma::trans(X) * (y - logitc(X * betas)) + Z * betas);

    // calculate updated gradients
    grads = gradient_MM_logisticc(betas, X, y, lam, alpha, gamma, vec);
    iteration += 1;

  }


  return List::create(Named("coefficients") = betas,
                      Named("total.iterations") = iteration,
                      Named("gradient") = grads);
}

