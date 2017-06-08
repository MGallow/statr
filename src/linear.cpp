// Matt Galloway

#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "MM.h"
#include "SVD.h"

using namespace Rcpp;




//' @title Linearc (c++)
//' @description Computes the linear regression coefficient estimates (ridge and bridge penalization and weights, optional)
//' @param X matrix
//' @param y matrix
//' @param lam optional tuning parameter for ridge regularization term. Defaults to 'lam = 0'
//' @param alpha optional tuning parameter for bridge regularization term. Defaults to "alpha = 1.5"
//' @param penalty choose from c("none", "ridge", "bridge"). Defaults to "none"
//' @param weights optional vector of weights for weighted least squares
//' @param intercept add column of ones if not already present. Defaults to TRUE
//' @param kernel use linear kernel to compute ridge regression coefficeients. Defaults to TRUE when p >> n (for "SVD")
//' @param method optimization algorithm. Choose from "SVD" or "MM". Defaults to "SVD"
//' @param tol tolerance - used to determine algorithm convergence for "MM". Defaults to 10^-5
//' @param maxit maximum iterations for "MM". Defaults to 10^5
//' @param vec optional vector to specify which coefficients will be penalized
//' @param init optional initialization for MM algorithm
//' @return returns the coefficient estimates
//' @export
//' @examples
//' Weighted ridge regression
//' library(dplyr)
//' X = dplyr::select(iris, -c(Species, Sepal.Length))
//' y = dplyr::select(iris, Sepal.Length)
//' linearc(X, y, lam = 0.1, penalty = "ridge", weights = rep(1:150), vec = c(0,1,1,1))
//'
//' Kernelized ridge regression
//' linearc(X, y, lam = 0.1, penalty = "ridge", kernel = T, vec = c(0,1,1,1))
//'
// [[Rcpp::export]]
List linearc(const arma::mat& X, const arma::colvec& y, double lam = 0, double alpha = 1.5, std::string penalty = "none", arma::colvec weights = 0, bool intercept = true, bool kernel = false, std::string method = "SVD", double tol = 1e-5, double maxit = 1e5, arma::colvec vec = 0, arma::colvec init = 0) {

  // checks
  int n = X.n_rows, p = X.n_cols;
  if (weights.n_rows == 1){
    weights = arma::ones<arma::colvec>(n);
  }
  if (weights.n_rows != n){
    stop("weights are wrong length!");
  }
  if (lam < 0){
    stop("lam must be nonnegative!");
  }


  // if SVD...
  arma::colvec betas = arma::ones<arma::colvec>(p);
  arma::colvec grads = arma::ones<arma::colvec>(p);
  int iterations = 1;
  if (method == "SVD"){

    // execute SVD script
    List linear = SVDc(X, y, lam, weights, intercept, kernel);

    betas = as<NumericVector>(linear["coefficients"]);
    grads = as<NumericVector>(linear["gradient"]);
  }


  // if MM algorithm...
  if (method == "MM"){

    // change gamma parameter, if needed
    int iterations = 1;
    int gamma = 1;
    if (penalty == "bridge"){
      gamma = 0;
    }

    // execute MM_linearc
    List linear = MM_linearc(X, y, lam, alpha, gamma, weights, intercept, tol, maxit, vec, init);

    iterations = as<int>(linear["total.iterations"]);
    betas = as<NumericVector>(linear["coefficients"]);
    grads = as<NumericVector>(linear["gradient"]);
    if ((iterations == maxit) & (penalty != "bridge")) {
      Rprintf("Algorithm did not converge...Try IRLS\n");
    }
  }


  // return list of coefficients
  return List::create(Named("coefficients") = betas,
                      Named("gradient") = grads);
}

