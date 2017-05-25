// Matt Galloway

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

//' @title Linearc (c++)
//' @description Computes the linear regression coefficient estimates (ridge-penalization and weights, optional)
//' @param X matrix
//' @param y matrix
//' @param lam optional tuning parameter for ridge regularization term. Defaults to 'lam = 0'
//' @param weights optional vector of weights for weighted least squares
//' @param intercept add column of ones if not already present. Defaults to TRUE
//' @param kernel use linear kernel to compute ridge regression coefficeients. Defaults to true when p >> n
//' @return returns the coefficient estimates
//' @export
//' @examples
//' Weighted ridge regression
//' library(dplyr)
//' X = dplyr::select(iris, -c(Species, Sepal.Length))
//' y = dplyr::select(iris, Sepal.Length)
//' linearc(X, y, lam = 0.1, weights = rep(1:150))
//'
//' Kernelized ridge regression
//' linearc(X, y, lam = 0.1, kernel = T)
//'
// [[Rcpp::export]]
List linearc(const arma::mat& X, const arma::colvec& y, double lam = 0, arma::colvec weights = 0, bool intercept = true, bool kernel = false) {

  // checks
  int n = X.n_rows, p = X.n_cols;
  if (weights.n_rows == 1){
    weights = arma::ones<arma::colvec>(n);
  }
  if (weights.n_rows != n){
    Rcout << "weights are wrong length!\n";
    return 0;
  }
  if (lam < 0){
    Rcout << "lam must be nonnegative!\n";
    return 0;
  }
  if (kernel & (lam == 0)){
    Rcout << "must specify lam to use kernel!\n";
    return 0;
  }


  // initialization
  arma::mat X_ = X;
  arma::colvec y_ = y;
  arma::rowvec X_bar = arma::ones<arma::rowvec>(p);
  double y_bar = 0;
  if (intercept){
    // if first column of ones, then remove it
    if (arma::all(X.col(0) == 1)){
      X_ = X.submat(0, 1, n - 1, p - 1);
      p = X_.n_cols;
    }

    // center the data
    X_bar = arma::trans(weights)*X_/arma::sum(weights);
    y_bar = arma::sum(weights % y_)/arma::sum(weights);

    X_ = X_ - arma::ones<arma::colvec>(n)*X_bar;
    y_ = y_ - y_bar;
  }

  // incorporate weights into "new data"
  arma::mat root = arma::sqrt(weights);
  X_ = X_.each_col() % root;
  y_ = y_ % root;

  // run SVD
  arma::mat U;
  arma::colvec D;
  arma::mat V;
  svd_econ(U, D, V, X_);


  // if p > n, linear kernel ridge regression
  arma::colvec betas = arma::ones<arma::colvec>(p);
  if ((p > n) & kernel){

    // adjust d vector for regularization
    arma::colvec d = arma::zeros<arma::colvec>(p);
    for (int i = 0; i < p; i++){
      d[i] = 1/(D[i]*D[i] + lam);
    }

    // calculate beta estimates
    betas = arma::trans(X_)*U*arma::diagmat(d)*arma::trans(U)*y_;

  } else {
    // compute normal ridge regression

    // adjust d vector for regularization
    arma::colvec d = arma::zeros<arma::colvec>(p);
    for (int i = 0; i < p; i++){
      d[i] = D[i]/(D[i]*D[i] + lam);
    }

    // calculate beta estimates
    betas = V*arma::diagmat(d)*arma::trans(U)*y_;

  }

  // add intercept, if needed
  if (intercept){

    // calculate intercept
    arma::mat b1 = y_bar - X_bar*betas;
    betas.insert_rows(0, b1);

  }


  // return list of coefficients
  return List::create(Named("coefficients") = betas);
}

