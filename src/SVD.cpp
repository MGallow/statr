// Matt Galloway

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;


//' @title Gradient of Linear Regression (c++)
//' @description Computes the gradient of linear regression (optional ridge regularization term). This function is to be used with the 'SVDc' function.
//'
//' @param beta estimates (includes intercept)
//' @param X matrix
//' @param y response vector of 0,1
//' @param lam tuning parameter for ridge regularization term
//' @param weights option vector of weights for weighted least squares
//' @param intercept add column of ones if not already present. Defaults to TRUE
//' @return returns the gradient
//' @examples
//' gradient_linearc(betas, X, y, lam = 0.1, weights = rep(1,150), intercept = TRUE)
//'
// [[Rcpp::export]]
arma::colvec gradient_linearc(const arma::colvec& betas, const arma::mat& X, const arma::colvec& y, double lam = 0, const arma::colvec& weights = 0, bool intercept = true) {

  // do not penalize intercept
  arma::mat X_ = X;
  int n = X_.n_rows, p = X_.n_cols;
  arma::colvec vec = arma::ones<arma::colvec>(p);

  // add intercept, if needed
  if (intercept){
    if (p != betas.n_rows){
      arma::colvec X1 = arma::ones<arma::colvec>(n);
      X_.insert_cols(0, X1);
      p = X_.n_cols;
    }
    vec = arma::ones<arma::colvec>(p - 1);
    vec.insert_rows(0, arma::zeros<arma::colvec>(1));
  }

  // gradient for beta
  return -arma::trans(X_) * arma::diagmat(weights) * y + arma::trans(X_) * arma::diagmat(weights) * X_ * betas + lam * vec % betas;

}



//--------------------------------------------------------------------------------------------


//' @title Linear Singular Value Decomposition (c++)
//' @description Computes the logistic regression coefficient estimates using SVD. This function is to be used with the 'linearc' function.
//'
//' @param X matrix
//' @param y matrix
//' @param lam optional tuning parameter for ridge regularization term. Defaults to 'lam = 0'
//' @param weights optional vector of weights for weighted least squares
//' @param intercept add column of ones if not already present. Defaults to TRUE
//' @param kernel use linear kernel to compute ridge regression coefficeients. Defaults to TRUE when p >> n (for "SVD")
//' @return returns beta estimates (includes intercept) and gradients.
//' @examples
//' SVDc(X, y, lam = 0.1 weights = rep(1, 150))
//'
// [[Rcpp::export]]
List SVDc(const arma::mat& X, const arma::colvec& y, double lam = 0, arma::colvec weights = 0, bool intercept = true, bool kernel = false) {

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


  // initialization
  arma::mat X_ = X;
  arma::colvec y_ = y;
  arma::colvec betas = arma::ones<arma::colvec>(p);
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
  arma::svd_econ(U, D, V, X_);


  // if p > n, linear kernel ridge regression
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

  // calculate gradients
  arma::colvec grads = arma::ones<arma::colvec>(p);
  grads = gradient_linearc(betas, X, y, lam, weights, intercept);


  return List::create(Named("coefficients") = betas,
                      Named("gradient") = grads);
}


