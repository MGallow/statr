//// Matt Galloway

#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "IRLS.h"

using namespace Rcpp;


//' @title Predict Logistic Regression (c++)
//' @description Generates prediction for logistic regression
//'
//' @param betas matrix of coefficientts
//' @param X matrix of (new) observations
//' @param y matrix of response values 0,1
//' @return predictions and loss metrics
//' @export
//' @examples
//'
//' fitted = logisticr(X, y, lam = 0.1, penalty = 'ridge', method = 'MM')
//' predict_logisticr(fitted$coefficients, X)
//'
// [[Rcpp::export]]
List predict_logisticc(const arma::colvec& betas, const arma::mat& X, const arma::colvec& y = 0) {

    // checks
    if (y.n_rows > 1) {
        if (X.n_rows != y.n_rows)
            stop("X and y must have equal observations!");
    }

    // add intercept, if needed
    arma::mat X_ = X;
    int n = X_.n_rows, p = X_.n_cols;
    if (p != betas.n_rows) {
        arma::colvec X1 = arma::ones<arma::colvec>(n);
      X_.insert_cols(0, X1);
      p = X_.n_cols;
    }

    // fitted values
    arma::colvec fitted = logitc(X_ * betas);
    arma::colvec classes = arma::round(fitted);

    // if y, return MSE, misclassification
    double MSE = 0;
    double log_loss = 0;
    double misclassification = 0;
    if (y.n_rows > 1) {

        // calculate metrics
        MSE = arma::mean(arma::pow(y - fitted, 2));
        arma::colvec log_losses = -y % arma::log(fitted) - (1 - y) % arma::log(1 -
            fitted);
        log_losses.replace(arma::datum::nan, 0);
        log_loss = arma::sum(log_losses);
        double misclassification = 0;
        arma::colvec misclass = arma::zeros<arma::colvec>(n);
        for (int i = 0; i < p; i++){
          if (classes[i] != y[i]){
            misclass[i] = 1;
          }
        }
        misclassification = arma::mean(misclass);

    }


    return List::create(Named("fitted.values") = fitted,
                        Named("class") = classes,
                        Named("MSE") = MSE,
                        Named("log.loss") = log_loss,
                        Named("misclassification") = misclassification);

}



////------------------------------------------------------------------------------------




//' @title Predict Linear Regression
//' @description Generates prediction for linear regression
//'
//' @param betas 'linearr' object or matrix of betas
//' @param X matrix of (new) observations
//' @param y matrix of response values
//' @return predictions and loss metrics
//' @export
//' @examples
//'
//' fitted = linearr(X, y, penalty = "ridge")
//' predict_linearr(fitted$coefficients, X)
//'
// [[Rcpp::export]]
List predict_linearc(const arma::colvec& betas, const arma::mat& X, const arma::colvec& y = 0) {

  // checks
  if (y.n_rows > 1) {
    if (X.n_rows != y.n_rows)
      stop("X and y must have equal observations!");
  }

  // add intercept, if needed
  arma::mat X_ = X;
  int n = X_.n_rows, p = X_.n_cols;
  if (p != betas.n_rows) {
    arma::colvec X1 = arma::ones<arma::colvec>(n);
    X_.insert_cols(0, X1);
    p = X_.n_cols;
  }

  // fitted values
  arma::colvec fitted = X_ * betas;
  arma::colvec classes = arma::round(fitted);

  // if y, return MSE, misclassification
  double MSE = 0;
  double RSS = 0;
  if (y.n_rows > 1) {

    // calculate metrics
    MSE = arma::mean(arma::pow(y - fitted, 2));
    RSS = arma::sum(arma::pow(y - fitted, 2));

  }


  return List::create(Named("fitted.values") = fitted,
                      Named("MSE") = MSE,
                      Named("RSS") = RSS);
}

