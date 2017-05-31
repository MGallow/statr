// Matt Galloway

#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "logistic.h"

using namespace Rcpp;




//' @title K fold (c++)
//' @description creates vector of shuffled indices
//'
//' @param n number of eleemtns
//' @param K number of folds
//' @return returns vector
//' @examples
//' kfold(10, 3)
//'
// [[Rcpp::export]]
arma::vec kfold(int n, int K){

  // create sequence 1:n
  arma::vec indices = arma::linspace<arma::vec>(1, n, n);

  // assign number fold
  for (int i = 0; i < n; i ++){
    indices[i] = i % K;
  }

  // shuffle indices
  indices = arma::shuffle(indices);

  return indices;

}



//--------------------------------------------------------------------------------------------


//' @title CV Logisticc (c++)
//' @description Computes the coefficient estimates for logistic regression. ridge regularization and bridge regularization optional.
//'
//' @param X matrix or data frame
//' @param y matrix or vector of response values 0,1
//' @param lam vector of tuning parameters for ridge regularization term. Defaults to `lam = 0`
//' @param alpha vector of tuning parameters for bridge regularization term. Defaults to 'alpha = 1.5'
//' @param penalty choose from c('none', 'ridge', 'bridge'). Defaults to 'none'
//' @param vec optional vector to specify which coefficients will be penalized
//' @param intercept Defaults to TRUE
//' @param method optimization algorithm. Choose from 'IRLS' or 'MM'. Defaults to 'IRLS'
//' @param tol tolerance - used to determine algorithm convergence. Defaults to 1e-5
//' @param maxit maximum iterations. Defaults to 1e5
//' @param init initialization of betas
//' @param K specify number of folds in cross validation, if necessary
//'
//' @return returns best lambda, best alpha, cv.errors
//' @export
//' @examples
//' CV Logistic Regression
//'
// [[Rcpp::export]]
List CV_logisticc(const arma::mat& X, const arma::colvec& y, const arma::colvec& lam = 0, const arma::colvec& alpha = 0, std::string penalty = "none", bool intercept = true, std::string method = "IRLS", double tol = 1e-5, double maxit = 1e4, arma::colvec vec = 0, int K = 3) {


  // initialization
  int n = X.n_rows;
  arma::mat CV_errors = arma::zeros<arma::mat>(lam.n_rows, alpha.n_rows);

  // designate folds and shuffle -- ensures randomized folds
  arma::vec folds = kfold(n, K);

  // parse data into folds and perform CV
  for (int k = 0; k < K; k++){

    // separate into training and testing data
    arma::uvec index = arma::find(folds == k);
    arma::uvec index_ = arma::find(folds != k);

    arma::mat X_train = X.rows(index);
    arma::mat X_test = X.rows(index_);

    arma::colvec y_train = y.rows(index);
    arma::colvec y_test = y.rows(index_);

    // loop over all tuning parameters
    arma::mat CV_error = arma::zeros<arma::mat>(lam.n_rows, alpha.n_rows);
    for (int i = 0; i < lam.n_rows; i++){
      for (int j = 0; j < alpha.n_rows; j++){

        // set temporary tuning parameters
        double lam_ = lam[i];
        double alpha_ = alpha[j];

        // generate new coefficients
        List logistic = logisticc(X_train, y_train, lam_, alpha_, penalty, intercept, method, tol, maxit, vec);
        arma::colvec betas = as<NumericVector>(logistic["coefficients"]);

        // generate predictions
        arma::colvec fitted = X_test * betas;

        // errors
        double error = arma::sum(arma::pow(fitted - y_test, 2));

        // input error to CV_error
        CV_error(i, j) = error;
      }
    }

    // append CV errors
    CV_errors += CV_error;

  }

  // determine optimal tuning parameters
  arma::uword ind = CV_errors.index_min();
  int lam_ind = ind % CV_errors.n_rows;
  int alpha_ind = floor(ind/CV_errors.n_rows);
  double best_lam = lam[lam_ind];
  double best_alpha = alpha[alpha_ind];


  // return list of coefficients
  return List::create(Named("best.lam") = best_lam,
                      Named("best.alpha") = best_alpha,
                      Named("cv.errors") = CV_errors);
}

