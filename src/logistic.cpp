// Matt Galloway

#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "IRLS.h"
#include "MM.h"

using namespace Rcpp;

//' @title Logistic Regression (c++)
//' @description Computes the coefficient estimates for logistic regression. ridge regularization and bridge regularization optional.
//'
//' @param X matrix
//' @param y matrix or vector of response values 0,1
//' @param lam optional tuning parameter for ridge regularization term. Defaults to `lam = 0`
//' @param alpha optional tuning parameter for bridge regularization term. Defaults to 'alpha = 1.5'
//' @param penalty choose from c('none', 'ridge', 'bridge'). Defaults to 'none'
//' @param intercept Defaults to TRUE
//' @param method optimization algorithm. Choose from 'IRLS' or 'MM'. Defaults to 'IRLS'
//' @param tol tolerance - used to determine algorithm convergence. Defaults to 1e-5
//' @param maxit maximum iterations. Defaults to 1e5
//' @param vec optional vector to specify which coefficients will be penalized
//' @param init optional initialization for MM algorithm
//' @return returns beta estimates (includes intercept), total iterations, and gradients.
//' @export
//' @examples
//' Logistic Regression
//' library(dplyr)
//' X = as.matrix(dplyr::select(iris, -Species))
//' y = as.matrix(dplyr::select(iris, Species))
//' y = ifelse(y == 'setosa', 1, 0)
//' logisticc(X, y, vec = c(0,1,1,1))
//'
//' ridge Logistic Regression with IRLS
//' logisticc(X, y, lam = 0.1, penalty = 'ridge', vec = c(0,1,1,1))
//'
//' ridge Logistic Regression with MM
//' logisticc(X, y, lam = 0.1, penalty = 'ridge', method = 'MM', vec = c(0,1,1,1))
//'
//' bridge Logistic Regression
//' logisticc(X, y, lam = 0.1, alpha = 1.5, penalty = 'bridge', method = "MM", vec = c(0,1,1,1))
//'
// [[Rcpp::export]]
List logisticc(const arma::mat& X, const arma::colvec& y, double lam = 0, double alpha = 1.5, std::string penalty = "none", bool intercept = true, std::string method = "IRLS", double tol = 1e-5, double maxit = 1e5, arma::colvec vec = 0, arma::colvec init = 0) {

  // checks
  int p = X.n_cols;
  if ((alpha >= 2 | alpha <= 1)){
    stop("alpha must be between 1 and 2!");
  }
  if (lam < 0){
    stop("lam must be nonnegative!");
  }


  // if IRLS algorithm...
  arma::colvec betas = arma::ones<arma::colvec>(p);
  arma::colvec grads = arma::ones<arma::colvec>(p);
  int iterations = 1;
  if (method == "IRLS") {

    // execute IRLS script
    List logistic = IRLSc(X, y, lam, penalty, intercept, tol, maxit, vec, init);

    iterations = as<int>(logistic["total.iterations"]);
    betas = as<NumericVector>(logistic["coefficients"]);
    grads = as<NumericVector>(logistic["gradient"]);
    if (iterations == maxit) {
      Rprintf("Algorithm did not converge...Try MM\n");
    }
  }

  // if MM algorithm...
  if (method == "MM") {

    // change gamma parameter, if needed
    int gamma = 1;
    if (penalty == "bridge"){
      gamma = 0;
    }

    // execute MM script
    List logistic = MMc(X, y, lam, alpha, gamma, intercept, tol, maxit, vec, init);

    iterations = as<int>(logistic["total.iterations"]);
    betas = as<NumericVector>(logistic["coefficients"]);
    grads = as<NumericVector>(logistic["gradient"]);
    if ((iterations == maxit) & (penalty != "bridge")) {
      Rprintf("Algorithm did not converge...Try IRLS\n");
    }
  }

  return List::create(Named("coefficients") = betas,
                      Named("total.iterations") = iterations,
                      Named("gradient") = grads);

}


