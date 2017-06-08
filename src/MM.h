#ifndef MM_H
#define MM_H

#include <RcppArmadillo.h>
#include <Rcpp.h>

Rcpp::List MMc(const arma::mat& X, const arma::colvec& y, double lam = 0, double alpha = 1.5, double gamma = 1, bool intercept = true, double tol = 1e-5, double maxit = 1e5, const arma::colvec& vec = 0, const arma::colvec& init = 0);

Rcpp::List MM_linearc(const arma::mat& X, const arma::colvec& y, double lam = 0, double alpha = 1.5, double gamma = 1, arma::colvec weights = 0, bool intercept = true, double tol = 1e-5, double maxit = 1e5, const arma::colvec& vec = 0, const arma::colvec& init = 0);

#endif
