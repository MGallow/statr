#ifndef LINEAR_H
#define LINEAR_H

#include <RcppArmadillo.h>
#include <Rcpp.h>

Rcpp::List linearc(const arma::mat& X, const arma::colvec& y, double lam = 0, double alpha = 1.5, std::string penalty = "none", arma::colvec weights = 0, bool intercept = true, bool kernel = false, std::string method = "SVD", double tol = 1e-5, double maxit = 1e5, arma::colvec vec = 0, arma::colvec init = 0);

#endif