#ifndef LOGISTIC_H
#define LOGISTIC_H

#include <RcppArmadillo.h>
#include <Rcpp.h>

Rcpp::List logisticc(const arma::mat& X, const arma::colvec& y, double lam = 0, double alpha = 1.5, std::string penalty = "none", bool intercept = true, std::string method = "IRLS", double tol = 1e-5, double maxit = 1e5, arma::colvec vec = 0, arma::colvec init = 0);

#endif