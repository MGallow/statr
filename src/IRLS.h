#ifndef IRLS_H
#define IRLS_H

#include <RcppArmadillo.h>
#include <Rcpp.h>

arma::colvec logitc(const arma::colvec& u);

Rcpp::List IRLSc(const arma::mat& X, const arma::colvec& y, double lam = 0, std::string penalty = "none", bool intercept = true, double tol = 1e-5, double maxit = 1e5, const arma::colvec& vec = 0, const arma::colvec& init = 0);


#endif
