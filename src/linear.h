#ifndef LINEAR_H
#define LINEAR_H

#include <RcppArmadillo.h>
#include <Rcpp.h>

Rcpp::List linearc(const arma::mat& X, const arma::colvec& y, double lam = 0, arma::colvec weights = 0, bool intercept = true, bool kernel = false);

#endif