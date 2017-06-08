#ifndef PREDICT_H
#define PREDICT_H

#include <RcppArmadillo.h>
#include <Rcpp.h>

Rcpp::List predict_linearc(const arma::colvec& betas, const arma::mat& X, const arma::colvec& y = 0);

Rcpp::List predict_logisticc(const arma::colvec& betas, const arma::mat& X, const arma::colvec& y = 0);

#endif
