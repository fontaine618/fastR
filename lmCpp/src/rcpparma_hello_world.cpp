// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::colvec lm_cpp_fit(
        arma::mat X,
        arma::colvec Y,
        double eps,
        int max_iter,
        double lr
){
    int p = X.n_cols;
    int n = X.n_rows;
    arma::colvec beta = arma::zeros<arma::colvec>(p);
    arma::colvec g = arma::zeros<arma::colvec>(p);
    for(int i=0; i<max_iter; i++){
        g = X.t() * (Y - X * beta) / (double)n;
        if(norm(g, 2) * lr < eps) break;
        beta += lr * g;
    }
    return(beta);
}
