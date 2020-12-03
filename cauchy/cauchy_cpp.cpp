#include<Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
NumericVector dcauchy_cpp(
    const NumericVector x,
    const double loc,
    const double scale
){
    NumericVector d = x-loc;
    d = scale / (M_PI * (d*d + scale*scale));
    return d;
}