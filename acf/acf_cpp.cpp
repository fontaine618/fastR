#include<Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
NumericVector acf_cpp(
    const NumericVector x,
    const int lag
){
    int n = x.length();
    NumericVector ac(lag+1);
    ac[0] = 1.;
    double mu = mean(x);
    double sig = sd(x);
    NumericVector xc = (x-mu)/sig;
    for(int k=1; k<=lag; k++){
        ac[k] = sum(xc[Range(0,(n-k-1))]*xc[Range(k, n-1)]) / (n-1);
    }
    return ac;
}