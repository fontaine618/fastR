#include<Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
List gmm_em_cpp(
    const NumericVector x,
    NumericVector mu,
    NumericVector sig,
    NumericVector pi,
    int max_iter,
    const double eps
){
    const int n = x.length();
    const int k = mu.length();
    NumericMatrix p(n, k);
    double llk = -9.9e30;
    double prev_llk;
    double tmp;
    List out;

    for (size_t l = 0; l < max_iter; l++)
    {
        // E Step
        for (size_t j = 0; j < k; j++)
        {
            p(_, j) = dnorm(x, mu[j], sig[j]) * pi[j];
        }
        // Check convergence
        prev_llk = llk;
        llk = 0.;
        for (size_t i = 0; i < n; i++)
        {
            tmp = sum(p(i, _));
            llk += log(tmp);
            p(i, _) = p(i, _ ) / tmp;
        }
        llk = llk / n;
        if (llk - prev_llk < eps){
            out["max_iter"] = l + 1;
            break;
        }       
        //M Step
        for (size_t j = 0; j < k; j++)
        {
            tmp = sum(p(_, j));
            mu[j] = sum(p(_, j) * x) / tmp;
            sig[j] = sqrt(sum(p(_, j) * (x-mu[j]) * (x-mu[j])) / tmp);
            pi[j] = tmp / n;
        }

    }   

    out["llk"] = llk;
    return(out);
}