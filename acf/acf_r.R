setwd("./acf")

acf_c = function(x, lag){
    ac = acf(x, lag, plot=FALSE)$acf
    return(as.vector(ac))
}

acf_r = function(x, lag.max=10){
    n = length(x)
    ac = rep(0, lag.max+1)
    ac[1] = 1
    mu = mean(x)
    sig = sd(x)
    xc = (x-mu)/sig
    for(k in seq(lag.max)){
        ac[k+1] = sum(xc[1:(n-k)]*xc[(k+1):n]) / (n-1)
    }
    return(ac)
}

Rcpp::sourceCpp("acf_cpp.cpp")

system("R CMD SHLIB acf_f.f90")
dyn.load("acf_f.so")

acf_f = function(x, lag){
    n = length(x)
    ac = .Fortran("acf", PACKAGE="acf_f",
        n=as.integer(n),
        x=as.double(x),
        lag=as.integer(lag),
        ac=double(lag+1)
    )$ac
    return(ac)
}

set.seed(1)
x = rnorm(1e7)
lag = 1000

bench::mark(
    acf_c(x, lag),
    acf_r(x, lag),
    acf_cpp(x, lag),
    acf_f(x, lag)
)