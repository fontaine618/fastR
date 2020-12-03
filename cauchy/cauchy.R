setwd("./cauchy")

dcauchy_r = function(x, loc, scale){
    return(scale / (pi * ((x-loc)^2 + scale^2)))
}

Rcpp::sourceCpp("cauchy_cpp.cpp")

system("R CMD SHLIB cauchy_f.f90")
dyn.load("cauchy_f.so")

dcauchy_f = function(x, loc, scale){
    n = length(x)
    d = .Fortran("dcauchy",
        n=as.integer(n),
        x=as.double(x),
        loc=as.double(loc),
        scale=as.double(scale),
        d=double(n)
    )$d
    return(d)
}

set.seed(1)
x = rnorm(1e8, 0, 1)
loc = 0
scale = 1

bench::mark(
    dcauchy(x, loc, scale),
    dcauchy_r(x, loc, scale),
    dcauchy_cpp(x, loc, scale),
    dcauchy_f(x, loc, scale)
)
