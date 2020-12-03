setwd("./lm")

lm_r = function(X, Y, eps=1e-10, max_iter=1e5, lr=0.01){
    p = ncol(X)
    n = nrow(X)
    beta = matrix(1, p, 1)
    for(i in seq(max_iter)){
        g = t(X) %*% (Y - X %*% beta) / n
        if(sqrt(sum(g^2)) * lr < eps) break
        beta = beta + lr * g
    }
    return(as.vector(beta))
}

Rcpp::sourceCpp("lm.cpp")
lm_cpp = function(X, Y, eps=1e-10, max_iter=1e5, lr=0.01){
    beta = lm_cpp_fit(X, Y, eps, max_iter, lr)
    return(as.vector(beta))
}

system("R CMD SHLIB lm_f.f90")
dyn.load("lm_f.so")
lm_f = function(X, Y, eps=1e-10, max_iter=1e5, lr=0.01){
    p = ncol(X)
    n = nrow(X)
    beta = .Fortran(
        "lm_f",
        n = as.integer(n),
        p = as.integer(p),
        X = as.double(X),
        Y = as.double(Y),
        eps = as.double(eps),
        max_iter = as.integer(max_iter),
        lr = as.double(lr),
        beta = double(p)
    )$beta
    return(as.vector(beta))
}


library(lmFortran)
library(lmCpp)

set.seed(1)
n = 1e3
p = 1e1
X = matrix(runif(n*p), n, p)
lr = 0.01
beta = matrix(1, p, 1)
Y = X %*% beta + rnorm(n)

bench::mark(
    lm_cpp(X, Y, lr=lr),
    lm_f(X, Y, lr=lr)
)