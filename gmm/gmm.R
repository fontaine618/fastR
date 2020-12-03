setwd("./gmm")

set.seed(1)
n = 1e6
k = 3
mu = seq(k)*5
sig = rep(1, k)
z = sample.int(k, n, TRUE)
x = rnorm(n)
x = x * sig[z] + mu[z]

system("R CMD SHLIB gmm_em_fortran.f90")
dyn.load("gmm_em_fortran.so")
source("gmm_em.R")

library(Rcpp)
sourceCpp("gmm_em_cpp.cpp")


bench::mark(
    GMM.EM.R(x, 3),
    GMM.EM.Fortran(x, 3),
    GMM.EM.Cpp(x, 3)
)

