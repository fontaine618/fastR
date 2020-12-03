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
