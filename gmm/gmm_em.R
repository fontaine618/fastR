GMM.init = function(x, k){
  r = range(x)
  mu0 = seq(from=r[1], to=r[2], length.out=k)
  sig0 = rep(diff(r), k) / (2*k)
  return(list(mu0=mu0, sig0=sig0))
}

GMM.EM.R = function(x, k=2, eps=1.0e-8, max_iter=1000){
  init = GMM.init(x, k)
  mu = init$mu0; sig = init$sig0; pi = rep(1/k, k)
  n = length(x)
  p = matrix(0, n, k)
  llk = -Inf
  for(i in seq(max_iter)){
    # E STEP
    for(j in seq(k)){
      p[, j] = dnorm(x, mu[j], sig[j]) * pi[j]
    }
    # CHECK CONVERGENCE
    prev_llk = llk
    llk = mean(log(apply(p, 1, sum)))
    if(llk - prev_llk < eps) break
    p = sweep(p, 1, apply(p, 1, sum), "/")
    # M STEP
    pi = apply(p, 2, mean)
    mu = apply(sweep(p, 1, x, "*"), 2, sum) / apply(p, 2, sum)
    c = outer(x, mu, "-")
    sig2 = apply(p * c^2, 2, sum) / apply(p, 2, sum)
    sig = sqrt(sig2)
  }
  return(c(i, llk))
}

GMM.EM.Fortran = function(x, k=2, eps=1.0e-8, max_iter=1000){
  init = GMM.init(x, k)
  mu = init$mu0; sig = init$sig0; pi = rep(1/k, k)
  n = length(x)
  out = .Fortran(
    "gmm_em", PACKAGE="gmm_em_fortran",
    n=as.integer(n),
    k=as.integer(k),
    x=as.double(x),
    mu=as.double(mu),
    sig2=as.double(sig^2),
    pi=as.double(pi),
    llk=double(1),
    max_iter=as.integer(max_iter),
    eps=as.double(eps)
  )
  return(c(out$max_iter, out$llk))
}

GMM.EM.Cpp = function(x, k=2, eps=1.0e-8, max_iter=1000){
  init = GMM.init(x, k)
  mu = init$mu0; sig = init$sig0; pi = rep(1/k, k)
  n = length(x)
  out = gmm_em_cpp(x, mu, sig, pi, max_iter, eps)
  return(c(out$max_iter, out$llk))
}