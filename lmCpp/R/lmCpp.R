lm_cpp = function(X, Y, eps=1e-10, max_iter=1e5, lr=0.01){
  beta = lm_cpp_fit(X, Y, eps, max_iter, lr)
  return(as.vector(beta))
}