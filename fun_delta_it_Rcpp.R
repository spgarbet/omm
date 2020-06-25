
sourceCpp("src/delta_it.cpp")

delta_it <- function(mm, mm.lag, gamma.mat, tol=1e-4, maxit=200, trace=0)
  delta_it_cpp(mm, mm.lag, gamma.mat, tol, maxit, trace) 
