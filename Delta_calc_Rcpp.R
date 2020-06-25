
sourceCpp("src/Delta_calc.cpp")

Delta_calc <- function(mm, mm.lag, gamma.mat, tol=1e-4, maxit=10000, trace=0)
  Delta_calc_cpp(mm, mm.lag, gamma.mat, tol, maxit, trace) 
