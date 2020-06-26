

sourceCpp("../src/hmat_calc.cpp")

hmat_calc <- function(Delta.vec, gamma.mat)
  hmat_calc_cpp(Delta.vec, gamma.mat) 
