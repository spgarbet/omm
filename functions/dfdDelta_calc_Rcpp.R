
sourceCpp("../src/dfdDelta_calc.cpp")

dfdDelta_calc <- function(mm.lag, hmat)
  dfdDelta_calc_cpp(mm.lag, hmat) 
