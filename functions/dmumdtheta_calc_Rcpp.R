
library(Rcpp)


sourceCpp("../src/dmumdtheta_calc.cpp")

dmumdtheta_calc <- function(dpi.dtheta)
  dmumdtheta_calc_cpp(dpi.dtheta) 
