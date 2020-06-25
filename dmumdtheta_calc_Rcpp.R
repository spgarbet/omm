# Rcpp interfers with inline and crashes R
# Doing this last seems to work
library(Rcpp)

sourceCpp("src/dmumdtheta_calc.cpp")

dmumdtheta_calc <- function(dpi.dtheta)
  dmumdtheta_calc_cpp(dpi.dtheta) 
