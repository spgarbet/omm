
sourceCpp("src/dpidtheta_calc.cpp")

dpidtheta_calc <- function(deta.k.dtheta, dpi.deta)
  dpidtheta_calc_cpp(deta.k.dtheta, dpi.deta) 
