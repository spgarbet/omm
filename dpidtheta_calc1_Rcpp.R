

sourceCpp("src/dpidtheta_calc1.cpp")

dpidtheta_calc1 <- function(diagmtx, dpi.deta, xit)
  dpidtheta_calc1_cpp(diagmtx, dpi.deta, xit) 
