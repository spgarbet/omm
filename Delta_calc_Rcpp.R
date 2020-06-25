# Rcpp interfers with inline and crashes R
# Doing this last seems to work
library(Rcpp)
# For Local compile
Sys.setenv(PKG_LIBS="-llapack")
# For ACCRE compile
#Sys.setenv(PKG_LIBS="-L${MKLROOT}/lib/intel64 -lmkl_rt -lpthread -lm -ldl")
#Sys.setenv(PKG_CXXFLAGS="-I${MKLROOT}/include")

sourceCpp("~/rsch/Git/omm-master/src/Delta_calc.cpp")

Delta_calc <- function(mm, mm.lag, gamma.mat, tol=1e-4, maxit=10000, trace=0)
  Delta_calc_cpp(mm, mm.lag, gamma.mat, tol, maxit, trace) 
