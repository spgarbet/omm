# Rcpp interfers with inline and crashes R
# Doing this last seems to work
library(Rcpp)
# For Local compile
Sys.setenv(PKG_LIBS="-llapack")
# For ACCRE compile
#Sys.setenv(PKG_LIBS="-L${MKLROOT}/lib/intel64 -lmkl_rt -lpthread -lm -ldl")
#Sys.setenv(PKG_CXXFLAGS="-I${MKLROOT}/include")

sourceCpp("~/rsch/Git/omm-master/src/dfdDelta_calc.cpp")

dfdDelta_calc <- function(mm.lag, hmat)
  dfdDelta_calc_cpp(mm.lag, hmat) 
