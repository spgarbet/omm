# Rcpp interfers with inline and crashes R
# Doing this last seems to work
library(Rcpp)
# For Local compile
Sys.setenv(PKG_LIBS="-llapack")
# For ACCRE compile
#Sys.setenv(PKG_LIBS="-L${MKLROOT}/lib/intel64 -lmkl_rt -lpthread -lm -ldl")
#Sys.setenv(PKG_CXXFLAGS="-I${MKLROOT}/include")

sourceCpp("~/rsch/Git/OMM/src/hmat_calc.cpp")

hmat_calc <- function(Delta.vec, gamma.mat)
  hmat_calc_cpp(Delta.vec, gamma.mat) 
