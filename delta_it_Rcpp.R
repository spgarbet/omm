# Rcpp interfers with inline and crashes R
# Doing this last seems to work
library(Rcpp)
# For Local compile
Sys.setenv(PKG_LIBS="-llapack")
# For ACCRE compile
#Sys.setenv(PKG_LIBS="-L${MKLROOT}/lib/intel64 -lmkl_rt -lpthread -lm -ldl")
#Sys.setenv(PKG_CXXFLAGS="-I${MKLROOT}/include")

sourceCpp("src/delta_it.cpp")

delta_it <- function(mm, mm.lag, gamma.mat, tol=1e-4, maxit=100, trace=0)
  delta_it_cpp(mm, mm.lag, gamma.mat, tol, maxit, trace) 
