# Rcpp interfers with inline and crashes R
# Doing this last seems to work
library(Rcpp)
# For Local compile
Sys.setenv(PKG_LIBS="-llapack")
# For ACCRE compile
#Sys.setenv(PKG_LIBS="-L${MKLROOT}/lib/intel64 -lmkl_rt -lpthread -lm -ldl")
#Sys.setenv(PKG_CXXFLAGS="-I${MKLROOT}/include")

sourceCpp("~/rsch/Git/omm-master/src/Delta_calcWith0TxProbs.cpp")

Delta_calcWith0TxProbs <- function(mm, mm.lag, gamma.mat, CalcTxMtx, tol=1e-4, maxit=10000, trace=0)
  Delta_calcWith0TxProbs_cpp(mm, mm.lag, gamma.mat, CalcTxMtx, tol, maxit, trace) 
