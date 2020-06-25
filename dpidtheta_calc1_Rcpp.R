# Rcpp interfers with inline and crashes R
# Doing this last seems to work
library(Rcpp)
# For Local compile
Sys.setenv(PKG_LIBS="-llapack")
# For ACCRE compile
#Sys.setenv(PKG_LIBS="-L${MKLROOT}/lib/intel64 -lmkl_rt -lpthread -lm -ldl")
#Sys.setenv(PKG_CXXFLAGS="-I${MKLROOT}/include")

sourceCpp("~/rsch/Git/omm-master/src/dpidtheta_calc1.cpp")

dpidtheta_calc1 <- function(diagmtx, dpi.deta, xit)
  dpidtheta_calc1_cpp(diagmtx, dpi.deta, xit) 
