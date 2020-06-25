# Rcpp interfers with inline and crashes R
# Doing this last seems to work
library(Rcpp)
# For Local compile
Sys.setenv(PKG_LIBS="-llapack")
# For ACCRE compile
#Sys.setenv(PKG_LIBS="-L${MKLROOT}/lib/intel64 -lmkl_rt -lpthread -lm -ldl")
#Sys.setenv(PKG_CXXFLAGS="-I${MKLROOT}/include")

sourceCpp("~/rsch/Git/omm-master/src/dhd_gamma_calc.cpp")

dhdgamma_mmlag_calc <- function(hmat, mm.lag, tmp.mat)
    dhd_gamma_calc(hmat, mm.lag, tmp.mat) 