
# For ACCRE compile
Sys.setenv(PKG_LIBS="-L${MKLROOT}/lib/intel64 -lmkl_rt -lpthread -lm -ldl")
Sys.setenv(PKG_CXXFLAGS="-I${MKLROOT}/include")

source("FourStatesOneAbsorbing.R")


  #############################################################################
 ##
## ACCRE batch run
args <- commandArgs(trailingOnly=TRUE)
x    <- as.numeric(args[1])
set.seed(x)
simulation(x)
