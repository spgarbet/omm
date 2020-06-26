# For Local compile
Sys.setenv(PKG_LIBS="-llapack")

source("FourStatesOneAbsorbing.R")

  #############################################################################
 ##
## For exection on local desktop
library(parallel)
 
mclapply(1:1, mc.cores=8, function(x)
{
  set.seed(x)
  simulation(x)
})

