# This assumes that the root directory code is running from is at the same
# level as the functions directory. 
# E.g:
# project
# ├── functions
# │   ├── functions.R
# │   ├── hmat_calc_Rcpp.R
# │   └── ...
# ├── src
# │   ├── hmat_calc.cpp
# │   └── ...
# ├── sim1   # Root Directory, if in project setwd("sim1") to start
# │   ├── FourStatesOneAbsorbing.R
# │   ├── FourStatesOneAbsorbing.slurm
# │   ├── sim-accre.R
# │   ├── sim-local.R
# │   └── ...

# This file loads the current set of functions for OMM

library(Rcpp)

source("../functions/hmat_calc_Rcpp.R")
source("../functions/dfdDelta_calc_Rcpp.R")
source("../functions/Delta_calc_Rcpp.R")
source("../functions/Delta_calcWith0TxProbs_Rcpp.R")
source("../functions/dpidtheta_calc_Rcpp.R")
source("../functions/dpidtheta_calc1_Rcpp.R")
source("../functions/dmumdtheta_calc_Rcpp.R")
source("../functions/dhdgamma_mmlag_calc_Rcpp.R")
source("../functions/LikelihoodFunctions.R")