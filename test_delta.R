  ##############################################################################
 #
#   This file will explore the differences between
#   implementations of delta_it in R, Rcpp, and inline C.
#

print("Beginning Test")

source("fun_delta_it_R.R")
load("data/delta_it-input.Rdata") # Simulation data point
reference <- delta_it_r(mm, mm.lag, gamma.mat)
orig.time   <- system.time(replicate(10000, delta_it_r(mm, mm.lag, gamma.mat)))

source("fun_delta_it_Rcpp.R")
rcppv     <- delta_it_cpp(mm, mm.lag, gamma.mat, 1e-4, 100, 0)
if(abs(sum(reference - rcppv)) > 1e-9) warning("Rcpp Version Broken!")
rcpp.time <- system.time(
  replicate(10000, delta_it_cpp(mm, mm.lag, gamma.mat, 1e-4, 100, 0)))

print("R Code")
print(orig.time)

print("Rcpp")
print(rcpp.time)

cat(paste("A", round(orig.time[1]/rcpp.time[1]), " fold increase\n"))