  ##############################################################################
 #
#   This file will explore the differences between
#   implementations of delta_it in R, Rcpp, and inline C.
#

print("Beginning Test")

source("fun_delta_it_R.R")
source("fun_delta_it_c.R")

load("data/delta_it-input.Rdata") # Simulation data point

#save(gamma.mat, mm, mm.lag, nlevels, file="~/delta_it-input.Rdata", version=2)

reference <- delta_it_r(mm, mm.lag, gamma.mat)
inlinec   <- delta_it_c(mm, mm.lag, gamma.mat)


# Let's test it
if(abs(sum(reference - inlinec)) > 1e-9)
  warning("C Version Broken!")

orig.time   <- system.time(replicate(5000, delta_it_r(mm, mm.lag, gamma.mat)))
inline.time <- system.time(replicate(5000, delta_it_c(mm, mm.lag, gamma.mat)))

source("fun_delta_it_Rcpp.R")

rcppv     <- delta_it_cpp(mm, mm.lag, gamma.mat, 1e-4, 100, 0)
if(abs(sum(reference - rcppv)) > 1e-9)
  warning("Rcpp Version Broken!")
rcpp.time <- system.time(replicate(5000, delta_it_cpp(mm, mm.lag, gamma.mat, 1e-4, 100, 0)))

print("R Code")
print(orig.time)
print("inline C")
print(inline.time)
print("Rcpp")
print(rcpp.time)