
library(inline)
inner_c <- paste0(readLines("src/delta_it.c"), collapse="\n")

# This defines the interface between R and C
inner   <- cfunction(
  signature( k        = "integer", # Categories
             gammamat = "double",  # (k-1, k)
             mm       = "double",  # k-1 (okay to pass full k)
             mmlag    = "double",  # k
             tol      = "double",  # Convergence criteria
             maxit    = "integer", # iteration cutoff
             trace    = "integer", # debugging level, 0 or 1
    
             # These are local variables, R is allocating
             hmat     = "double",  # (k-1, k)
             dfdDelta = "double",  # (k-1, k-1)
             del      = "double",  # k-1
             fDelta   = "double",  # k-1
             Deltavec = "double",  # k-1
             ipiv     = "integer", # k-1
             work     = "double"   # (k-1, k-1)
           ),
           inner_c,
           includes   = c("#include <R_ext/Lapack.h>","#include <R_ext/BLAS.h>"),
           language   = "C",
  
           # Local Compile
           libargs    = "-llapack",
           # -----
           
           # ACCRE
           #libargs = "-L${MKLROOT}/lib/intel64 -lmkl_rt -lpthread -lm -ldl",
           #cppargs = "-I${MKLROOT}/include",
           # -----
  
           convention = ".C")

# Let's put an R wrapper on it that does the memory
delta_it_c <- function(mm,              # row of marginal (mean) probabilities
                                        # by response category k=1...K 
                       mm.lag,          # lagged value of mm
                       gamma.mat,       # Transition log odds ratios
                       k          = length(mm.lag),
                       tol        = 1e-4,
                       maxit      = 100,
                       trace      = 0)
{
  # Using inline all inputs are outputs, but we only want one
  inner(k, gamma.mat, mm, mm.lag, tol, maxit, trace,
        rep(0, k*(k-1)), rep(0, (k-1)*(k-1)), rep(0,k-1), rep(0,k-1), rep(0, k-1), rep(0, k-1), rep(0, (k-1)*(k-1)))$Deltavec
}