mm= c(1, 0, 0, 0)
mm.lag= c(1, 0, 0, 0)
gamma.mat = matrix(c(10.452300, -4.999838, -4.999635,  5.452912,  7.861847, -3.060784,  5.011114,  1.091025,  8.134455,  0.000000,  0.000000,  0.000000), nrow=3)

source("fun_delta_it_Rcpp.R")

delta_it(mm, mm.lag, gamma.mat, trace=1, maxit=1000)