#setwd("problem")

source("Functions.R")

 

alpha     <-c(-1.5,-0.5,0.5)
beta      <- c(.25,.25,.25)
gamma.mat <- rbind( c(6,4,2),
                    c(2,4,2),
                    c(-2,0,2))

params    <- c(alpha, beta, c(t(gamma.mat)))

load("ShawnDat.rda")

model <- optim( params, logLikeCalc, yval=ShawnDat$y, x=model.matrix(~0+tx*t, data=ShawnDat), id=ShawnDat$id, control=list(trace=0, maxit=10000), hessian=TRUE)

diag(solve(model$hessian))
