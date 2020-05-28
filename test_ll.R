source("delta_it_Rcpp.R")
a
gamma.mat <- rbind( c(3,2,1),
                    c(1,2,1),
                    c(-0.5,1,2))
params    <- c(alpha, beta, c(t(gamma.mat)))

Y         <- GenDatOMTM1(id        = id,
                         XMat      = XMat,
                         alpha     = alpha,
                         beta      = beta,
                         gamma.mat = gamma.mat)

   
print("Generated Data")

#mod  <- optim( params, logLikeCalc, yval=Y, x=XMat, id=id, control=list(trace=0, maxit=10000), hessian=TRUE)

### logLikeCalc does not use the gradient

print(system.time(nlm(logLikeCalc, params,  yval=Y, x=XMat, id=id, stepmax=.1, iterlim=250, check.analyticals = TRUE, print.level=0, hessian=TRUE)))

### logLikeCalc3 uses the gradient

print(system.time(nlm(logLikeCalc3, params,  yval=Y, x=XMat, id=id, stepmax=.1, iterlim=250, check.analyticals = TRUE, print.level=0, hessian=TRUE)))
