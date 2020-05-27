source("Functions.R")



N         <- 25
mi        <- 20
id        <- rep(1:N, each=mi)
tx        <- rep(rep(c(0,1), each=mi), N/2)
t         <- rep(seq(0,1,length.out=mi), N)
XMat      <- cbind(tx,t,tx*t)
alpha     <-c(-2,-1,1)
beta      <- c(0.25,0.5,-.25)
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
