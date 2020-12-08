source("../functions/FunctionsRcpp.R")
source("../functions/LikelihoodFunctions.R")

################################################
################################################
################################################

simulation <- function(iter){
    
est <- covar <- conv <- list()
est1 <- covar1 <- conv1 <- list()
est2 <- covar2 <- conv2 <- list()
est3 <- covar3 <- conv3 <- list()
    
N         <- 1000
mi        <- 20
id        <- rep(1:N, each=mi)
tx        <- rep(rep(c(0,1), each=mi), N/2)
t         <- rep(seq(0,1,length.out=mi), N)
#XMat      <- cbind(tx,t,tx*t)
XMat      <- cbind(tx,t)

alpha     <-c(-3.25,-.5, 3)
#beta      <- c(0.25,0.5,0.25)
beta      <- c(0,0.5)

beta.ppo  <- c(-0.5, -2.5)

gamma.mat <- rbind( c(20,2,1),
                    c(-5,5,2),
                    c(-5,1,5))

inits        <- c(alpha, beta, rep(0,9))
params       <- c(alpha, beta, c(t(gamma.mat)))
params1      <- c(alpha, beta, beta.ppo, c(t(gamma.mat)))
ppo.mat.indx <- matrix(c(2,2,3,2), byrow=TRUE, ncol=2)

for (blah in 1:4){
Y       <- GenDatOMTM1(id=id, XMat=XMat, alpha=alpha, beta=beta, gamma.mat=gamma.mat)
print(date())
mod          <- nlm(logLikeCalc4, params,  yval=Y, x=XMat, id=id, UseGrad=TRUE, stepmax=.25, iterlim=250, check.analyticals = FALSE, print.level=0, hessian=FALSE)
print(date())
mod1         <- nlm(logLikeCalc4, params1,  yval=Y, x=XMat, id=id, UseGrad=TRUE, ppo.mat.indx=ppo.mat.indx, stepmax=.25, iterlim=250, check.analyticals = FALSE, print.level=0, hessian=FALSE)
print(date())

Y1       <- GenDatOMTM1.ppo(id=id, XMat=XMat, alpha=alpha, beta=beta, gamma.mat=gamma.mat, ppo.mat.indx=ppo.mat.indx, beta.ppo=beta.ppo)
print(date())
mod2          <- nlm(logLikeCalc4, params,  yval=Y1, x=XMat, id=id, UseGrad=TRUE, stepmax=.25, iterlim=250, check.analyticals = FALSE, print.level=0, hessian=FALSE)
print(date())
mod3         <- nlm(logLikeCalc4, params1,  yval=Y1, x=XMat, id=id, UseGrad=TRUE, ppo.mat.indx=ppo.mat.indx, stepmax=.25, iterlim=250, check.analyticals = FALSE, print.level=0, hessian=FALSE)
print(date())

library(VGAM)
mod.vgam1 <- vgam(Y1 ~ t, cumulative(reverse=FALSE, parallel=FALSE))
coef(mod.vgam1)

## Calculate Variance-Covariance Matrices
print(date())
mod.covar  = CalcVarCov(MOD=mod,  epsilon=1e-7,  yval=Y,  x=XMat,  id=id, UseGrad=TRUE)
mod1.covar = CalcVarCov(MOD=mod1, epsilon=1e-7,  yval=Y, x=XMat, id=id, UseGrad=TRUE, ppo.mat.indx=ppo.mat.indx)
mod2.covar = CalcVarCov(MOD=mod2, epsilon=1e-7,  yval=Y1, x=XMat, id=id, UseGrad=TRUE)
mod3.covar = CalcVarCov(MOD=mod3, epsilon=1e-7,  yval=Y1, x=XMat, id=id, UseGrad=TRUE, ppo.mat.indx=ppo.mat.indx)
print(date())

est[[blah]] <- mod$estimate
covar[[blah]] <- mod.covar
conv[[blah]] <- mod$code
est1[[blah]] <- mod1$estimate
covar1[[blah]] <- mod1.covar
conv1[[blah]] <- mod1$code
est2[[blah]] <- mod2$estimate
covar2[[blah]] <- mod2.covar
conv2[[blah]] <- mod2$code
est3[[blah]] <- mod3$estimate
covar3[[blah]] <- mod3.covar
conv3[[blah]] <- mod3$code
}

Results <- list(est=est, covar=covar, conv=conv,
                est1=est1, covar1=covar1, conv1=conv1,
                est2=est2, covar2=covar2, conv2=conv2,
                est3=est3, covar3=covar3, conv3=conv3)
save(Results, file=paste('output/i',iter,'.rda', sep=""), version=2)
}

