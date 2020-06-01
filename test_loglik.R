source("fun_delta_it_Rcpp.R")
source("fun_loglik.R")
source("generate_data.R")

set.seed(1)

N         <- 26
mi        <- 20
id        <- rep(1:N, each=mi)
tx        <- rep(rep(c(0,1), each=mi), N/2)
t         <- rep(seq(0,1,length.out=mi), N)
XMat      <- cbind(tx,t,tx*t)


alpha     <- c(-2,-1, 1)
beta      <- c(0.25, 0.5, -.25)
gamma.mat <- rbind( c( 3,  2,  1),
                    c( 1,  2,  1),
                    c(-0.5,1,  2))
params    <- c(alpha, beta, c(t(gamma.mat)))

Y         <- GenDatOMTM1(id        = id,
                         XMat      = XMat,
                         alpha     = alpha,
                         beta      = beta,
                         gamma.mat = gamma.mat)

   
print("Generated Data") 

# 1.6 is the time to beat
print(system.time(ref<-nlm(logLikeCalc,   params,      
                           yval=Y,        x=XMat,      id=id,
                           stepmax=.1,    iterlim=250, check.analyticals=TRUE,
                           hessian=TRUE)))

#print(system.time(nlm(logLikeCalcWithGrad, params,  yval=Y, x=XMat, id=id, stepmax=.1, iterlim=250, check.analyticals = TRUE, print.level=0, hessian=TRUE)))

source("fun_loglik_experimental.R")


    States    <- sort(unique(Y))
    K         <- length(States)
    K1        <- K-1
    uid       <- unique(id)
    N         <- length(uid)
    npar      <- length(params)
    
print(system.time(exp<-nlm(logLikeCalcExp,params,      
                           yval=Y,        x=XMat,      id=id,
                           States=States, K=K, K1=K1,  uid=uid,
                           N=N, npar=npar,
                           stepmax=.1,    iterlim=250, check.analyticals=TRUE,
                           hessian=TRUE)))


if(sum(abs(ref$estimate - exp$estimate)) > 1e-6)
  warning("Experimental Results Deviate Too Much")