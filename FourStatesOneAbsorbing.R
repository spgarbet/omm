library(Rcpp)

source("hmat_calc_Rcpp.R")
source("dfdDelta_calc_Rcpp.R")
source("Delta_calc_Rcpp.R")
source("Delta_calcWith0TxProbs_Rcpp.R")
source("dpidtheta_calc_Rcpp.R")
source("dpidtheta_calc1_Rcpp.R")
source("dmumdtheta_calc_Rcpp.R")
source("dhdgamma_mmlag_calc_Rcpp.R")
source("LikelihoodFunctions.R")

simulation <- function(iter)
{

N         <- 250
mi        <- 25
id        <- rep(1:N, each=mi)
tx        <- rep(rep(c(0,1), each=mi), N/2)
t         <- rep(seq(0,1,length.out=mi), N)
XMat      <- cbind(tx,t,tx*t)
alpha     <-c(-3,-1.5, 0)
beta      <- c(0.25,0.25,.25)

gamma.mat <- rbind( c(8,3,1),
                    c(1,5,1),
                    c(-1,1,5))

inits    <- c(alpha, beta, rep(0,9))
params    <- c(alpha, beta, c(t(gamma.mat)))

est <- covar <- conv <- list()
est2 <- covar2 <- conv2 <- list()

startdt <- date()
for (i in 1:5){ print(c(i, date()))
    Y           <- GenDatOMTM1(id=id,
                               XMat=XMat, 
                               alpha=alpha, 
                               beta=beta, 
                               gamma.mat=gamma.mat,
                               tx=tx,
                               t=t)
    
    L    <- length(Y)
    dup  <- duplicated(id)
    Y1   <- Y
    drop <- rep(0, length(Y))
    for (ww in 2:L){ if (Y1[ww-1]==1 & dup[ww]){ Y1[ww]=1
    drop[ww]=1}}
    
    Y2    <- Y[drop==0]
    XMat2 <- XMat[drop==0,]
    id2   <- id[drop==0]
    
    ###########################################  
    Prof <- c(NA)
    mod  <- nlm(logLikeCalc4, params,  yval=Y, x=XMat, id=id, UseGrad=TRUE, stepmax=1, iterlim=250, check.analyticals = FALSE, print.level=0, hessian=FALSE)
    
    npar        <- length(mod$estimate)
    Hessian.eps <- 1e-7
    eps.mtx     <- diag(rep(Hessian.eps, npar))
    grad.at.max <- mod$gradient
    ObsInfo.tmp <- ObsInfo <- matrix(NA, npar, npar)
    
    ## Observed Information
    for (j in 1:npar){
        temp            <- logLikeCalc4(mod$estimate+eps.mtx[j,], y=Y, x=XMat, id=id, ProfileCol = Prof)
        ObsInfo.tmp[j,] <- (attr(temp,"gradient")-grad.at.max)/(Hessian.eps)
    }
    for (m in 1:npar){
        for (n in 1:npar){ ObsInfo[m,n] <-  (ObsInfo.tmp[m,n]+ObsInfo.tmp[n,m])/2}}
    Est <- mod$estimate[-Prof]
    Covar <- solve(ObsInfo[-Prof,-Prof])
    
    date()
    
    ########################################### 
    TransIndMtx     = matrix(1,3,4)
    TransIndMtx[,1] = 0
    Prof2           = c(7,10,13)
    mod2 <- nlm(logLikeCalc4, mod$estimate,  yval=Y2, x=XMat2, id=id2, UseGrad=TRUE, stepmax=1, iterlim=250, 
                check.analyticals = FALSE, print.level=0, hessian=FALSE, TransIndMtx=TransIndMtx, ProfileCol=Prof2)
    
    npar        <- length(mod2$estimate)
    Hessian.eps <- 1e-7
    eps.mtx     <- diag(rep(Hessian.eps, npar))
    grad.at.max <- mod2$gradient
    ObsInfo.tmp <- ObsInfo <- matrix(NA, npar, npar)
    
    ## Observed Information
    for (j in 1:npar){
        temp            <- logLikeCalc4(mod2$estimate+eps.mtx[j,], y=Y2, x=XMat2, id=id2, ref.muc=3, ProfileCol = Prof2)
        ObsInfo.tmp[j,] <- (attr(temp,"gradient")-grad.at.max)/(Hessian.eps)
    }
    for (m in 1:npar){
        for (n in 1:npar){ ObsInfo[m,n] <-  (ObsInfo.tmp[m,n]+ObsInfo.tmp[n,m])/2}}
    Est2 <- mod2$estimate[-Prof2]
    Covar2 <- solve(ObsInfo[-Prof2,-Prof2])
    
    est[[i]] <- Est
    covar[[i]] <- Covar
    conv[[i]] <- mod$code
    est2[[i]] <- Est2
    covar2[[i]] <- Covar2
    conv2[[i]] <- mod2$code
    
}

Results <- list(est=est, covar=covar, conv=conv,
                est2=est2, covar2=covar2, conv2=conv2)
save(Results, file=paste('output/i',iter,'.Rdata', sep=""))
}
