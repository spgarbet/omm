path <- "~/rsch/MMOrdinal/Sims/TrueModels/"
setwd(path)
source(paste(path, "Functions.R", sep=""))
library(rms)
N         <- 200
mi        <- 30
id        <- rep(1:N, each=mi)
tx        <- rep(rep(c(0,1), each=mi), N/2)
t         <- rep(seq(0,(mi-1),length.out=mi), N)
XMat      <- cbind(tx,t,tx*t)
alpha     <-c(-2,-1,1)
beta      <- c(0.25,0.05,.025)
gamma.mat <- rbind( c(3,2,1),
                    c(1,2,1),
                    c(0,1,2))
params    <- c(alpha, beta, c(t(gamma.mat)))
ests.abs  <- hess.abs <- convergence.abs <- list()
ests.drop <- hess.drop <- convergence.drop <- list()

date()
Y           <- GenDatOMTM1(id=id,
                    XMat=XMat, 
                    alpha=alpha, 
                    beta=beta, 
                    gamma.mat=gamma.mat)
startdt <- date()
model.abs  <- optim( params, logLikeCalc, yval=Y, x=XMat, id=id, control=list(trace=2, maxit=10000), hessian=TRUE)
enddt <- date()

date()
lp        <- XMat %*% beta
gamma.mat0 <- cbind(gamma.mat, 0)
K1        <- length(alpha)
K         <- K1+1
start     <- 1

for (i in unique(id)){
    idi <- id[id==i]       
    lpi <- lp[id==i]
    txi <- tx[id==i]
    ti  <- t[id==i]
    mi  <- length(idi)
    
    
    ## marginal mean / probability matrix and lagged values
    ## rows correspond to time
    ## columns correspond to outcome category
    lpi.mat <- matrix(NA, mi, K1)
    for (j in 1:mi){ for( k in 1:K1){ lpi.mat[j,k] <- lpi[j]+alpha[k] }}
    
    cprobi.mat <- cbind(0, expit(lpi.mat), 1) ## cumuulative prob matrix
    probi.mat  <- matrix(NA, mi, K) ## multinomial probabilities
    for (k in 2:(K+1)){
        probi.mat[,(k-1)] <- cprobi.mat[,k]-cprobi.mat[,(k-1)]
    }
    
    ### This value has any impact based on the way I am generating the data.
    #probi.matlag <- rbind(c(0, rep(1/(K-2),K-2), 0), probi.mat[-mi,])
    #probi.matlag <- rbind(c(1/k, K), probi.mat[-mi,])
    probi.matlag <- rbind(probi.mat[1,], probi.mat[-mi,])
    
    Deltait <- NULL
    for (j in 1:mi){ Deltait <- rbind(Deltait, c(findDeltait(mm=probi.mat[j,], mm.lag=probi.matlag[j,], gamma.mat=gamma.mat, K=K)))}
    for (j in 1:mi){ Deltait <- rbind(Deltait, c(delta_it_c(mm=probi.mat[j,], mm.lag=probi.matlag[j,], gamma.mat=gamma.mat0)))}
    j <- 1
    findDeltait(mm=probi.mat[j,], mm.lag=probi.matlag[j,], gamma.mat=gamma.mat, K=K)
    nlevels=K
 
findDeltait <- function(mm,       # row of marginal (mean) probabilities of each response category k=1...K 
                        mm.lag,       # lagged value of mm
                        gamma.mat,    # Transition log odds ratios
                        K){     # K
    K1 <- K-1
    mm          <- mm[-K]
    Delta.vec   <- rep(0, K1)
    del         <- rep(1,K1)
    while (max(abs(del))>1e-4){
        
        Delta.mat   <- matrix(rep(Delta.vec,each=K), ncol=K, byrow=TRUE)
        hmat.num    <- exp(Delta.mat+gamma.mat)
        hmat.denom  <- 1+ colSums(hmat.num)
        hmat        <- sweep(hmat.num,2,hmat.denom,"/")
        fDelta      <- hmat %*% mm.lag - mm
        df.dDelta   <- matrix(0, K1, K1)
        
        ## diagonal for df.dDelta
        for (l in 1:K1){
            df.dDelta[l,l] <- sum(hmat[l,]*(1-hmat[l,])*mm.lag)
        }
        ## upper triangle for df.dDelta
        for (l in 1:(K1-1)){
        for (m in (l+1):K1){
            df.dDelta[l,m] <- -sum(hmat[l,]*hmat[m,]*mm.lag)
        }}
        ## upper triangle for df.dDelta
        df.dDelta[lower.tri(df.dDelta)] <-  df.dDelta[upper.tri(df.dDelta)]

        del <- solve(df.dDelta) %*% fDelta
        Delta.vec  <- Delta.vec - del
    }
    Delta.vec
}

















































## Force an Absorbing State
L          <- length(Y)
Ylag       <- c(10, Y[-L])
dup        <- duplicated(id)
Ylag[!dup] <- 10
Absorbed <- rep(0, L)
for (b in 2:L){ if (!is.na(Ylag[b]) & dup[b]==TRUE & (Ylag[b]==1 | Absorbed[(b-1)]==1)) 
                         Absorbed[b] <- 1 }

print(c("abs start", date()))
## Model keeping followup after absorbing state
YwithAbs   <- ifelse(Absorbed==0, Y, Absorbed)

Dat <- data.frame(id, YwithAbs, tx, t)
for(time in unique(t)){
    time <- 1
    mod <- orm(YwithAbs ~ tx, data=Dat, subset=t==time)
}



model.abs  <- optim( params, logLikeCalc, yval=YwithAbs, x=XMat, id=id, control=list(trace=0, maxit=10000), hessian=TRUE)

print(c("drop start", date()))
## Model dropping data after absorbing state
Absorbed1  <- ifelse(Absorbed==1 & t>.25, 1, 0)
id1        <- id[Absorbed1!=1]
XMat1      <- XMat[Absorbed1!=1,]
Y1         <- Y[Absorbed1!=1]
model.drop <- optim( params, logLikeCalc, yval=Y1, x=XMat1, id=id1, control=list(trace=0, maxit=10000), hessian=TRUE)
print(c("drop end", date()))

Calc.TransProbs(Y=Y, id=id, prob=T)
Calc.TransProbs(Y=Y1, id=id1, prob=T)
Calc.TransProbs(Y=YwithAbs, id=id, prob=T)

#table(Y, t)

ests.abs[[i]] <- model.abs$par
hess.abs[[i]] <- model.abs$hessian
convergence.abs[[i]] <- model.abs$convergence

ests.drop[[i]] <- model.drop$par
hess.drop[[i]] <- model.drop$hessian
convergence.drop[[i]] <- model.drop$convergence

}

Results <- list(ests.abs=ests.abs, hess.abs=hess.abs, convergence.abs=convergence.abs,
                ests.drop=ests.drop, hess.drop=hess.drop, convergence.drop=convergence.drop)
save(Results, file=paste(path,'/output/i',ITER,'.Rdata', sep=""))

#blah2 <- nlm( logLikeCalc, PARAMS, yval=Y, x=X, id=ID, print.level=2 )

