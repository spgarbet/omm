####################################################
####################################################
####################################################################################
####################################################################################
####################################################################################
### Solving for the (K-1)-vector of Delta_{itk}... at time t
findDeltait <- function(mm,           # row of marginal (mean) probabilities of each response category k=1...K 
                        mm.lag,       # lagged value of mm
                        gamma.mat,    # Transition log odds ratios
                        nlevels){     # K
    nlevelsmin1 <- nlevels-1
    mm          <- mm[-nlevels]
    Delta.vec   <- rep(0, nlevelsmin1)
    del         <- rep(1,nlevelsmin1)
    while (max(del)>1e-4){
        
        Delta.mat       <- matrix(rep(Delta.vec,each=nlevels), ncol=nlevels, byrow=TRUE)
        hmat.num        <- exp(Delta.mat+gamma.mat)
        hmat.denom      <- 1+ colSums(hmat.num)
        hmat            <- sweep(hmat.num,2,hmat.denom,"/")
        fDelta    <- hmat %*% mm.lag - mm
        df.dDelta <- matrix(0, nlevelsmin1, nlevelsmin1)
        for (i in 1:nlevelsmin1){
            tmp           <- matrix(0, nlevelsmin1, nlevels)
            tmp[i,]       <- 1
            new.hmat      <- tmp-hmat
            new.hmat      <- new.hmat*tmp ## this differs from LD2007, but I think it is right
            dh.dDelta     <- new.hmat*hmat
            df.dDelta[i,] <- dh.dDelta %*% mm.lag
        }
        del <- solve(df.dDelta) %*% fDelta
        Delta.vec  <- Delta.vec - del
    }
    Delta.vec
}
####################################################################################
####################################################################################
####################################################################################

expit <- function(x){expx <- exp(x)
expx/(1+expx)}

N           <- 1000
nlevels     <- 5 
nlevelsmin1 <- nlevels-1
mi          <- 11  
alpha       <- c(-2,  -1,  0, 1.5)
beta        <- c(-0.25, -0.1, -0.1)

id          <- rep(1:N, each=mi)
tx          <- rep(rep(c(0,1), each=mi), N/2)
t           <- rep(seq(0,2,length.out=mi), N)

## This creates a transition log odds ratio matrix that comes close to 
## forcing an absorbing state (K=1) but it is not deterministic.
## I am not sure it makes sense in a setting of a deterministic 
## absorbing state.  
gamma.mat <- rbind( c( 10,  0,   0,   0,   0),
                    c(-25,  1.5, 0,   0,   0),
                    c(-25,  0,   1.5, 0,   0),
                    c(-25,  0,   0,   1.5, 0))

date()
## linear predictor without the intercepts X^t_it Beta, Sec 1.1
lp          <- cbind(tx,t,tx*t) %*% beta
y <- yval <- NULL

for (i in unique(id)){
    idi <- id[id==i]       
    lpi <- lp[id==i]
    txi <- tx[id==i]
    ti  <- t[id==i]

    ## marginal mean / probability matrix and lagged values
    ## rows correspond to time
    ## columns correspond to outcome category
    lpi.mat <- matrix(NA, mi, nlevelsmin1)
    for (j in 1:mi){ for( k in 1:nlevelsmin1){ lpi.mat[j,k] <- lpi[j]+alpha[k] }}

    cprobi.mat <- cbind(0, expit(lpi.mat), 1) ## cumulative prob matrix
    probi.mat  <- matrix(NA, mi, nlevels) ## multinomial probabilities
    for (k in 2:(nlevels+1)){
        probi.mat[,(k-1)] <- cprobi.mat[,k]-cprobi.mat[,(k-1)]
    }

    ### equal probability for all categories but none in K=1 and K=K at baseline
    ### I'm not sure if this value has any impact at all, based on the way I am generating
    ### the data.  I need to doublecheck.
    probi.matlag <- rbind(c(0, rep(1/(nlevels-2),nlevels-2), 0), probi.mat[-mi,])

    ## Calculate Deltait across all timepoints.  
    Deltait <- NULL
    for (j in 1:mi)
    {
      Deltait <- rbind(Deltait, c(findDeltait(mm=probi.mat[j,], mm.lag=probi.matlag[j,], gamma.mat=gamma.mat, nlevels=nlevels)))
    }

    # use marginal probs to simulate Yi1
    yi <- yilag <- c(rmultinom(1,1,probi.mat[1,]))
    
    ## sim Yi2... Yimi
    for (j in 2:mi){
        tmp      <- exp(Deltait[j,] + c(gamma.mat %*% yilag))
        tmp1     <- 1+sum(tmp)
        prob.y.c <- c(tmp,1)/tmp1 
        yilag    <- c(rmultinom(1,1,prob.y.c))
        yi       <- rbind(yi, yilag)
    }
    yival <- rep(10, nrow(yi))
    for(j in 1:nrow(yi)){ yival[j] <- which(yi[j,]==1) }
    y <- rbind(y, yi)
    yval <- c(yval, yival)
}
date()    

## some testing.
 
library(rms) 
## since the mean and dependence model parameters are orthogonal, Franks function
## should repreduce the mean model
TestDat <- data.frame(id, tx, t, yval)
TestDat <- subset(TestDat, subset=t>=0)
mod <- orm(yval~tx*t, data=TestDat)
    

####################################################################################
####################################################################################
####################################################################################
#
# Create simple test case to optimize on

j <- 10
mm <- probi.mat[j,]
mm.lag=probi.matlag[j,]
gamma.mat
nlevels

save(mm, mm.lag, gamma.mat, nlevels, file="delta_it-input.Rdata")

