
# Uncomment to test random gammas
delta_it_r <- function(mm,            # row of marginal (mean) probabilities of each response category k=1...K 
                       mm.lag,        # lagged value of mm
                       gamma.mat,     # Transition log odds ratios
                       K=length(mm),
                       ...)           # K
{
    K1          <- K-1
    mm          <- mm[-K]
    Delta.vec   <- rep(0, K1)
    del         <- rep(1, K1)
    while (max(abs(del))>1e-4)
    {
        Delta.mat   <- matrix(rep(Delta.vec,each=K), ncol=K, byrow=TRUE)
        hmat.num    <- exp(Delta.mat+gamma.mat)
        hmat.denom  <- 1+ colSums(hmat.num)
        hmat        <- sweep(hmat.num,2,hmat.denom,"/")
        fDelta      <- hmat %*% mm.lag - mm
        df.dDelta   <- matrix(0, K1, K1)
        
        ## diagonal for df.dDelta
        for (l in 1:K1)
          df.dDelta[l,l] <- sum(hmat[l,]*(1-hmat[l,])*mm.lag)

        ## upper triangle for df.dDelta
        for (l in 1:(K1-1))
        {
          for (m in (l+1):K1)
          {
            df.dDelta[l,m] <- -sum(hmat[l,]*hmat[m,]*mm.lag)
          }
        }
        ## upper triangle for df.dDelta
        df.dDelta[lower.tri(df.dDelta)] <-  t(df.dDelta)[lower.tri(df.dDelta)]

        del        <- solve(df.dDelta) %*% fDelta
        Delta.vec  <- Delta.vec - del
    }
    Delta.vec
}

delta_it <- delta_it_r