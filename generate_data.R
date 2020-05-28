
expit <- plogis

GenDatOMTM1 <- function(id, 
                        XMat, #### design matrix that should NNOT include any intercepts
                        alpha, 
                        beta, 
                        gamma.mat){
    
    lp        <- XMat %*% beta
    gamma.mat0 <- cbind(gamma.mat, 0)
    K1        <- length(alpha)
    K         <- K1+1
    start     <- 1
    
    ## linear predictor without the intercepts
    #y <- yval <- NULL
    yval <- rep(NA, length(lp))
    
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
        probi.matlag <- rbind(probi.mat[1,], probi.mat[-mi,])
        
        ## Calculate Deltait across all timepoints.  
        Deltait <- NULL
        #for (j in 1:mi){ Deltait <- rbind(Deltait, c(findDeltait(mm=probi.mat[j,], mm.lag=probi.matlag[j,], gamma.mat=gamma.mat0, K=K)))}
        for (j in 1:mi){ Deltait <- rbind(Deltait, c(delta_it(mm=probi.mat[j,], mm.lag=probi.matlag[j,], gamma.mat=gamma.mat0, trace=1)))}
        
        # use marginal probs to simulate Yi1
        
        yi <- yilag <- c(rmultinom(1,1,probi.mat[1,]))
        
        ## sim Yi2... Yimi
        for (j in 2:mi){
            tmp      <- exp(Deltait[j,] + c(gamma.mat0 %*% yilag))
            tmp1     <- 1+sum(tmp)
            prob.y.c <- c(tmp,1)/tmp1 
            yilag    <- c(rmultinom(1,1,prob.y.c))
            yi       <- rbind(yi, yilag)
        }
        yival <- rep(10, nrow(yi))
        for(j in 1:nrow(yi)){ yival[j] <- which(yi[j,]==1) }
        #y    <- rbind(y, yi)
        #yval <- c(yval, yival)
        yval[c(start:(start+mi-1))] <- yival
        start <- start+mi

    }
    yval
}


Calc.TransProbs <- function(Y, id, prob=T){
    L <- length(Y)
    tmp   <- split(Y, id)
    X     <- do.call("rbind", tmp)
    tt <- table( c(X[,-ncol(X)]), c(X[,-1]) )
    if(prob) tt <- tt / rowSums(tt)
    out <- list(Tx.Mtx= tt, Marg.Prob=table(Y)/L)
    out
}