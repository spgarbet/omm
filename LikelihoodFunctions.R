
####################################################################################
### Solving for the (K-1)-vector of Delta_{itk}... at time t
findDeltait <- function(mm,       # row of marginal (mean) probabilities of each response category k=1...K 
                        mm.lag,       # lagged value of mm
                        gamma.mat,    # Transition log odds ratios
                        K){     # K
    K1 <- K-1
    mm          <- mm[-K]
    Delta.vec   <- rep(0, K1)
    del         <- rep(1, K1)
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
        ## lower triangle for df.dDelta
        df.dDelta[lower.tri(df.dDelta)] <-  t(df.dDelta)[lower.tri(df.dDelta)]
        #print(df.dDelta)
        del <- solve(df.dDelta) %*% fDelta
        Delta.vec  <- Delta.vec - del
    }
    Delta.vec
}


####################################################################################
####################################################################################
####################################################################################
expit <- plogis

####################################################################################
####################################################################################
####################################################################################

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
        for (j in 1:mi){ Deltait <- rbind(Deltait, c(Delta_calc(mm=probi.mat[j,], mm.lag=probi.matlag[j,], gamma.mat=gamma.mat0)))}
        
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


Calc.TransProbs <- function(Y, id){
    
    L <- length(Y)
    Y1 <- c(NA,Y[-L])
    dup <- duplicated(id)
    Y <- Y[dup]
    Y1 <- Y1[dup]
    
    tab <- table(Y,Y1)
    joint <- tab/L
    trans <- sweep(tab, 2, colSums(tab), "/")
 
    out <- list(Tx.Mtx.Col2Row=trans, Joint.Mtx=joint, Marg.Prob=table(Y)/L)
    out
}

####################################################################################
####################################################################################
####################################################################################

logLikeCalc <- function(params, yval,x, id){
    
    States  <- sort(unique(yval))
    K       <- length(States)
    K1      <- K-1
    uid     <- unique(id)
    N       <- length(uid)
    npar    <- length(params)
    alpha.e <- params[1:K1]
    beta.e  <- params[(K1+1):(npar-(K1^2))]
    gamma.e <- matrix(params[(npar+1-(K1^2)):npar], K1, K1, byrow=TRUE)
    gamma.mtx <- cbind(gamma.e, 0)
    li1 <- li2 <- rep(0,N)
    #print(params)
    for (i in 1:N){
     #   print(uid[i])
        yival <- yval[id==uid[i]]
        xi    <- x[id==uid[i],]
        mi    <- nrow(xi)
      #  print(mi)
        ZiMat   <- YiMat <- cprobi.m <- probi.m <- matrix(NA, mi, K)
        logLi1  <- logLi2 <- 0
        logLij2 <- rep(0, mi)
        
        for (k in 1:K){ ZiMat[,k] <- as.integer(yival<= States[k])
        YiMat[,k] <- as.integer(yival== States[k])
        }
        
        for (k in 1:K1){cprobi.m[,k] <- expit(alpha.e[k] + xi %*% beta.e)}
        cprobi.m[,K] <- 1
        probi.m[,1] <- cprobi.m[,1]
        for (k in 2:K){probi.m[,k] <- cprobi.m[,(k)]- cprobi.m[,(k-1)]}
        
        ## LogLi1 only comes from the marginal portion of the model (first observation)
        for (k in 1:K1){ #print(cprobi.m[1,]);print(probi.m[1,])
            logLi1 <- logLi1 +  ZiMat[1,k]    *log(cprobi.m[1,k]     / probi.m[1,(k+1)]) - 
                ZiMat[1,(k+1)]*log(cprobi.m[1,(k+1)] / probi.m[1,(k+1)])}
        
        ## logLi2 comes from observations 2 to mi
        YiMat2 <- YiMat[,-K]
        
        for (j in 2:mi){
            #Deltaij2       <- findDeltait(mm=probi.m[j,], mm.lag=probi.m[(j-1),], gamma.mat=gamma.mtx, K=K)
            Deltaij2       <- Delta_calc(mm=probi.m[j,], mm.lag=probi.m[(j-1),], gamma.mat=gamma.mtx)
            
            #lpij.c         <- Deltaij2 + t(gamma.e) %*% YiMat2[(j-1),]
            lpij.c         <- Deltaij2 + gamma.mtx[,yival[(j-1)]] #%*% YiMat2[(j-1),]
            
            probij.cK      <- 1/(1+sum(exp(lpij.c)))
            logLij2[j] <- YiMat2[j,] %*% lpij.c + log(probij.cK)
            #logLij2[j]     <- lpij.c[yival[j]] + log(probij.cK)
            
        }
        li1[i] <- logLi1
        li2[i] <- sum(logLij2)
    }
    #print(sum(li1))
    #print(sum(li2))
    
    li <- -1*(sum(li1)+sum(li2)) 
    li
}

logLikeScoreCalc <- function(params, yval,x, id, ProfileCol=NA){
    
    
    logLike <- logLikeCalc(params, yval,x, id)
    
    npar <- length(params)
    param <- params
    Deriv <- sapply(1:npar,  function(rr)
    {
        grad(function(mmm) { new.param <- param
        new.param[rr] <- mmm
        logLikeCalc(new.param, yval, x, id)
        },
        param[rr],
        method="Richardson")
    }
    )
    Deriv[ProfileCol] <- 0
    out <- logLike
    attr(out,"gradient") <- Deriv
    out
}



ScoreOnlyCalc <- function(params, yval,x, id, ProfileCol=NA){
    
    
    npar <- length(params)
    param <- params
    Deriv <- sapply(1:npar,  function(rr)
    {
        grad(function(mmm) { new.param <- param
        new.param[rr] <- mmm
        logLikeCalc(new.param, yval, x, id)
        },
        param[rr],
        method="simple")
    }
    )
    Deriv[ProfileCol] <- 0
    Deriv
}

logLikeCalc3 <- function(params, yval,x, id, ProfileCol=NA){
    
    #params <- params
    #yval=Y
    #x=XMat
    #id=id
    
    States    <- sort(unique(yval))
    K         <- length(States)
    K1        <- K-1
    uid       <- unique(id)
    N         <- length(uid)
    Ntot      <- length(Y)
    npar      <- length(params)
    alpha.e   <- params[1:K1]
    beta.e    <- params[(K1+1):(npar-(K1^2))]
    gamma.e   <- matrix(params[(npar+1-(K1^2)):npar], K1, K1, byrow=TRUE)
    gamma.mtx <- cbind(gamma.e, 0)
    
    alpha.len <- K1
    beta.len  <- length(beta.e)
    gamma.len <- K1^2
    
    ## Z and Y matrices
    ZMat <- YMat <- matrix(0, Ntot, K)
    for (k in 1:K){ YMat[,k] <- ifelse(yval==k, 1, 0)
                    ZMat[,k] <- ifelse(yval<=k, 1, 0)
    }
    
    ## marginal cumulative probs and state probs
    cprob.m <- prob.m <- matrix(NA, Ntot, K)
    for (k in 1:K1){
        cprob.m[,k] <- expit(alpha.e[k] + x %*% beta.e)
    }
    cprob.m[,K] <- 1
    prob.m[,1]  <- cprob.m[,1]
    for (k in 2:K){
        prob.m[,k] <- cprob.m[,(k)]- cprob.m[,(k-1)]
    }
    one.min.cprob.m <- 1-cprob.m
    one.min.prob.m  <- 1-prob.m
    
    ## Phi (we only use the first observation for each subject)
    phi.k <- matrix(NA, Ntot, K1)
    for (k in 1:K1){
        phi.k[,k] <- log(cprob.m[,k]/(prob.m[,(k+1)]))
    }
    g.phi.k <- log(1+exp(phi.k))
    exp.g.phi.k <- exp(g.phi.k)
    
    ## dPhi/deta (we only use the first observation for each subject)
    ## dpi/deta
    dphi.k.deta.k <- dphi.k.deta.k1 <- dpi.k.deta.k <- matrix(NA, Ntot, K1)
    for (k in 1:K1){ 
        dphi.k.deta.k[,k]  <-  exp.g.phi.k[,k]*one.min.cprob.m[,k]  
        dphi.k.deta.k1[,k] <-  -exp.g.phi.k[,k]*one.min.cprob.m[,(k+1)]
        dpi.k.deta.k[,k]   <-  cprob.m[,k]*one.min.cprob.m[,k]
    }
    
    ## Some useful matrices to be used later        
    tmp1 <- diag(1,K1)
    tmp1a <- matrix(0, K1, K1)
    tmp2 <- rbind(tmp1[-1,], rep(0,K1))
    #tmp2 <- rbind(rep(0,K1),
    #             tmp1[-K1,])
    tmp4 <- matrix(0, K1, K1^2)
    
    ## To be used in the gradient calculation
    tmp.mat <- matrix(0, K1, K1^2)
    for (ppp in 1:K1){ 
        tmp.mat[ppp,(c(1:K1)+(ppp-1)*K1)] <- 1 
    }
    
    li1 <- li2 <- rep(0,N)
    Grad1Vec <- Grad2Vec <- rep(0,npar)
    Grad1Mat <- Grad2Mat <- matrix(0,N, npar)
    
    for (i in 1:N){
        
        yival            <- yval[id==uid[i]]
        xi               <- x[id==uid[i],]
        mi               <- nrow(xi)
        ZiMat            <- ZMat[id==uid[i],]
        YiMat            <- YMat[id==uid[i],]
        cprobi.m         <- cprob.m[id==uid[i],]
        probi.m          <- prob.m[id==uid[i],]
        one.min.cprobi.m <- one.min.cprob.m[id==uid[i],]
        phii.k           <- phi.k[id==uid[i],]
        g.phii.k         <- g.phi.k[id==uid[i],]
        dphii.k.deta.k   <- dphi.k.deta.k[id==uid[i],]
        dphii.k.deta.k1  <- dphi.k.deta.k1[id==uid[i],]
        dpii.k.deta.k    <- dpi.k.deta.k[id==uid[i],]
        
        logLi1  <- logLi2 <- 0
        logLij2 <- rep(0, mi)
        
        ## LogLi1 and dLogLi1/d(alpha,beta)
        ## LogLi1 only comes from the marginal portion of the model (first observation)
        for (k in 1:K1){ #print(cprobi.m[1,]);print(probi.m[1,])
            logLi1 <- logLi1 +  ZiMat[1,k]*phii.k[1,k] - ZiMat[1,(k+1)]*g.phii.k[1,k]
        }
        
        ## LogLi1/dtheta        
        tmp3           <- matrix(rep(xi[1,],K1), nrow=K1, byrow=TRUE)              
        deta.k.dtheta  <- cbind(tmp1,tmp3,tmp4)
        deta.k1.dtheta <- cbind(tmp2,tmp3,tmp4)
        
        Grad1iVec <- Grad2iVec <- rep(0, npar)
        
        for (k in 1:K1){ 
            Grad1iVec <- Grad1iVec+(ZiMat[1,k]-(ZiMat[1,(k+1)]*cprobi.m[1,k]/cprobi.m[1,(k+1)]))*(dphii.k.deta.k[1,k]*deta.k.dtheta[k,] + dphii.k.deta.k1[1,k]*deta.k1.dtheta[k,])
        }   
        
        ## logLi2 comes from observations 2 to mi
        YiMat2 <- YiMat[,-K]
        #print(i)
        #print(probi.m[1,])
        #print(probi.m[(j-1),])
        
        ##############################
        ##############################
        ##############################
        
        for (j in 2:mi){
            mm           <- probi.m[j,]
            mm.lag       <- probi.m[(j-1),]
            Yit          <- YiMat2[j,]
            Yit1         <- YiMat2[(j-1),]
            xit          <- xi[j,]
            xit1         <- xi[(j-1),]
            dpi.deta     <- dpii.k.deta.k[j,]
            dpi.deta.lag <- dpii.k.deta.k[(j-1),]
            #print(mm)
            #print(mm.lag)
            #print(c(gamma.mtx))
            #Deltaij2   <- findDeltait(mm=probi.m[j,], mm.lag=probi.m[(j-1),], gamma.mat=gamma.mtx, K=K)
            Deltaij2    <- Delta_calc(mm=mm, mm.lag=mm.lag, gamma.mat=gamma.mtx)
            #lpij.c     <- Deltaij2 + t(gamma.e) %*% YiMat2[(j-1),]
            lpij.c      <- Deltaij2 + gamma.mtx[,yival[(j-1)]] #%*% YiMat2[(j-1),]
            probij.cK   <- 1/(1+sum(exp(lpij.c)))
            logLij2[j]  <- YiMat2[j,] %*% lpij.c + log(probij.cK)
            #logLij2[j] <- lpij.c[yival[j]] + log(probij.cK)
            muit.c      <- exp(lpij.c)*probij.cK
            
            tmp3   <- matrix(rep(xit,K1), nrow=K1, byrow=TRUE)              
            tmp3l  <- matrix(rep(xit1,K1), nrow=K1, byrow=TRUE)              
            
            deta.k.dtheta      <- cbind(tmp1,tmp3)
            deta.k.dtheta.lag  <- cbind(tmp1,tmp3l)
            
            Delta.mat   <- matrix(rep(Deltaij2,each=K), ncol=K, byrow=TRUE)
            hmat.num    <- exp(Delta.mat+gamma.mtx)
            hmat.denom  <- 1+ colSums(hmat.num)
            hmat        <- sweep(hmat.num,2,hmat.denom,"/")
            
        
            df.dDelta   <- matrix(0, K1, K1)
            ## Left side of system of equations that solve for dDelta/dtheta
            for (l in 1:K1){df.dDelta[l,l] <- sum(hmat[l,]*(1-hmat[l,])*mm.lag)}
            ## upper triangle for df.dDelta
            for (l in 1:(K1-1)){
                for (m in (l+1):K1){
                    df.dDelta[l,m] <- -sum(hmat[l,]*hmat[m,]*mm.lag)
                }}
            ## lower triangle for 
            df.dDelta[lower.tri(df.dDelta)] <-  t(df.dDelta)[lower.tri(df.dDelta)]
            
            ## right side of system of equations that solve for dDelta/dtheta (only alpha and beta)
            
            dpi.k.dtheta  <- sweep(deta.k.dtheta, 1,dpi.deta,"*")
            dpi.k1.dtheta <- rbind(0,dpi.k.dtheta[-K1,])
            dmum.k.dtheta <- dpi.k.dtheta - dpi.k1.dtheta
            
            tmp.dpi           <- sweep(deta.k.dtheta.lag, 1, dpi.deta.lag,"*")
            dpi.k.dtheta.lag  <- rbind(tmp.dpi, 0)
            dpi.k1.dtheta.lag <- rbind(0,tmp.dpi)
            dmum.k.dtheta.lag <- dpi.k.dtheta.lag - dpi.k1.dtheta.lag
            h.dmum.lag.dtheta <- hmat %*% dmum.k.dtheta.lag
        
            rhs.alpha.beta    <- dmum.k.dtheta - h.dmum.lag.dtheta
            
            ## right side of system of equations that solve for dDelta/dtheta (only gamma)
            hmat1        <- hmat[,-K]
            mm.lag1      <- mm.lag[-K]
            hmat2        <- matrix(rep(t(hmat1),K1), nrow=K1, byrow=TRUE)
            mm.lag.mat   <- matrix(rep(mm.lag1, K1^2), nrow=K1, byrow=TRUE)
            hmat1.K1times <- matrix(rep(hmat1, K1), nrow=K1)
            
            hmat3            <- tmp.mat-hmat1.K1times
            dh.dgamma.mm.lag <- -1*hmat2*mm.lag.mat*hmat3
            
            rhs <- cbind(rhs.alpha.beta, dh.dgamma.mm.lag)
            ddelta.dtheta <- matrix(0,K1, npar)
            for (mmm in 1:npar){
                ddelta.dtheta[,mmm] <- solve(df.dDelta, rhs[,mmm])
            }
            #ddelta.dtheta1 <- ddelta.dtheta 
            Grad2iVec <- Grad2iVec + (Yit-muit.c) %*% ddelta.dtheta + c(rep(0, alpha.len+beta.len), c(rep(Yit-muit.c, each=K1) * rep(Yit1, K1)))
        }
        li1[i]       <- -1*logLi1
        li2[i]       <- -1*sum(logLij2)
        Grad1Mat[i,] <- -1*Grad1iVec
        Grad2Mat[i,] <- -1*Grad2iVec
        
    }
    #print(sum(li1))
    #print(sum(li2))
    
    Grad1Vec <- apply(Grad1Mat,2,sum)
    Grad2Vec <- apply(Grad2Mat,2,sum)
    GradVec  <- Grad1Vec+Grad2Vec 
    GradVec[ProfileCol] <- 0
    #print(Grad1Vec)
    #print(Grad2Vec)
    li       <- sum(li1)+sum(li2)
    out      <- sum(li)
    attr(out,"gradient") <- GradVec
    out
}


logLikeCalc4 <- function(params, yval, x, id, ProfileCol=NA, ref.muc=NA, UseGrad=TRUE){
    
  #  params <- inits
  #  yval=Y
  #  x=XMat
  #  id=id
  # ref.muc <- NA
  #   # 
    States    <- sort(unique(yval))
    K         <- length(States)
    K1        <- K-1
    uid       <- unique(id)
    N         <- length(uid)
    Ntot      <- length(Y)
    npar      <- length(params)
    alpha.e   <- params[1:K1]
    beta.e    <- params[(K1+1):(npar-(K1^2))]
    gamma.e   <- matrix(params[(npar+1-(K1^2)):npar], K1, K1, byrow=TRUE)
    gamma.mtx <- cbind(gamma.e, 0)
    
    alpha.len <- K1
    beta.len  <- length(beta.e)
    gamma.len <- K1^2
    
    ## Z and Y matrices
    ZMat <- YMat <- matrix(0, Ntot, K)
    for (k in 1:K){ YMat[,k] <- ifelse(yval==k, 1, 0)
    ZMat[,k] <- ifelse(yval<=k, 1, 0)
    }
    
    ## marginal cumulative probs and state probs
    cprob.m <- prob.m <- matrix(NA, Ntot, K)
    for (k in 1:K1){
        cprob.m[,k] <- expit(alpha.e[k] + x %*% beta.e)
    }
    cprob.m[,K] <- 1
    prob.m[,1]  <- cprob.m[,1]
    for (k in 2:K){
        prob.m[,k] <- cprob.m[,(k)]- cprob.m[,(k-1)]
    }
    one.min.cprob.m <- 1-cprob.m
    one.min.prob.m  <- 1-prob.m
    
    ## Phi (we only use the first observation for each subject)
    phi.k <- matrix(NA, Ntot, K1)
    for (k in 1:K1){
        phi.k[,k] <- log(cprob.m[,k]/(prob.m[,(k+1)]))
    }
    g.phi.k <- log(1+exp(phi.k))
    exp.g.phi.k <- exp(g.phi.k)
    
    ## dPhi/deta (we only use the first observation for each subject)
    ## dpi/deta
    dphi.k.deta.k <- dphi.k.deta.k1 <- dpi.k.deta.k <- matrix(NA, Ntot, K1)
    for (k in 1:K1){ 
        dphi.k.deta.k[,k]  <-  exp.g.phi.k[,k]*one.min.cprob.m[,k]  
        dphi.k.deta.k1[,k] <-  -exp.g.phi.k[,k]*one.min.cprob.m[,(k+1)]
        dpi.k.deta.k[,k]   <-  cprob.m[,k]*one.min.cprob.m[,k]
    }
    
    ## Some useful matrices to be used in gradient calculations       
    tmp1 <- diag(1,K1)
    tmp1a <- matrix(0, K1, K1)
    tmp2 <- rbind(tmp1[-1,], rep(0,K1))
    #tmp2a <- rbind(rep(0,K1),
     #            tmp1[-K1,])
    tmp4 <- matrix(0, K1, K1^2)
    
    ## To be used in the gradient calculation
    tmp.mat <- matrix(0, K1, K1^2)
    for (ppp in 1:K1){ 
        tmp.mat[ppp,(c(1:K1)+(ppp-1)*K1)] <- 1 
    }
    
    li1 <- li2 <- rep(0,N)
    Grad1Vec <- Grad2Vec <- rep(0,npar)
    Grad1Mat <- Grad2Mat <- matrix(0,N, npar)
    
    for (i in 1:N){
        
        yival            <- yval[id==uid[i]]
        xi               <- x[id==uid[i],]
        mi               <- nrow(xi)
        ZiMat            <- ZMat[id==uid[i],]
        YiMat            <- YMat[id==uid[i],]
        cprobi.m         <- cprob.m[id==uid[i],]
        probi.m          <- prob.m[id==uid[i],]
        one.min.cprobi.m <- one.min.cprob.m[id==uid[i],]
        phii.k           <- phi.k[id==uid[i],]
        g.phii.k         <- g.phi.k[id==uid[i],]
        dphii.k.deta.k   <- dphi.k.deta.k[id==uid[i],]
        dphii.k.deta.k1  <- dphi.k.deta.k1[id==uid[i],]
        dpii.k.deta.k    <- dpi.k.deta.k[id==uid[i],]
        
        logLi1  <- logLi2 <- 0
        logLij2 <- rep(0, mi)
        
        ## LogLi1 and dLogLi1/d(alpha,beta)
        ## LogLi1 only comes from the marginal portion of the model (first observation)
        for (k in 1:K1){ #print(cprobi.m[1,]);print(probi.m[1,])
            logLi1 <- logLi1 +  ZiMat[1,k]*phii.k[1,k] - ZiMat[1,(k+1)]*g.phii.k[1,k]
        }
        
        ## LogLi1/dtheta        
        tmp3           <- matrix(rep(xi[1,],K1), nrow=K1, byrow=TRUE)              
        deta.k.dtheta  <- cbind(tmp1,tmp3,tmp4)
        deta.k1.dtheta <- cbind(tmp2,tmp3,tmp4)
        
        Grad1iVec <- Grad2iVec <- rep(0, npar)
        
        for (k in 1:K1){ 
            Grad1iVec <- Grad1iVec+(ZiMat[1,k]-(ZiMat[1,(k+1)]*cprobi.m[1,k]/cprobi.m[1,(k+1)]))*(dphii.k.deta.k[1,k]*deta.k.dtheta[k,] + dphii.k.deta.k1[1,k]*deta.k1.dtheta[k,])
        }   
        
        ## logLi2 comes from observations 2 to mi
        YiMat2 <- YiMat[,-K]
        #print(i)
        #print(probi.m[1,])
        #print(probi.m[(j-1),])
        
        ##############################
        ##############################
        ##############################
        ############ Altered for new ref.muc (reference state in conditional model)
        if (!is.na(ref.muc)){
            probi.m <- cbind(probi.m[,-ref.muc], probi.m[,ref.muc])
            YiMat.tmp <- cbind(YiMat[,-ref.muc], YiMat[,ref.muc])
            YiMat2 <- YiMat.tmp[,-K] 
            yival2 <- ifelse(yival<ref.muc, yival,
                         ifelse(yival==ref.muc, K, yival-1))
            yival <- yival2
        }
        ##############################
        ##############################
        ##############################
        #print(i)
        
        for (j in 2:mi){
            mm           <- probi.m[j,]
            mm.lag       <- probi.m[(j-1),]
            Yit          <- YiMat2[j,]
            Yit1         <- YiMat2[(j-1),]
            xit          <- xi[j,]
            xit1         <- xi[(j-1),]
            dpi.deta     <- dpii.k.deta.k[j,]
            dpi.deta.lag <- dpii.k.deta.k[(j-1),]
            #print(mm)
            #print(mm.lag)
            #print(c(gamma.mtx))
            #Deltaij2   <- findDeltait(mm=probi.m[j,], mm.lag=probi.m[(j-1),], gamma.mat=gamma.mtx, K=K)
            Deltaij2    <- Delta_calc(mm=mm, mm.lag=mm.lag, gamma.mat=gamma.mtx)
            #lpij.c     <- Deltaij2 + t(gamma.e) %*% YiMat2[(j-1),]
            lpij.c      <- Deltaij2 + gamma.mtx[,yival[(j-1)]] #%*% YiMat2[(j-1),]
            probij.cK   <- 1/(1+sum(exp(lpij.c)))
            logLij2[j]  <- YiMat2[j,] %*% lpij.c + log(probij.cK)
            #logLij2[j] <- lpij.c[yival[j]] + log(probij.cK)
            muit.c      <- exp(lpij.c)*probij.cK
            
            # Delta.mat   <- matrix(rep(Deltaij2,each=K), ncol=K, byrow=TRUE)
            # hmat.num    <- exp(Delta.mat+gamma.mtx)
            # hmat.denom  <- 1+ colSums(hmat.num)
            # hmat        <- sweep(hmat.num,2,hmat.denom,"/")
            hmat <- hmat_calc(Deltaij2, gamma.mtx)

            # df.dDelta   <- matrix(0, K1, K1)
            # ## Left side of system of equations that solve for dDelta/dtheta
            # for (l in 1:K1){df.dDelta[l,l] <- sum(hmat[l,]*(1-hmat[l,])*mm.lag)}
            # ## upper triangle for df.dDelta
            # for (l in 1:(K1-1)){
            #     for (m in (l+1):K1){
            #         df.dDelta[l,m] <- -sum(hmat[l,]*hmat[m,]*mm.lag)
            #     }}
            # ## lower triangle for
            # df.dDelta[lower.tri(df.dDelta)] <-  t(df.dDelta)[lower.tri(df.dDelta)]
            df.dDelta <- dfdDelta_calc(mm.lag=mm.lag, hmat=hmat)
            ## right side of system of equations that solve for dDelta/dtheta (only alpha and beta)

            #dpi.k.dtheta  <- sweep(deta.k.dtheta, 1,dpi.deta,"*")
            #dpi.k1.dtheta <- rbind(0,dpi.k.dtheta[-K1,])
            #dmum.k.dtheta <- dpi.k.dtheta - dpi.k1.dtheta
            
            #tmp3   <- matrix(rep(xit,K1), nrow=K1, byrow=TRUE)
            #tmp3l  <- matrix(rep(xit1,K1), nrow=K1, byrow=TRUE)
            
            #deta.k.dtheta      <- cbind(tmp1,tmp3)
            #deta.k.dtheta.lag  <- cbind(tmp1,tmp3l)
            
            #tmp.dpi.lag       <- dpidtheta_calc(deta.k.dtheta.lag, dpi.deta.lag)
            tmp.dpi.lag      <- dpidtheta_calc1(tmp1, dpi.deta.lag, xit1)
            #tmp.dpi.lag       <- sweep(deta.k.dtheta.lag, 1, dpi.deta.lag,"*")
            dpi.k.dtheta.lag  <- rbind(tmp.dpi.lag, 0)
            dpi.k1.dtheta.lag <- rbind(0,tmp.dpi.lag)
            dmum.k.dtheta.lag <- dpi.k.dtheta.lag - dpi.k1.dtheta.lag
            
            ########################################################################
            ########################################################################
            ## for new ref.muc #############################
            
            #tmp.dpi           <- dpidtheta_calc(deta.k.dtheta, dpi.deta)
            tmp.dpi      <- dpidtheta_calc1(tmp1, dpi.deta, xit)
            #tmp.dpi           <- sweep(deta.k.dtheta, 1,dpi.deta,"*")
            dpi.k.dtheta      <- rbind(tmp.dpi, 0)
            dpi.k1.dtheta     <- rbind(0,tmp.dpi)
            dmum.k.dtheta.tmp <- dpi.k.dtheta - dpi.k1.dtheta
            if (is.na(ref.muc)){dmum.k.dtheta  <- dmum.k.dtheta.tmp[-K,]
            }else{             dmum.k.dtheta  <- dmum.k.dtheta.tmp[-ref.muc,]}
            
            if (!is.na(ref.muc)) dmum.k.dtheta.lag <- rbind(dmum.k.dtheta.lag[-ref.muc,], dmum.k.dtheta.lag[ref.muc,])
            h.dmum.lag.dtheta     <- hmat %*% dmum.k.dtheta.lag
            ########################################################################
            ########################################################################
            ########################################################################
            
            rhs.alpha.beta    <- dmum.k.dtheta - h.dmum.lag.dtheta

            ## right side of system of equations that solve for dDelta/dtheta (only gamma)
            hmat1        <- hmat[,-K]
            mm.lag1      <- mm.lag[-K]
            hmat2        <- matrix(rep(t(hmat1),K1), nrow=K1, byrow=TRUE)
            mm.lag.mat   <- matrix(rep(mm.lag1, K1^2), nrow=K1, byrow=TRUE)
            hmat1.K1times <- matrix(rep(hmat1, K1), nrow=K1)

            hmat3            <- tmp.mat-hmat1.K1times
            dh.dgamma.mm.lag <- -1*hmat2*mm.lag.mat*hmat3

            rhs <- cbind(rhs.alpha.beta, dh.dgamma.mm.lag)
            #ddelta.dtheta <- matrix(0,K1, npar)
            #for (mmm in 1:npar){
            #    ddelta.dtheta[,mmm] <- solve(df.dDelta, rhs[,mmm])
            #}
            ddelta.dtheta <- solve(df.dDelta, rhs)
            
            #ddelta.dtheta1 <- ddelta.dtheta
            Grad2iVec <- Grad2iVec + (Yit-muit.c) %*% ddelta.dtheta + c(rep(0, alpha.len+beta.len), c(rep(Yit-muit.c, each=K1) * rep(Yit1, K1)))
        }
        li1[i]       <- -1*logLi1
        li2[i]       <- -1*sum(logLij2)
        Grad1Mat[i,] <- -1*Grad1iVec
        Grad2Mat[i,] <- -1*Grad2iVec

    }

    Grad1Vec <- apply(Grad1Mat,2,sum)
    Grad2Vec <- apply(Grad2Mat,2,sum)
    GradVec  <- Grad1Vec+Grad2Vec
    GradVec[ProfileCol] <- 0

    li       <- sum(li1)+sum(li2)
    out      <- sum(li)
    if (UseGrad==TRUE) attr(out,"gradient") <- GradVec
    out
}

cluster.summary <- function( id, x, fun ){
    xlist <- split( x, id )
    nj <- unlist( lapply( xlist, length ) )
    xj <- unlist( lapply( xlist, fun) )
    xsummary <- rep( xj, nj )
    xsummary}

