
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
    
    lp         <- XMat %*% beta
    gamma.mat0 <- cbind(gamma.mat, 0)
    K1         <- length(alpha)
    K          <- K1+1
    start      <- 1
    
    ## linear predictor without the intercepts
    #y <- yval <- NULL
    yval <- rep(NA, length(lp))
    
    for (i in unique(id)){
        idi <- id[id==i]       
        lpi <- lp[id==i]
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



GenDatOMTM1.ppo <- function(id, 
                        XMat, #### design matrix that should NNOT include any intercepts
                        alpha, 
                        beta, 
                        gamma.mat,
                        ppo.mat.indx,
                        beta.ppo){
    
    # id = id
    # XMat=XMat
    # alpha=alpha
    # beta=beta
    # gamma.mat=gamma.mat
    # ppo.mat.indx = matrix(c(1,2,1,3), byrow=TRUE, ncol=2)
    # beta.ppo = c(0.15, 0.15)
    # 
    x = XMat
    n.ppo = nrow(ppo.mat.indx)
    ppo.k = ppo.mat.indx[,1]
    ppo.x = ppo.mat.indx[,2]
    
    lp         <- x %*% beta
    gamma.mat0 <- cbind(gamma.mat, 0)
    K1         <- length(alpha)
    K          <- K1+1
    start      <- 1
    
    yval <- rep(NA, length(lp))
    for (i in unique(id)){
        idi <- id[id==i]       
        lpi <- lp[id==i]
        mi  <- length(idi)
        xi = x[id==i,]
        
        ## marginal mean / probability matrix and lagged values
        ## rows correspond to time
        ## columns correspond to outcome category
        lpi.mat <- matrix(NA, mi, K1)
        for (j in 1:mi){ 
            xit = xi[j,]
            for( k in 1:K1){ tmp.1        = matrix(0, 1, n.ppo)
                             blah         = which(ppo.k==k)
                             tmp.1[blah]  = xit[ppo.x[blah]]
                             lpi.mat[j,k] = lpi[j]+alpha[k] + tmp.1 %*% beta.ppo }}
        
        cprobi.mat <- cbind(0, expit(lpi.mat), 1) ## cumulative prob matrix
        probi.mat  <- matrix(NA, mi, K)           ## multinomial probabilities
        for (k in 2:(K+1)){
            probi.mat[,(k-1)] <- cprobi.mat[,k]-cprobi.mat[,(k-1)]
        }
        
        ### This value has any impact based on the way I am generating the data.
        probi.matlag <- rbind(probi.mat[1,], probi.mat[-mi,])
        
        ## Calculate Deltait across all timepoints.  
        Deltait <- NULL
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
        yival                          <- rep(1e6, nrow(yi))
        for(j in 1:nrow(yi)){ yival[j] <- which(yi[j,]==1) }
        yval[c(start:(start+mi-1))]    <- yival
        start                          <- start+mi
    }
    yval
}

Calc.TransProbs <- function(Y, id, digs){
    
    L <- length(Y)
    Y1 <- c(NA,Y[-L])
    dup <- duplicated(id)
    Y <- Y[dup]
    Ylag <- Y1[dup]
    L <- length(Y)
    
    tab <- table(Y,Ylag)
    joint <- tab/L
    trans <- sweep(tab, 2, colSums(tab), "/")
 
    out <- list(Tx.Mtx.Col2Row=round(trans,digs), Joint.Mtx=round(joint, digs), Marg.Prob=round(table(Y)/L, digs))
    out
}

CreateTranIndMtx <- function(Y, id){
    
    K <- length(unique(Y))
    L <- length(Y)
    Y1 <- c(NA,Y[-L])
    dup <- duplicated(id)
    Y <- Y[dup]
    Ylag <- Y1[dup]
    L <- length(Y)
    
    AnyTrans <- 1*(table(Y,Ylag)>0)[-K,]
    AnyTrans
}

####################################################################################



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



logLikeCalc4 <- function(params, 
                         yval, 
                         x, 
                         ppo.mat.indx=NULL, # this is a matrix of 1s and 0s.  Right now it only permits non-ppo for all xs in the kth state
                         id, 
                         ProfileCol=NA, 
                         ref.muc=NA,
                         TransIndMtx=NA, ## used in case there are 0 probability transition.  It should be a K-1 by K matrix with 1s in all cells except for those that correspond to 0 probability transtions.
                         UseGrad=TRUE,
                         pick3=FALSE){

   # params  <- c(alpha, beta, c(.1,.1,.1), c(t(gamma.mat)))
   # params <- params1 #mod1$estimate
   # yval=Y
   # x=XMat
   # id=id
   # ProfileCol=NA
   # ref.muc=NA
   # TransIndMtx=NA ## used in case there are 0 probability transition.  It should be a K-1 by K matrix with 1s in all cells except for those that correspond to 0 probability transtions.
   # UseGrad=TRUE

    
    ## Just want to check likelihood contributions after Y enters the absorbing state
    uid1 <- unique(id[yval==1])[1:3]
    
    States    = sort(unique(yval))
    K         = length(States)
    K1        = K-1
    uid       = unique(id)
    N         = length(uid)
    Ntot      = length(yval)
    npar      = length(params)
    alpha.e   = params[1:K1]
    beta.e    = params[(K1+1):(npar-(K1^2))]
    gamma.e   = matrix(params[(npar+1-(K1^2)):npar], K1, K1, byrow=TRUE)
    gamma.mtx = cbind(gamma.e, 0)
    alpha.len = K1
    beta.len  = length(beta.e)
    gamma.len = K1^2
    
    ## matrices and indices for partial PO model
    if (!is.null(ppo.mat.indx)){
        x.ppo = x[,c(ppo.mat.indx[,2])]
        x.ppo0 = x.ppo*0
        n.ppo = ncol(x.ppo)
        ppo.k = ppo.mat.indx[,1]
        ppo.x = ppo.mat.indx[,2]
        tmp3a = matrix(0, K1, n.ppo)
    }
    
    if (!is.matrix(TransIndMtx)){ TransIndMtx = matrix(1, K1, K)}
    
    ## Z and Y matrices
    ZMat = YMat = matrix(0, Ntot, K)
    for (k in 1:K){ YMat[,k] = ifelse(yval==k, 1, 0)
                    ZMat[,k] = ifelse(yval<=k, 1, 0)
    }
    
    ## Marginal cumulative probs and state probs
    ## The else part regards partial PO
    cprob.m = prob.m = matrix(NA, Ntot, K)
    for (k in 1:K1){
        if (is.null(ppo.mat.indx)){cprob.m[,k]   = expit(alpha.e[k] + x %*% beta.e)
        }else{                     x.ppo.k       = x.ppo0
                                   aaa           = which(ppo.k==k)
                                   x.ppo.k[,aaa] = x.ppo[,aaa]
                                   cprob.m[,k]   = expit(alpha.e[k] + cbind(x,x.ppo.k) %*% beta.e)
        }
    }
    cprob.m[,K] = 1
    prob.m[,1]  = cprob.m[,1]
    for (k in 2:K){
        prob.m[,k] = cprob.m[,(k)]- cprob.m[,(k-1)]
    }
    one.min.cprob.m = 1-cprob.m
    one.min.prob.m  = 1-prob.m
    
    ## Phi and g(Phi) (we only use the first observation for each subject)
    phi.k = matrix(NA, Ntot, K1)
    for (k in 1:K1){
        phi.k[,k] = log(cprob.m[,k]/(prob.m[,(k+1)]))
    }
    g.phi.k     = log(1+exp(phi.k))
    exp.g.phi.k = exp(g.phi.k)
    
    ## dPhi/deta (we only use the first observation for each subject)
    ## dpi/deta
    dphi.k.deta.k = dphi.k.deta.k1 = dpi.k.deta.k = matrix(NA, Ntot, K1)
    for (k in 1:K1){ 
        dphi.k.deta.k[,k]  =  exp.g.phi.k[,k]*one.min.cprob.m[,k]  
        dphi.k.deta.k1[,k] =  -exp.g.phi.k[,k]*one.min.cprob.m[,(k+1)]
        dpi.k.deta.k[,k]   =  cprob.m[,k]*one.min.cprob.m[,k]
    }
    
    ## Some useful matrices to be used in gradient calculations       
    tmp1  = diag(1,K1)
    tmp2  = rbind(tmp1[-1,], rep(0,K1))
    tmp4  = matrix(0, K1, K1^2)
    
    ## To be used in the gradient calculation
    tmp.mat = matrix(0, K1, K1^2)
    for (ppp in 1:K1){ 
        tmp.mat[ppp,(c(1:K1)+(ppp-1)*K1)] = 1 
    }
    
    li1 = li2 = rep(0,N)
    Grad1Vec = Grad2Vec = rep(0,npar)
    Grad1Mat = Grad2Mat = matrix(0,N, npar)
    
    for (i in 1:N){
       
        yival            = yval[id==uid[i]]
        xi               = x[id==uid[i],]
        mi               = length(yival) 
        ZiMat            = ZMat[id==uid[i],]
        YiMat            = YMat[id==uid[i],]
        cprobi.m         = cprob.m[id==uid[i],]
        probi.m          = prob.m[id==uid[i],]
        one.min.cprobi.m = one.min.cprob.m[id==uid[i],]
        phii.k           = phi.k[id==uid[i],]
        g.phii.k         = g.phi.k[id==uid[i],]
        dphii.k.deta.k   = dphi.k.deta.k[id==uid[i],]
        dphii.k.deta.k1  = dphi.k.deta.k1[id==uid[i],]
        dpii.k.deta.k    = dpi.k.deta.k[id==uid[i],]
        
        
        if (mi==1){ xi    = matrix(xi,nrow=1)
                    ZiMat = matrix(ZiMat,nrow=1)
                    YiMat = matrix(YiMat,nrow=1)
                    cprobi.m = matrix(cprobi.m, nrow=1)
                    probi.m  = matrix(probi.m, nrow=1)
                    phii.k           = matrix(phii.k, nrow=1)
                    g.phii.k         = matrix(g.phii.k, nrow=1)
                    dphii.k.deta.k   = matrix(dphii.k.deta.k, nrow=1)
                    dphii.k.deta.k1  = matrix(dphii.k.deta.k1, nrow=1)
                    dpii.k.deta.k    = matrix(dpii.k.deta.k, nrow=1)
        }
        
        logLi1    = logLi2 = 0
        logLij2   = rep(0, mi)
        Grad1iVec = Grad2iVec = rep(0, npar)
        
        
        ## LogLi1 and dLogLi1/d(alpha,beta)
        ## LogLi1 only comes from the marginal portion of the model (first observation)
        for (k in 1:K1){ #print(cprobi.m[1,]);print(probi.m[1,])
            logLi1 = logLi1 +  ZiMat[1,k]*phii.k[1,k] - ZiMat[1,(k+1)]*g.phii.k[1,k]
        }
        
        ## LogLi1/dtheta 
        ## the else term involves the partial proportional odds part of the model
        xi1   = xi[1,]
        tmp3  = matrix(rep(xi1,K1), nrow=K1, byrow=TRUE)
        #tmp3 = matrix(c(1:9),3,3, byrow=TRUE)
      
        if (is.null(ppo.mat.indx)){ deta.k.dtheta  = cbind(tmp1,tmp3,tmp4)
                                    deta.k1.dtheta = cbind(tmp2,tmp3,tmp4)
        }else{ 
            tmp3b = tmp3a
            for (rrr in 1:n.ppo){ tmp3b[ppo.k[rrr],rrr] = xi1[ppo.x[rrr]]}
            deta.k.dtheta  = cbind(tmp1,tmp3,tmp3b, tmp4)
            deta.k1.dtheta = cbind(tmp2,tmp3,tmp3b, tmp4)
        }
    
        for (k in 1:K1){ 
            Grad1iVec = Grad1iVec+(ZiMat[1,k]-(ZiMat[1,(k+1)]*cprobi.m[1,k]/cprobi.m[1,(k+1)]))*(dphii.k.deta.k[1,k]*deta.k.dtheta[k,] + dphii.k.deta.k1[1,k]*deta.k1.dtheta[k,])
        }   
        
        if (mi>1){
            ## logLi2 comes from observations 2 to mi
            YiMat2 = YiMat[,-K]
        
            ############ Altered for new ref.muc (reference state in conditional model)
            if (!is.na(ref.muc)){
                probi.m   = cbind(probi.m[,-ref.muc], probi.m[,ref.muc])
                YiMat.tmp = cbind(YiMat[,-ref.muc], YiMat[,ref.muc])
                YiMat2    = YiMat.tmp[,-K] 
                yival2    = ifelse(yival<ref.muc, yival,
                             ifelse(yival==ref.muc, K, yival-1))
                yival     = yival2
            }
        
            for (j in 2:mi){
                mm           = probi.m[j,]
                mm.lag       = probi.m[(j-1),]
                Yit          = YiMat2[j,]
                Yit1         = YiMat2[(j-1),]
                xit          = xi[j,]
                xit1         = xi[(j-1),]
                dpi.deta     = dpii.k.deta.k[j,]
                dpi.deta.lag = dpii.k.deta.k[(j-1),]
                
                #############################################################################################################
                #Deltaij2    = Delta_calc(mm=mm, mm.lag=mm.lag, gamma.mat=gamma.mtx)
                Deltaij2    = Delta_calcWith0TxProbs(mm=mm, mm.lag=mm.lag, gamma.mat=gamma.mtx, CalcTxMtx=TransIndMtx)
                #############################################################################################################
                
                #lpij.c     = Deltaij2 + t(gamma.e) %*% YiMat2[(j-1),]
                lpij.c      = Deltaij2 + gamma.mtx[,yival[(j-1)]] #%*% YiMat2[(j-1),]
                probij.cK   = 1/(1+sum(exp(lpij.c)))
                
                logLij2[j]  = Yit %*% lpij.c + log(probij.cK)
                muit.c      = exp(lpij.c)*probij.cK
            
                #############################################################################################################
                #hmat = hmat_calc(Deltaij2, gamma.mtx)
                hmat = hmat_calc(Deltaij2, gamma.mtx)*TransIndMtx
                
                #############################################################################################################
                
                # df.dDelta   = matrix(0, K1, K1)
                # ## Left side of system of equations that solve for dDelta/dtheta
                # for (l in 1:K1){df.dDelta[l,l] = sum(hmat[l,]*(1-hmat[l,])*mm.lag)}
                # ## upper triangle for df.dDelta
                # for (l in 1:(K1-1)){
                #     for (m in (l+1):K1){
                #         df.dDelta[l,m] = -sum(hmat[l,]*hmat[m,]*mm.lag)
                #     }}
                # ## lower triangle for
                # df.dDelta[lower.tri(df.dDelta)] =  t(df.dDelta)[lower.tri(df.dDelta)]
                df.dDelta = dfdDelta_calc(mm.lag=mm.lag, hmat=hmat)
                
                ## right side of system of equations that solve for dDelta/dtheta (only alpha and beta)
                ## the else part is for partial PO
                if (is.null(ppo.mat.indx)){
                    dpidtheta.lag      = dpidtheta_calc1(tmp1, dpi.deta.lag, xit1)
                }else{
                    tmp3b = tmp3a
                    for (rrr in 1:n.ppo){ tmp3b[ppo.k[rrr],rrr] = xit1[ppo.x[rrr]]
                    }
                    dpidtheta.lag      = dpidtheta_calc2(tmp1, dpi.deta.lag, xit1,tmp3b)
                }
                ########################################################################
                dmumdtheta.lag = dmumdtheta_calc(dpidtheta.lag)
                ## for new ref.muc #############################
                if (!is.na(ref.muc)) dmumdtheta.lag = rbind(dmumdtheta.lag[-ref.muc,], dmumdtheta.lag[ref.muc,])
                h.dmumdtheta.lag     = hmat %*% dmumdtheta.lag

                ######################################################################## 
                ## right side of system of equations that solve for dDelta/dtheta (only alpha and beta)
                ## the else part is for partial PO
                if (is.null(ppo.mat.indx)){
                    dpidtheta         = dpidtheta_calc1(tmp1, dpi.deta, xit)
                }else{
                    tmp3b = tmp3a
                    for (rrr in 1:n.ppo){ tmp3b[ppo.k[rrr],rrr] = xit[ppo.x[rrr]]}
                    dpidtheta         = dpidtheta_calc2(tmp1, dpi.deta, xit,tmp3b)
                }
                
                dmumdtheta = dmumdtheta_calc(dpidtheta)
                ## for new ref.muc #############################
                if (is.na(ref.muc)){dmumdtheta  = dmumdtheta[-K,]
                }else{              dmumdtheta  = dmumdtheta[-ref.muc,]}
            
                ########################################################################
                rhs.alpha.beta  = dmumdtheta - h.dmumdtheta.lag
                dhdgamma.mm.lag = dhdgamma_mmlag_calc(hmat, mm.lag, tmp.mat)
                rhs             = cbind(rhs.alpha.beta, dhdgamma.mm.lag)
                ddelta.dtheta   = solve(df.dDelta, rhs)
            
                GradiAdd = (Yit-muit.c) %*% ddelta.dtheta + c(rep(0, alpha.len+beta.len), c(rep(Yit-muit.c, each=K1) * rep(Yit1, K1)))
                Grad2iVec = Grad2iVec + GradiAdd
                }
        }
        # if (min(yival)==1){
        #     print(uid[i])
        #     print(yival)
        #     print(logLij2)
        #     break("stop")
        # }
        
        li1[i]       = -1*logLi1
        li2[i]       = -1*sum(logLij2)
        Grad1Mat[i,] = -1*Grad1iVec
        Grad2Mat[i,] = -1*Grad2iVec
        
        if (pick3==TRUE){
            if (i %in% uid1){
                print(round(c(i,logLi1),4))
                print(round(rbind(yival, logLij2),4))
            }}

    }

    Grad1Vec = apply(Grad1Mat,2,sum)
    Grad2Vec = apply(Grad2Mat,2,sum)
    GradVec  = Grad1Vec+Grad2Vec
    GradVec[ProfileCol] = 0

    li       = sum(li1)+sum(li2)
    out      = sum(li)
    if (UseGrad==TRUE) attr(out,"gradient") = GradVec
    out
}

CalcVarCov <- function(MOD, epsilon, yval, x, id, ProfileCol=NA, ref.muc=NA, TransIndMtx=NA, UseGrad=TRUE,  ppo.mat.indx=NULL){
    npar        <- length(MOD$estimate)
    eps.mtx     <- diag(rep(epsilon, npar))
    grad.at.max <- MOD$gradient
    ObsInfo.tmp <- ObsInfo <- matrix(NA, npar, npar)
    
    ## Observed Information
    for (j in 1:npar){
        temp            <- logLikeCalc4(MOD$estimate+eps.mtx[j,], yval=yval, x=x, id=id, ProfileCol = ProfileCol, ref.muc=ref.muc, TransIndMtx=TransIndMtx, UseGrad=UseGrad,  ppo.mat.indx=ppo.mat.indx)
        ObsInfo.tmp[j,] <- (attr(temp,"gradient")-grad.at.max)/(epsilon)
    }
    for (m in 1:npar){ for (n in 1:npar){ ObsInfo[m,n] <-  (ObsInfo.tmp[m,n]+ObsInfo.tmp[n,m])/2}}
    VarCov <- solve(ObsInfo)
    VarCov
}


cluster.summary = function( id, x, fun ){
    xlist = split( x, id )
    nj = unlist( lapply( xlist, length ) )
    xj = unlist( lapply( xlist, fun) )
    xsummary = rep( xj, nj )
    xsummary}

