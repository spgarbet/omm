library(inline)
  
##############################################################################
##
## Now for the C translation
##
inner_c <- '
  int    km1;      // K - 1
  int    i,j;      // i denotes row, j denotes column, 
  int    itr;      // Track iterations
  double temp;     // Useful double for temp storage
  double maxslope; // Finds maximum slope
  
  itr=0;
  maxslope=(*tol)+0.1; // Make sure first iteration happens
  
  km1 = *k - 1; // Inline has all inputs as pointers to memory
                // *k accesses the value at that pointer
                // Unfortunately this makes "*" context dependent in C/C++
  while(maxslope > *tol)
  {
    if(++itr > *maxit) error("Maximum Iterations Exceeded");
    
     if(*trace > 0) Rprintf("Iteration %d, ", itr);
    
    // Delta.mat <- matrix(rep(Delta.vec,each=nlevels), ncol=nlevels, byrow=TRUE)
    // hmat.num <- exp(Delta.mat+gamma.mat)
    // hmat.denom <- 1+colSums(hmat.num)
    // hmat <- sweep(hmat.num,2,hmat.denom,"/")
    for(j=0; j<(*k); ++j)  // j is the column, outer-loop due to column-wise nature
    { 
      temp = 1.0; // Denominator for column
      for(i=0; i<km1; ++i) // i is the row
      {
        // matrix memory is a column-wise in layout
        // vector sequential in memory.
        // i.e., all matrices appear as as.vector(matrix) in memory
        // location = row + col*nrow 
        hmat[i+j*km1] = exp(Deltavec[i]+gammamat[i+j*km1]);
        temp += hmat[i+j*km1];
      }
      for(i=0; i<km1; ++i) hmat[i+j*km1] /= temp;
    }
 
    // fDelta <- hmat %*% mm.lag - mm
    for(i=0; i<km1; ++i) fDelta[i] = -mm[i]; // BLAS will modify this vector
    i    = 1;   // Each entry in vector is one apart, BLAS allows arbitrary spacing
    temp = 1.0; // Operations are scaled by 1.0
    dgemv_("N", &km1, k, &temp, hmat, &km1, mmlag, &i, &temp, fDelta, &i);
    // Ref: http://www.netlib.org/lapack/explore-html/d7/d15/group__double__blas__level2_gadd421a107a488d524859b4a64c1901a9.html#gadd421a107a488d524859b4a64c1901a9
  
    // df.dDelta <- rep(0, nlevelsmin1)
    // memset is a very fast shortcut that needs the number
    // of bytes to overwrite with a byte (0). 
    // IEEE754 all zero bytes is a zero for a double.
    // sizeof returns the number of bytes of a double,
    memset(dfdDelta, 0, km1*sizeof(double));
  
    // hmat <- (1-hmat)*hmat 
    for(i=0; i<km1*(*k); ++i)
    {
      temp    = hmat[i];      // Modification in place is not possible in C/C++
      hmat[i] = temp*(1-temp);// if it references itself, so using temp
    }
  
    //for(i in 1:nlevelsmin1)
    for(i=0; i<km1; ++i)
    {
      //for(j in 1:nlevels)
      for(j=0; j<(*k); ++j)
      {
        // df.dDelta[i] <- df.dDelta[i] + hmat[i,j] * mm.lag[j]
        dfdDelta[i] += hmat[i+j*km1] * mmlag[j];
      }
    }

    // These operations are all elementwise, so one loop
    maxslope = 0.0;
    if(*trace > 0) Rprintf("  del: ");
    for(i=0; i<km1; ++i)
    {
      // del <- 1/df.dDelta * fDelta
      del[i] = 1/dfdDelta[i] * fDelta[i];
  
      // Delta.vec  <- Delta.vec - del
      Deltavec[i] -= del[i];
      
      // This recreates the original
      // if(del[i] > maxslope) maxslope = del[i];
      
      // This I think is the correct version
      if(fabs(del[i]) > maxslope) maxslope = fabs(del[i]); 
      
      if(*trace > 0) Rprintf("%le ", del[i]);
    }
    if(*trace > 0) Rprintf("\\n");

  }
'

# This defines the interface between R and C
inner   <- cfunction(
  signature( k        = "integer", # Categories
             gammamat = "double",  # (k-1, k)
             mm       = "double",  # k-1 (okay to pass full k)
             mmlag    = "double",  # k
             tol      = "double",  # Convergence criteria
             maxit    = "integer", # iteration cutoff
             trace    = "integer", # debugging level, 0 or 1
    
             # These are local variables, R is allocating
             hmat     = "double",  # (k-1, k)
             dfdDelta = "double",  # k-1
             del      = "double",  # k-1
             fDelta   = "double",  # k-1
             Deltavec = "double"   # k-1
           ),
           inner_c,
           includes   = c("#include <R_ext/Lapack.h>","#include <R_ext/BLAS.h>"),
           language   = "C",
           convention = ".C")

# Let's put an R wrapper on it that does the memory
delta_it_c <- function(mm,              # row of marginal (mean) probabilities
                                        # by response category k=1...K 
                       mm.lag,          # lagged value of mm
                       gamma.mat,       # Transition log odds ratios
                       k          = length(mm),
                       tol        = 1e-4,
                       maxit      = 100,
                       trace      = 0)
{
  # Using inline all inputs are outputs, but we only want one
  inner(k, gamma.mat, mm, mm.lag, tol, maxit, trace,
        rep(0, k*k-1), rep(0, k-1), rep(0,k-1), rep(0,k-1), rep(0, k-1))$Deltavec
}

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
    while (max(abs(del))>1e-4){
        
        Delta.mat   <- matrix(rep(Delta.vec,each=nlevels), ncol=nlevels, byrow=TRUE)
        hmat.num    <- exp(Delta.mat+gamma.mat)
        hmat.denom  <- 1+ colSums(hmat.num)
        hmat        <- sweep(hmat.num,2,hmat.denom,"/")
        fDelta      <- hmat %*% mm.lag - mm
        df.dDelta   <- matrix(0, nlevelsmin1, nlevelsmin1)
        
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
    #id        <- rep(1:N, each=mi)
    K1        <- length(alpha)
    K         <- K1+1
    
    ## linear predictor without the intercepts
    #lp          <- cbind(tx,t,tx*t) %*% beta
    y <- yval <- NULL
    
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
        probi.matlag <- rbind(c(0, rep(1/(K-2),K-2), 0), probi.mat[-mi,])
        
        ## Calculate Deltait across all timepoints.  
        Deltait <- NULL
        #for (j in 1:mi){ Deltait <- rbind(Deltait, c(findDeltait(mm=probi.mat[j,], mm.lag=probi.matlag[j,], gamma.mat=gamma.mat, K=K)))}
        for (j in 1:mi){ Deltait <- rbind(Deltait, c(delta_it_c(mm=probi.mat[j,], mm.lag=probi.matlag[j,], gamma.mat=gamma.mat0)))}
        
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
        y    <- rbind(y, yi)
        yval <- c(yval, yival)
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


####################################################################################
####################################################################################
####################################################################################

logLikeCalc <- function(params, yval,x, id){
    N       <- nrow(x)
    States  <- sort(unique(yval))
    K       <- length(States)
    K1      <- K-1
    uid     <- unique(id)
    npar    <- length(params)
    alpha.e <- params[1:K1]
    beta.e  <- params[(K1+1):(npar-(K1^2))]
    gamma.e <- matrix(params[(npar+1-(K1^2)):npar], K1, K1, byrow=TRUE)
    gamma.mtx <- cbind(gamma.e, 0)
    li1 <- li2 <- rep(0,N)
    print(params)
    for (i in unique(id)){
        
        yival <- yval[id==i]
        xi    <- x[id==i,]
        mi    <- nrow(xi)
        
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
            #Deltaij2       <- findDeltait(mm=probi.m[j,], mm.lag=probi.m[(j-1),], gamma.mat=gamma.mtx, nlevels=K)
            Deltaij2       <- delta_it_c(mm=probi.m[j,], mm.lag=probi.m[(j-1),], gamma.mat=gamma.mtx)
            
            #lpij.c         <- Deltaij2 + t(gamma.e) %*% YiMat2[(j-1),]
            lpij.c         <- Deltaij2 + gamma.mtx[,yival[(j-1)]] #%*% YiMat2[(j-1),]
            
            probij.cK      <- 1/(1+sum(exp(lpij.c)))
            logLij2[j] <- YiMat2[j,] %*% lpij.c + log(probij.cK)
            #logLij2[j]     <- lpij.c[yival[j]] + log(probij.cK)
            
        }
        li1[i] <- logLi1
        li2[i] <- sum(logLij2)
    }
    li <- -1*(sum(li1)+sum(li2)) 
    li
}








