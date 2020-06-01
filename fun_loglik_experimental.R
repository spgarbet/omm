  ##############################################################################
 #
#   Log Likelihood calculation for Ordinal Transition Matrix 
#   Experimental version with no gradient
#

logLikeCalcExp <- function(params, yval, x, id, States, K, K1, uid, N, npar)
{
    alpha.e   <- params[1:K1]
    beta.e    <- params[(K1+1):(npar-(K1^2))]
    gamma.e   <- matrix(params[(npar+1-(K1^2)):npar], K1, K1, byrow=TRUE)
    gamma.mtx <- cbind(gamma.e, 0)
    
    li1 <- li2 <- rep(0,N)

    logLi1  <- logLi2 <- 0
    logLij2 <- rep(0, mi)
    
    -1*sum(sapply(uid, FUN=function(i)
    {
      # Using split here was slower
      yival   <- yval[id==i]
      xi      <- x[id==i,]
      mi      <- nrow(xi)
      

        
      ZiMat   <- YiMat <- cprobi.m <- probi.m <- matrix(NA, mi, K)

        for (k in 1:K)
        {
          ZiMat[,k] <- as.integer(yival<= States[k])
          YiMat[,k] <- as.integer(yival== States[k])
        }
      
      # This was slower
      #ZiMat <- sapply(States, function(s) as.integer(yival <= s))
      #YiMat <- sapply(States, function(s) as.integer(yival == s))
        
        for (k in 1:K1) cprobi.m[,k] <- expit(alpha.e[k] + xi %*% beta.e)
        
        cprobi.m[,K] <- 1
        probi.m[,1]  <- cprobi.m[,1]
        for (k in 2:K)  probi.m[,k] <- cprobi.m[,(k)]- cprobi.m[,(k-1)]
        
        ## LogLi1 only comes from the marginal portion of the model (first observation)
        for (k in 1:K1)
        # sum(sapply(...  is slower!
        #logLi1 <- sum(sapply(1:K1, function(k)
        {
          #print(cprobi.m[1,]);print(probi.m[1,])
          logLi1 <- 
            logLi1 +  
            ZiMat[1,k]*log(cprobi.m[1,k] / probi.m[1,(k+1)]) - 
            ZiMat[1,(k+1)]*log(cprobi.m[1,(k+1)] / probi.m[1,(k+1)])
        }#))
        
        ## logLi2 comes from observations 2 to mi
        YiMat2 <- YiMat[,-K]
        
        for (j in 2:mi)
        {
            Deltaij2       <- delta_it(mm=probi.m[j,], mm.lag=probi.m[(j-1),], gamma.mat=gamma.mtx)
            lpij.c         <- Deltaij2 + gamma.mtx[,yival[(j-1)]] 
            probij.cK      <- 1/(1+sum(exp(lpij.c)))
            logLij2[j]     <- YiMat2[j,] %*% lpij.c + log(probij.cK)
        }
        logLi1+sum(logLij2)
    }))
    #print(sum(li1))
    #print(sum(li2))
    
    #li <- -1*(sum(li1)+sum(li2)) 
    #li
}