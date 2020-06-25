library(Rcpp)

K   <- 4
Km1 <- 3

# Input
tmp.mat         <- matrix((12+27):13, nrow=3)
hmat            <- matrix(1:12, nrow=3)  # (k-1, k) 
mm.lag          <- 40:43

dhdgamma_r <- function(tmp.mat, hmat, mm.lag)
{
  hmat1           <- hmat[,-K]  # Drop last column
  hmat2           <- matrix(rep(t(hmat1),Km1), nrow=Km1, byrow=TRUE)
  hmat1.K1times   <- matrix(rep(hmat1, Km1), nrow=Km1)
  hmat3           <- tmp.mat-hmat1.K1times
  mm.lag1         <- mm.lag[-K]
  mm.lag.mat      <- matrix(rep(mm.lag1, Km1^2), nrow=Km1, byrow=TRUE)
  
  -1*hmat2*mm.lag.mat*hmat3
}
dhdgamma.mm.lag <- dhdgamma_r(tmp.mat, hmat, mm.lag)

#### Take 2, rearrange operations to be "C" friendly

dhdgamma.mm.lag.2 <- matrix(rep(0, 27), nrow=3)
for(i in 1:Km1)
{
  for(j in 1:Km1)
  {
    dhdgamma.mm.lag.2[, (i-1)+(j-1)*Km1+1] <- -hmat[j, i]*mm.lag[i]
    for(k in 1:Km1)
    {
      dhdgamma.mm.lag.2[k, (i-1)+(j-1)*Km1+1] <- dhdgamma.mm.lag.2[k, (i-1)+(j-1)*Km1+1] *
        (tmp.mat[k, (i-1)+(j-1)*Km1+1] - hmat[k, i])
    }
  }
}

if(sum(abs(dhdgamma.mm.lag - dhdgamma.mm.lag.2)) > 1e-9)
  stop("Rearranging operations failed")

#### Take 3, convert to Rcpp
sourceCpp("src/dhd_gamma_calc.cpp")
dhdgamma.mm.lag.3 <- dhd_gamma_calc(hmat, mm.lag, tmp.mat)
if(sum(abs(dhdgamma.mm.lag - dhdgamma.mm.lag.3)) > 1e-9)
  stop("Rcpp version failed")


#### Time comparison
rtime <- system.time(replicate(2500000, function() dhdgamma_r(tmp.mat, hmat, mm.lag)))
ctime <- system.time(replicate(2500000, function() dhd_gamma_calc(hmat, mm.lag, tmp.mat)))
print(rtime)
print(ctime)

