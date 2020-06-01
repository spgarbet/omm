#include <Rcpp.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>

using namespace Rcpp;

// [[Rcpp::export()]]
NumericMatrix hmat_calc_cpp(
  NumericVector Deltavec,    // k-1
  NumericMatrix gammamat // (k-1, k)
)
{
  int    km1=Deltavec.length();
  int    k=km1+1;      // K - 1
  int    i,j;        // i denotes row, j denotes column,
  double temp;
  
  NumericMatrix hmat(km1, k);
    // Delta.mat   <- matrix(rep(Delta.vec,each=K), ncol=K, byrow=TRUE)
    // hmat.num    <- exp(Delta.mat+gamma.mat)
    // hmat.denom  <- 1+ colSums(hmat.num)
    // hmat        <- sweep(hmat.num,2,hmat.denom,"/")
    for(j=0; j<k; ++j)  // j is the column, outer-loop due to column-wise nature
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
 
  return(hmat);
}
