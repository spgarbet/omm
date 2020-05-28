#include <Rcpp.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>

using namespace Rcpp;


// [[Rcpp::export]]
NumericMatrix LUfactor(NumericMatrix A)
{
  int r=A.nrow();
  int c=A.ncol();
  int lda=1;
  IntegerMatrix ipiv(r, c);
  int info;
  
  dgetrf_(&r, &c, &A[0], &lda, &ipiv[0], &info);
  
  return A;
}
