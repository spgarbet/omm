#include <Rcpp.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>

using namespace Rcpp;

// In place hmat calculation, not exported to R
void hmat_calc(
  NumericMatrix hmat,      // (k-1, k)  OUT
  NumericVector Deltavec,  // k-1       IN
  NumericMatrix gammamat   // (k-1, k)  IN
)
{
  int    km1=Deltavec.length();
  int    k=km1+1; 
  int    row,col;
  double denom;

  // Delta.mat   <- matrix(rep(Delta.vec,each=K), ncol=K, byrow=TRUE)
  // hmat.num    <- exp(Delta.mat+gamma.mat)
  // hmat.denom  <- 1 + colSums(hmat.num)
  // hmat        <- sweep(hmat.num,2,hmat.denom,"/")
  for(col=0; col<k; ++col)  // outer-loop due to column-wise nature of calc
  { 
    denom = 1.0; // Denominator for column
    for(row=0; row<km1; ++row) 
    {
      // matrix memory is a column-wise in layout
      // vector sequential in memory.
      // i.e., all matrices appear as as.vector(matrix) in memory
      // location = row + col*nrow 
      hmat[row+col*km1] = exp(Deltavec[row]+gammamat[row+col*km1]);
      denom += hmat[row+col*km1];
    }
    for(row=0; row<km1; ++row) hmat[row+col*km1] /= denom;
  }
}


// Exported to R, allocates memory for hmat
// [[Rcpp::export()]]
NumericMatrix hmat_calc_rcpp(// (k-1, k)
  NumericVector Deltavec,    //  k-1
  NumericMatrix gammamat     // (k-1, k)
)
{
  NumericMatrix hmat(gammamat.nrow(), gammamat.ncol());
  
  hmat_calc(hmat, Deltavec, gammamat);
  
  return(hmat);
}

void dfd_delta_calc(
  NumericMatrix df_dDelta, // (k-1, k-1)  OUT
  NumericMatrix hmat,      // (k-1, k)    IN
  NumericVector fDelta,    // k-1         IN
  NumericVector mm,        // k           IN
  NumericVector mm_lag,    // k           IN
  int           trace      //             IN
)
{
  int           k   = hmat.ncol();
  int           km1 = hmat.nrow();   
  int           i;
  int           row;
  int           col;
  double        temp;
  IntegerVector ipiv(km1);      // Matrix pivot for LU decomposition
  NumericMatrix work(km1, km1); // Working memory needed by LAPACK

  // df.dDelta   <- matrix(0, K1, K1)
  // memset is a very fast shortcut that needs the number
  // of bytes to overwrite with a byte (0). 
  // IEEE754 all zero bytes is a zero for a double.
  // sizeof returns the number of bytes of a double,
  memset(&df_dDelta[0], 0, km1*km1*sizeof(double));
  
  // for (l in 1:K1)
  //   df.dDelta[l,l] <- sum(hmat[l,]*(1-hmat[l,])*mm.lag)
  for(row=0; row<km1; ++row)
  {
    for(col=0; col<k; ++col)
    {
      // Use packed "U" col major storage, 
      // Bytes are stored using triangle number offsets
      // https://software.intel.com/content/www/us/en/develop/documentation/mkl-developer-reference-c/top/lapack-routines/matrix-storage-schemes-for-lapack-routines.html
      df_dDelta[row+row*(row+1)/2] += hmat[row+col*km1] * (1-hmat[row+col*km1])*mm_lag[col];
    }
  }

  // ## upper triangle for df.dDelta
  // for (l in 1:(K1-1))
  // {
  //   for (m in (l+1):K1)
  //   {
  //     df.dDelta[l,m] <- -sum(hmat[l,]*hmat[m,]*mm.lag)
  //   }
  // }
  for(row=0; row<(km1-1); ++row)   // for (l in 1:(K1-1))
  {
    for(col=row+1; col<km1; ++col)   //   for (m in (l+1):K1)
    {
      for(i=0; i<k; ++i)  // Elements to sum
      {
        df_dDelta[row+col*(col+1)/2] -= hmat[row+i*km1]*hmat[col+i*km1]*mm_lag[i];
      }
    }
  }
  
  //del <- solve(df.dDelta) %*% fDelta
  // Solve inverse of real symmetric indefinite matrix in packed storage
  // https://www.math.utah.edu/software/lapack/lapack-d/dsptrf.html
  dsptrf_("U", &km1, &df_dDelta[0], &ipiv[0], &i);
  if(i != 0) Rcpp::stop("Bunch-Kaufman factorization failed: info %d", i);
  // https://www.math.utah.edu/software/lapack/lapack-d/dsptri.html
  dsptri_("U", &km1, &df_dDelta[0], &ipiv[0], &work[0], &i);
  if(i != 0) Rcpp::stop("Inversion failed, DSPTRI Info %d", i);
  
  if(trace > 1)
  {
    Rprintf("\n  dfd.Delta:\n");
    for(col=0; col<km1; ++col)
    {
      Rprintf("    ");
      for(row=0; row<km1; ++row)
      {
        Rprintf("%13.6le ", df_dDelta[row + col*km1]); 
      }
      Rprintf("\n");
    }
  }
}

// [[Rcpp::export()]]
NumericVector delta_it_cpp(
  NumericVector mm,       // k
  NumericVector mm_lag,   // k
  NumericMatrix gammamat, // (k-1, k)
  double        tol,
  int           maxit,
  int           trace
)
{
  int    k=mm.length();
  int    km1=k-1;      // K - 1
  int    i;           
  int    row, col;
  int    itr;          // Track iterations
  double one=1.0;
  double zero=0.0;
  double maxslope;     // Finds maximum slope
  
  NumericMatrix hmat(km1, k);
  NumericMatrix df_dDelta(km1, km1);
  NumericVector del(km1);
  NumericVector fDelta(km1);
  NumericVector Deltavec(km1);
  
  // Double check assumptions about size
  if(mm.length() != mm_lag.length())
    Rcpp::stop("mm size != mm.lag size\n");
  if(gammamat.nrow() != km1)
    Rcpp::stop("gamma rows incorrect size compared with mm\n");
  if(gammamat.ncol() != k)
    Rcpp::stop("gamma columns incorrect size compared with mm\n");
  
  itr=0;
  maxslope=tol+0.1; // Make sure first iteration happens
                
  while(maxslope > tol)
  {
    if(++itr > maxit) Rcpp::stop("Maximum Iterations Exceeded");
    
    if(trace > 0) Rprintf("Iteration %d, ", itr);

    // hmat [OUT], Deltavec[IN], gammamat[IN]
    hmat_calc(hmat, Deltavec, gammamat);
        
    // fDelta <- hmat %*% mm.lag - mm
    for(i=0; i<km1; ++i) fDelta[i] = -mm[i]; // BLAS will modify this vector
    i    = 1;   // Each entry in vector is one apart, BLAS allows arbitrary spacing
    dgemv_("N", &km1, &k, &one, &hmat[0], &km1, &mm_lag[0], &i, &one, &fDelta[0], &i);
    // Ref: http://www.netlib.org/lapack/explore-html/d7/d15/group__double__blas__level2_gadd421a107a488d524859b4a64c1901a9.html#gadd421a107a488d524859b4a64c1901a9
  
    // df_dDelta [OUT], hmat [IN], fDelta [IN], mm [IN], mm_lag[IN], trace[IN]
    dfd_delta_calc(df_dDelta, hmat, fDelta, mm, mm_lag, trace);
    
    // fDelta <- hmat %*% mm.lag - mm
    i=1;       // Elements are next to each other (i.e. no comb like gaps)
    // This is the  "%*% fDelta" piece with a triangular packed matrix
    // Result is left in del
    //http://www.netlib.org/lapack/explore-html/d7/d15/group__double__blas__level2_gab746575c4f7dd4eec72e8110d42cefe9.html#gab746575c4f7dd4eec72e8110d42cefe9
    dspmv_("U", &km1, &one, &df_dDelta[0], &fDelta[0], &i, &zero, &del[0], &i);
    
    // Modify Deltavec elementwise and find maximum absolute slope.
    maxslope = 0.0;
    if(trace > 0) Rprintf("  del: ");

    for(i=0; i<km1; ++i)
    {
      // Delta.vec  <- Delta.vec - del
      Deltavec[i] -= del[i];

      if(fabs(del[i]) > maxslope) maxslope = fabs(del[i]); 
      
      if(trace > 0) Rprintf("%le ", del[i]);
    }
    if(trace > 0) Rprintf("\n");

  }
  
  Deltavec.attr("hmat")      = hmat;
  
  // Unpack df_dDelta from Triangular Upper Packed to Full Symmetric (why no dtpttr?)
  for(col=km1-1; col >= 0; --col)
  {
    for(row=km1-1; row >= 0; --row)
    {
      df_dDelta[row+col*km1] = df_dDelta[row+col*(col+1)/2];
    }
  }
  for(col=0; col < (km1-1); ++col)
  {
    for(row=col+1; row < km1; ++row)
    {
      df_dDelta[row+col*km1] = df_dDelta[col+row*km1];
    }
  }

  Deltavec.attr("df.dDelta") = df_dDelta;
  return(Deltavec);
}


// [[Rcpp::export]]
NumericVector omm_loglikelihood(
    NumericVector params,
    NumericVector yval,
    NumericMatrix x,
    NumericVector id, // Must be clustered and grouped by id
    NumericMatrix z_mat,
    NumericMatrix y_mat,
    IntegerVector states,
    IntegerVector uid,
    bool          grad)
{
  int           k    = states.length();
  int           km1  = k-1;
  int           n    = uid.length();
  int           npar = params.length();
  int           i;
  int           offset = 0; // Loop through ids. 
  NumericVector gradient(npar);
  NumericVector ll(1);
  
  ll[0] = 0.0; // Initialize ll to zero
  
  //alpha.e <- params[1:K1]
  NumericVector alpha_e=params[Rcpp::Range(0, km1-1)];
  
  //beta.e  <- params[(K1+1):(npar-(K1^2))]
  i = npar - km1*km1 - 1;
  NumericVector beta_e =params[Rcpp::Range(km1, i)];
  
  //gamma.e   <- matrix(params[(npar+1-(K1^2)):npar], K1, K1, byrow=TRUE)
  NumericMatrix gamma(km1, k);
  memcpy(&gamma[0], &params[i+1], km1*km1*sizeof(double));
  //gamma.mtx <- cbind(gamma.e, 0)
  memset(&gamma[km1*km1], 0, km1*sizeof(double));
  
  
  ll[0] *= -1.0;
  if(grad) ll.attr("gradient") = gradient;
  return(ll);
}

