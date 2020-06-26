

sourceCpp("../src/Delta_calcWith0TxProbs.cpp")

Delta_calcWith0TxProbs <- function(mm, mm.lag, gamma.mat, CalcTxMtx, tol=1e-4, maxit=10000, trace=0)
  Delta_calcWith0TxProbs_cpp(mm, mm.lag, gamma.mat, CalcTxMtx, tol, maxit, trace) 
