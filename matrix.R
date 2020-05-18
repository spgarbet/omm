
library(inline)
# http://www.netlib.org/lapack/explore-html/db/dc9/group__single__blas__level3_gafe51bacb54592ff5de056acabd83c260.html#gafe51bacb54592ff5de056acabd83c260
blas <- '
  double alpha=1.0;
  dgemm_("N", "N", m, n, k, &alpha, a , m, b, k, &alpha, c, m);
'
mmul   <- cfunction(signature(
                      a="double",
                      b="double",
                      c="double",
                      m="integer",
                      n="integer",
                      k="integer"
                    ),
                    blas,
                    includes = c("#include <R_ext/BLAS.h>"),
                    language = "C",
                    convention = ".C")
x <- matrix(1:6,  nrow=3, ncol=2)
y <- matrix(7:14, nrow=2, ncol=4)
z <- rep(0, 12)
ans <- mmul(x, y, z, nrow(x), ncol(y), nrow(y))$c

# Did it work?
x %*% y
matrix(ans, nrow=3, ncol=4)

# Switch to LAPACK

lapack <- '
  dpptrf_("U", N, A, info);
'
ichol   <- cfunction(signature(
                      N="integer",
                      A="double",
                      info="integer"
                    ),
                    lapack,
                    includes = c("#include <R_ext/Lapack.h>"),
                    language = "C",
                    libargs="-llapack",
                    convention = ".C")

# Did it work?
ichol(2, c(2, -1, 2), 0)
chol(matrix(c(2,-1,-1,2), nrow=2))

