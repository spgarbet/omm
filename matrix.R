library(inline)
x <- as.numeric(1:10)
n <- as.integer(10)
## Not run:
## A simple Fortran example - n and x: assumed-size
code <- "
      integer i
      do 1 i=1, n(1)
    1 x(i) = x(i)**3
"
cubefn <- cfunction(signature(n="integer", x="numeric"), code, convention=".Fortran")
print(cubefn)

cubefn(n, x)$x


code4 <- "
  int i;
  for(i=0; i<*n; ++i)
    x[i] = pow(x[i],3);
"
cubefn4 <- cfunction(signature(n="integer", x="numeric"), code4, language = "C", convention = ".C")
cubefn4(20, 1:20)

code5 <- "
  x[*n,*n] = 2;
"
cubefn5 <- cfunction(signature(n="integer", x="matrix"), code5, language = "C", convention = ".C")
cubefn5(2, matrix(1:9, nrow=3))

library(inline)
# http://www.netlib.org/lapack/explore-html/db/dc9/group__single__blas__level3_gafe51bacb54592ff5de056acabd83c260.html#gafe51bacb54592ff5de056acabd83c260
lapack <- '
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
                    lapack,
                    includes = c("#include <R_ext/Lapack.h>","#include <R_ext/BLAS.h>"),
                    language = "C",
                    convention = ".C")
x <- matrix(1:6,  nrow=3, ncol=2)
y <- matrix(7:14, nrow=2, ncol=4)
z <- rep(0, 12)
ans <- mmul(x, y, z, nrow(x), ncol(y), nrow(y))$c
x %*% y
matrix(ans, nrow=3, ncol=4)
