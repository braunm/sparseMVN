# sparseMVN

The [sparseMVN](braunm.github.io/sparseMVN/) package provides standard multivariate normal (MVN) sampling and density algorithms that are optimized for sparse covariance and precision matrices.

# Installation

The package can be installed from CRAN:

```
install.packages("sparseMVN")
```
or from GitHub

```
install_github("braunm/sparseMVN")
```

# Using the package

To establish some definitions, let $x\sim MVN(\mu,\Sigma)$ be a k-dimensional MVN random variate, where $\mu$ is a mean vector, $\Sigma$ is a $k \times k$ covariance matrix, and its inverse, $\Sigma^{-1}$ is the precision matrix.

When the distinction is unimportant, $\Sigma^*$ will represent  either $\Sigma$ or $\Sigma^{-1}$.

The `rmvn.sparse()` function samples $x$ from this distribution, and  `dmvn.sparse()` computes its density.  They work much like `mvtnorm::rmvnorm` and `mvtnorm::dmvnorm`,  except that instead of providing the covariance matrix $\Sigma$ as a garden-variety, base R matrix, the user provides a sparse Cholesky factorization of either the covariance or the precision.

The reasons why this is a good strategy are detailed in the vignette, but put simply, sparseMVN does the following:

1.  It avoids storing redundant or duplicate data in memory.  Even in a dense symmetric matrix, all of the $\binom{k-1}{2}$ in the upper triangle are duplicated in the lower triangle.  Also, elements in a large sparse matrix are nearly all zeros, so why store them?  It is more efficient to compress only the nonzero values into a format that points to where in the matrix those values belong.  Then, we can use sparse linear algebra routines on those matrices (dense routines would perform a lot of multiply-by-zero operations, which are clearly wasteful).

2.  It lets the user factor a matrix just once, and reuse that factorization in subsequent calls to `rmvn.sparse()` and `dmvn.sparse()`.   Under the hood, every call to `rmvnorm` or `dmvnorm` involves factoring a matrix, which is repetitive if that matrix does not change.


Here's an example. Suppose mu and S are the mean and covariance of an MVN random variable, and that S is symmetric (obviously) and sparse, but stored as a typical base R matrix. First, the user has to coerce S into a `dsCMatrix` object, which is symmetric, sparse, and column-compressed. Then, the user creates a `CHMsuper` or `CHMsimpl` object containing information about the sparse Cholesky decomposition.

```
CH <- Cholesky(as(S, 'dsCMatrix'))
```
Note that you cannot use `chol` instead.  CH is not a matrix, but a structure containing information about a matrix factorization.



After that, getting MVN random samples, and computing the densities of MVN variates, is straightforward.

```
x <- rmvn.sparse(10, mu, CH, prec=FALSE) ## 10 random draws of x
d <- dmvn.sparse(x, mu, CH, prec=FALSE) ## densities of the 10 draws
```
Note that CH was computed only once, but used in two different function calls.


If S were a precision matrix instead, the Cholesky step is unchanged, but the `prec` argument in `rmvn.sparse` and `dmvn.sparse` would be TRUE.  By default, `dmvn.sparse` returns a log density, but that can be overridden with a `log=FALSE` argument.
