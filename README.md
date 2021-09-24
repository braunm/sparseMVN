# sparseMVN

The [sparseMVN](braunm.github.io/sparseMVN/) package provides standard multivariate normal (MVN) sampling and density algorithms that are optimized for sparse covariance and precision matrices.

# Installation

The package can be installed from CRAN:

```r
install.packages("sparseMVN")
```
or from GitHub

```r
install_github("braunm/sparseMVN")
```

# Using the package

To establish some definitions, let $$x\sim MVN(\mu,\Sigma)$$ be a k-dimensional MVN random variate, where $\mu$ is a mean vector, $\Sigma$ is a $k \times k$ covariance matrix, and its inverse, $\Sigma^{-1}$ is the precision matrix. When the distinction is unimportant, $\Sigma^*$ will represent  either $\Sigma$ or $\Sigma^{-1}$.




The `rmvn.sparse()` function samples $x$ from this distribution, and  `dmvn.sparse()` computes the density.  They work much like `mvtnorm::rmvnorm` and `mvtnorm::dmvnorm`,  except that:

1.  the `mvtnorm` functions require the user to provide the covariance matrix $\Sigma$ as a garden-variety, base R matrix with $k^2$ elements; while
2. the `sparseMVN` functions require  the user to provide a sparse Cholesky factorization of the matrix, but that matrix can be *either* the covariance $\Sigma$ *or* the precision $\Sigma$.

The option to work with the precision matrix is useful because sometimes, that's what the user has.  An example is a Laplace approximation of a log posterior density, where covariance of approximating MVN is the *inverse* of the Hessian.  Being able to use the Hessian directly as the precision matrix saves a matrix inversion step that can be computationally costly when $k$ is very large.


Operations on a dense and "unstructured"  representation of a covariance matrix  do not take symmetry and sparsity into account.

contains a lot of redundant information, as all  $k^2$ elements are explicitly included, but at most only $\binom{k+1}{2}$ of those elements are distinct. Further, if the matrix is sparse, most of those values are all zero anyway.  Not only do objects that store only the non-zero elements (and pointers to their locations in the matrix) require less memory, but linear algebra operations can be made much faster by taking the sparsity and symmetry into account.


do most of the work for ramdon sampling and computing densities.  But the user does have the opportunity/responsibility for doing a little bit of work in advance.  Instead of supplying a dense, symmetric covariance matrix, the user must precompute a Cholesky decomposition of a matrix,  (either a covariance or precision).  For example, suppose
```r
## mu is a numeric vector of length k
## Sigma is a k x k covariance matrix
CH <- Matrix::Cholesky(Sigma)
x <- sparse.rmvn(10, mu, CH, prec=FALSE) ## 10 random draws of x
d <- sparse.dmvn(A, mu, CH, prec=FALSE) ## densities of the 10 draws
```
