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

The `r sparse.rmvn()` and `r sparse.dmvn()` do most of the work for ramdon sampling and computing densities.  But the user does have the opportunity/responsibility for doing a little bit of work in advance.  Instead of supplying a dense, symmetric covariance matrix, the user must precompute a Cholesky decomposition of a matrix, which can be either a covariance or precision.  For example, suppose $$x\sim MVN(\mu,\Sigma)$$, where $$x\in\mathbb{R}^k$$, $$\mu\in\mathbb{R}^k$$, $$\Sigma\in\mathbb{R}^{k\times k}$$ and $$k$$ is very large.

```r
## mu is a numeric vector of length k
## Sigma is a k x k covariance matrix
CH <- Matrix::Cholesky(Sigma)
x <- sparse.rmvn(10, mu, CH, prec=FALSE) ## 10 random draws of x
d <- sparse.dmvn(A, mu, CH, prec=FALSE) ## densities of the 10 draws
```
