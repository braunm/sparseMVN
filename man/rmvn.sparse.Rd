% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rmvn-sparse.R
\name{rmvn.sparse}
\alias{rmvn.sparse}
\title{Sample from multivariate normal distribution}
\usage{
rmvn.sparse(n, mu, CH, prec = TRUE)
}
\arguments{
\item{n}{number of samples}

\item{mu}{mean (numeric vector)}

\item{CH}{An object of class dCHMsimpl or dCHMsuper that represents
the Cholesky factorization of either the precision (default) or covariance
matrix.  See details.}

\item{prec}{If TRUE, CH is the Cholesky decomposition of the precision
matrix.  If false, it is the decomposition for the covariance matrix.}
}
\value{
A matrix of samples from an MVN distribution (one in each row)
}
\description{
Efficient sampling and density calculation from a multivariate
normal,
when the covariance or precision matrix is sparse. These functions are
designed for MVN samples of very large dimension.
}
\section{Details}{

This function uses sparse matrix operations to sample from a multivariate normal distribution.  The user must compute
the Cholesky decomposition first, using the Cholesky function in the Matrix
package.  This function operates on a sparse symmetric matrix, and returns
an object of class dCHMsimpl or dCHMsuper (this depends on the algorithm
that was used for the decomposition).  This object contains information about
any fill-reducing permutations that were used to preserve sparsity.  The
rmvn.sparse and dmvn.sparse functions use this permutation information, even
if pivoting was turned off.
}

\examples{
   require(Matrix)
   m <- 20
   p <- 2
   k <- 4

## build sample sparse covariance matrix
   Q1 <- tril(kronecker(Matrix(seq(0.1,p,length=p*p),p,p),diag(m)))
   Q2 <- cbind(Q1,Matrix(0,m*p,k))
   Q3 <- rbind(Q2,cbind(Matrix(rnorm(k*m*p),k,m*p),Diagonal(k)))
   V <- tcrossprod(Q3)
   CH <- Cholesky(V)

   x <- rmvn.sparse(10,rep(0,p*m+k),CH, FALSE)
   y <- dmvn.sparse(x[1,],rep(0,p*m+k), CH, FALSE)

}
