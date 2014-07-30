\name{mvn.sparse}
\alias{dmvn.sparse}
\alias{rmvn.sparse}
\title{Sampling and log density functions for multivariate normal distribution, given a sparse covariance or precision matrix.}
\description{Efficient sampling and density calculation from a multivariate normal, when th covariance or precision matrix is sparse. These functions are designed for MVN samples of very large dimension.}
\usage{
dmvn.sparse(x, mu, CH, prec=TRUE)
rmvn.sparse(n, mu, CH, prec=TRUE)
}
\arguments{
	\item{x}{numeric matrix, where each row is an mvn sample.}
	\item{n}{number of draws}
	\item{mu}{mean (numeric vector)}
	\item{CH}{An object of class dCHMsimpl or dCHMsuper that represents
the Cholesky factorization of either the precision (default) or covariance matrix.  See details.}
	\item{prec}{If TRUE, CH is the Cholesky decomposition of the precision matrix.  If false, it is the decomposition for the covariance matrix.}
}
\details{
This function uses sparse matrix operations to sample from, or compute the log density of, a multivariate normal distribution  The user must compute the Cholesky decomposition first, using the Cholesky function in the Matrix package.  This function operates on a sparse symmetric matrix, and returns an object of class dCHMsimpl or dCHMsuper (this depends on the algorithm that was used for the decomposition).  This object contains information about any fill-reducing permutations that were used to preserve sparsity.  The rmvn.sparse and dmvn.sparse functions use this permutation information, even if pivoting was turned off.
}
\value{
For \code{rmvn.sparse}, a base R matrix with n rows, with each draw in each row. The number of columns is equal to the length of mu.

For \code{dmvn.sparse}, a vector with a length equal to the number of rows of \code{x}, containing the MVN log densities evaluated at each row of \code{x}.
}
\examples{
if (require(Matrix)) {
   m <- 50
   p <- 2
   k <- 6

## build sample covariance matrix
   Q1 <- tril(kronecker(Matrix(seq(0.1,p,length=p*p),p,p),diag(m)))
   Q2 <- cBind(Q1,Matrix(0,m*p,k))
   Q3 <- rBind(Q2,cBind(Matrix(rnorm(k*m*p),k,m*p),Diagonal(k)))
   V <- tcrossprod(Q3)
   CH <- Cholesky(V)

   ## assume mean at origin

   ## 10 samples from MVN
   x <- rmvn.sparse(10,rep(0,p*m+k),CH, FALSE)
   print(x)

   ## log density of first sample in x
   y <- dmvn.sparse(x[1,],rep(0,p*m+k), CH, FALSE)
   print(y)
}

}
