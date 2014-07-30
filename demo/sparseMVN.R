
N <- 100 ## number of samples
m <- 1000  ## number of blocks in sparse covariance matrix
p <- 3 ## size of each block
k <- 15  ## 
prec <- FALSE

## build block-arrow covariance/precision matrix

mu <- rep(0,p*m+k)  ## assume mean at origin
Q1 <- tril(kronecker(Matrix(seq(0.1,p,length=p*p),p,p),diag(m)))
Q2 <- cBind(Q1,Matrix(0,m*p,k))
Q3 <- rBind(Q2,cBind(Matrix(rnorm(k*m*p),k,m*p),Diagonal(k)))
CV <- tcrossprod(Q3)

tchol <- system.time(chol.CV <- Cholesky(CV))  ## creates a dCHMsimpl object

## draw random samples using rmvn.sparse
trand <- system.time(x.sp <- rmvn.sparse(N, mu, chol.CV,prec=prec))
cat("Time to sample ",N," times from rmvn.sparse: ",trand["elapsed"],"\n")


## computing log densities using dmvn.sparse
tdens <- system.time(d.sp <- dmvn.sparse(x.sp[1,], mu, chol.CV, prec=prec))
cat("Time to compute log density of first sample with dmvn.sparse: ",tdens["elapsed"],"\n")
