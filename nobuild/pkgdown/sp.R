### Some test code on manipulating sparse matrices
library(Matrix)

d <- 4 ## number of blocks
p <- 2 ## size of each block
k <- 2 ## number of variables

mu <- seq(-3,3,length=p*d+k)

w1 <-  seq(0.1,p,length=p*p)
w2a <- matrix(w1, p, p)
w2b <- Matrix(w1, p, p)

w3d <- kronecker(diag(d), w2a) ## dense
w3s <- kronecker(diag(d), w2b) ## sparse

Q1 <- tril(w3s)
##Q1 <- tril(kronecker(Matrix(seq(0.1,p,length=p*p),p,p),diag(d)))
Q2 <- cbind(Q1,Matrix(0,d*p,k))
Q3 <- rbind(Q2,cbind(Matrix(rnorm(k*d*p),k,d*p),Diagonal(k)))

## look at it dense and sparse
CV <- Matrix::tcrossprod(Q3)




chol.CV <- Matrix::Cholesky(CV)  ## creates a dCHMsimpl object

chol_mat <- chol(CV) ## returns dtCMatrix upper triangular, losing the specific Cholesky info
