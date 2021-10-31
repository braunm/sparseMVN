
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
w3z <- kronecker(w2b, diag(d))  ## sparse, banded



Q1 <- tril(w3z)
##Q1 <- tril(kronecker(Matrix(seq(0.1,p,length=p*p),p,p),diag(d)))
Q2 <- cbind(Q1,Matrix(0,d*p,k))
Q3 <- rbind(Q2,cbind(Matrix(rnorm(k*d*p),k,d*p),Diagonal(k)))

## look at it dense and sparse
CV <- Matrix::tcrossprod(Q3)




chol.CV <- Matrix::Cholesky(CV)  ## creates a dCHMsimpl object

chol_mat <- chol(CV) ## returns dtCMatrix upper triangular, losing the specific Cholesky info

stop()

k <- seq(10, 100, by=2)
R <- map_dfr(k, function(i) {
  a <- matrix(rnorm(i^2), i, i)
  b <- crossprod(a * lower.tri(a, diag=TRUE)) + .001 * diag(i)
  tm <- microbenchmark(solve(b), chol(b))
  res <- tm %>%
    summary() %>%
    as_tibble() %>%
    mutate(k=i)
  return(res)
})

P <- R %>%
  ggplot(aes(x=k, y=mean, color=expr)) +
  geom_line() +
  scale_x_continuous(breaks=seq(0, 100, by=10)) +
  scale_y_continuous(breaks=seq(0, 1000, by=50)) +
  coord_trans(y='log')
print(P)
