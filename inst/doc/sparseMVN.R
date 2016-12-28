## ----echo=FALSE----------------------------------------------------------
knitr::opts_chunk$set(collapse = FALSE, comment = "#", message=FALSE) #$
options(digits=4)
suppressMessages(library(dplyr))
suppressMessages(library(scales))
suppressMessages(library(trustOptim))
suppressMessages(library(xtable))
options(xtable.include.rownames=FALSE,
        xtable.booktabs=TRUE)

## ----echo=FALSE----------------------------------------------------------
require(sparseMVN)
require(Matrix)
require(mvtnorm)
N <- 5
k <- 2
p <- k ## dimension of mu
nv1 <- N*k+p
nels1 <- nv1^2
nnz1 <- N*k^2 + 2*p*N*k + p^2
nnz1LT <- N*k*(k+1)/2  + p*N*k + p*(p+1)/2
Q <- 1000
nv2 <- Q*k+p
nels2 <- nv2^2
nnz2 <- Q*k^2 + 2*p*Q*k + p^2
nnz2LT <- Q*k*(k+1)/2 + p*Q*k + p*(p+1)/2
options(scipen=999)

## ----blockarrow, echo=FALSE----------------------------------------------
Mat <- as(kronecker(diag(N), matrix(1, k, k)),"sparseMatrix")
Mat <- rBind(Mat, Matrix(1, p, N*k))
Mat <- cBind(Mat, Matrix(1, k*N+p, p))
printSpMatrix(as(Mat,"nMatrix"))

## ----echo=FALSE, results="hide"------------------------------------------
Mat2 <- as(kronecker(diag(Q),matrix(1,k,k)),"lMatrix") %>%
    rBind(Matrix(TRUE,p,Q*k)) %>%
    cBind(Matrix(TRUE, k*Q+p, p)) %>%
    as("dgCMatrix") %>%
    as("symmetricMatrix")
A2 <- as(Mat2,"matrix")

## ----eval=FALSE----------------------------------------------------------
#  rmvn.sparse(ndraws, mu, CH, prec=TRUE, log=TRUE)
#  dmvn.sparse(x, mu, CH, prec=TRUE, log=TRUE)

## ----results='hide'------------------------------------------------------
D <- binary.sim(N=50, k=2, T=50)
priors <- list(inv.Sigma=diag(2), inv.Omega=diag(2))
start <- rep(c(-1,1),51)
opt <- trust.optim(start,
                   fn=sparseMVN::binary.f,
                   gr=sparseMVN::binary.grad,
                   hs=sparseMVN::binary.hess,
                   data=D, priors=list(inv.Sigma=diag(2),
                                       inv.Omega=diag(2)),
                   method="Sparse",
                   control=list(function.scale.factor=-1))

## ------------------------------------------------------------------------
R <- 100
pm <- opt[["solution"]]
H <- -opt[["hessian"]]
CH <- Cholesky(H)
samples <- rmvn.sparse(R, pm, CH, prec=TRUE)

## ------------------------------------------------------------------------
logf <- dmvn.sparse(samples, pm, CH, prec=TRUE)

## ------------------------------------------------------------------------
Matrix::nnzero(H)
Hinv <- drop0(solve(H))
Matrix::nnzero(Hinv)

## ------------------------------------------------------------------------
logf_dense <- dmvnorm(samples, pm, as.matrix(Hinv), log=TRUE)
all.equal(logf, logf_dense)

## ------------------------------------------------------------------------
TM <- filter(runtimes, stringr::str_detect(step, "[rd]_")) %>%
    select(-p,-nels) %>%
    tidyr::gather(stat,value,c(mean,sd)) %>%
    reshape2::dcast(s+N+k+prec+nvars+nnz+pct.nnz~step+stat)

## ----echo=FALSE,results='asis'-------------------------------------------
filter(TM, !prec) %>%
    select(-prec) %>%
    xtable(digits=c(rep(0,6),rep(2,9))) %>%
    print(only.contents=TRUE,include.colnames=FALSE,size="small")

