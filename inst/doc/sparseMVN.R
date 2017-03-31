## ----echo=FALSE, cache=FALSE-----------------------------------------------
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(mvtnorm))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(trustOptim))
suppressPackageStartupMessages(library(xtable))
knitr::render_sweave()
knitr::opts_chunk$set(prompt=TRUE, cache=FALSE,error=FALSE,
                      comment="#", collapse=FALSE) #$
options(replace.assign=TRUE, width=77, prompt="R> ",
        digits=4, xtable.include.rownames=FALSE,
        xtable.booktabs=TRUE)
sanitize <- function(x) x

## ----echo=FALSE------------------------------------------------------------
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

## ----blockarrow, echo=FALSE, prompt=FALSE----------------------------------
Mat <- as(kronecker(diag(N), matrix(1, k, k)),"sparseMatrix")
Mat <- rBind(Mat, Matrix(1, p, N*k))
Mat <- cBind(Mat, Matrix(1, k*N+p, p))
printSpMatrix(as(Mat,"nMatrix"))

## ----banded, echo=FALSE----------------------------------------------------
Mat <- kronecker(Matrix(1, k, k), diag(N))
Mat <- rBind(Mat, Matrix(1, p, N * k))
Mat <- cBind(Mat, Matrix(1, k*N+p, p))
printSpMatrix(as(Mat,"nMatrix"))

## ----echo=FALSE, results="hide"--------------------------------------------
Mat2 <- as(kronecker(diag(Q),matrix(1,k,k)),"lMatrix") %>%
    rBind(Matrix(TRUE,p,Q*k)) %>%
    cBind(Matrix(TRUE, k*Q+p, p)) %>%
    as("dgCMatrix") %>%
    as("symmetricMatrix")
A2 <- as(Mat2,"matrix")

## ----eval=FALSE, prompt=FALSE----------------------------------------------
#  rmvn.sparse(n, mu, CH, prec=TRUE)
#  dmvn.sparse(x, mu, CH, prec=TRUE, log=TRUE)

## ----results='hide'--------------------------------------------------------
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

## --------------------------------------------------------------------------
R <- 100
pm <- opt[["solution"]]
H <- -opt[["hessian"]]
CH <- Cholesky(H)
samples <- rmvn.sparse(R, pm, CH, prec=TRUE)

## --------------------------------------------------------------------------
logf <- dmvn.sparse(samples, pm, CH, prec=TRUE)

## --------------------------------------------------------------------------
Matrix::nnzero(H)
Hinv <- drop0(solve(H))
Matrix::nnzero(Hinv)

## --------------------------------------------------------------------------
logf_dense <- dmvnorm(samples, pm, as.matrix(Hinv), log=TRUE)
all.equal(logf, logf_dense)

## ----echo=FALSE------------------------------------------------------------
load("runtimes.Rdata")
RT <- group_by(runtimes,s, N, k, prec, bench.expr) %>%
    summarize(mean = mean(bench.time/10^6),
              sd = sd(bench.time/10^6)) %>%
    ungroup() %>%
    tidyr::separate(bench.expr, c("step", "storage"), sep="_")

## ----echo=FALSE,results='asis'---------------------------------------------
filter(cases, !prec) %>%
    select(`$N$`=N, `$k$`=k, `$M$`=nvars, `$M^2$`=nels,
           nnz, `nnz (LT)`=nnzLT, `\\% nnz`=pct.nnz) %>%
    xtable::xtable(digits=c(rep(0,7),3)) %>%
    print(include.rownames=FALSE,
          floating=FALSE, sanitize.colnames.function=sanitize,
          format.args=list(big.mark=","))

## ----echo=FALSE------------------------------------------------------------
tmp <- filter(cases, prec) %>%
    select(N,k,nvars,nnz)
tim1 <- filter(RT, step %in% c("rand","density")) %>%
    mutate(type=ifelse(prec,"cov","prec")) %>%
    select(-s, -prec) %>%
    tidyr::gather(stat,value, c(mean,sd)) %>%
    reshape2::dcast(type+N+k~step+storage+stat) %>%
    left_join(tmp, ., by=c("N","k"))

## ----echo=FALSE,results='asis'---------------------------------------------
filter(tim1, type=="cov") %>%
    select(-type) %>%
    xtable(digits=c(rep(0,4),rep(0,9))) %>%
    print(only.contents=TRUE,include.colnames=FALSE,
          format.args=list(big.mark=","),hline.after=NULL)

## ----echo=FALSE,results='asis'---------------------------------------------
filter(tim1, type=="prec") %>%
    select(-type) %>%
    xtable(digits=c(rep(0,4),rep(0,9))) %>%
    print(only.contents=TRUE,include.colnames=FALSE,
          format.args=list(big.mark=","),hline.after=NULL)

## ----echo=FALSE,results='asis'---------------------------------------------
tmp <- filter(cases, prec) %>%
    select(N,k,nvars,pct.nnz)
filter(RT, step %in% c("chol","solve") & prec) %>%
    select(-s, -prec) %>%
    tidyr::gather(stat,value, c(mean,sd)) %>%
    reshape2::dcast(N+k~step+storage+stat) %>%
    left_join(tmp, ., by=c("N","k")) %>%
    xtable(digits=c(rep(0,4),3,rep(1,8))) %>%
    print(only.contents=TRUE,include.colnames=FALSE,
          format.args=list(big.mark=","),hline.after=NULL)

