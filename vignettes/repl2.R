library("sparseMVN")
library("microbenchmark")
library("Matrix")
library("mvtnorm")
library("dplyr")

##library(tidyr)
library(doParallel)
set.seed(123)
cores <- 10
registerDoParallel(cores=cores)

build_mat <- function(N, k) {

    ##   Q1 <- tril(kronecker(Matrix(seq(0.1,k,length=k*k),k,k),diag(N)))
    t1 <- exp(rnorm(k*k))
    Q1 <- tril(kronecker(Matrix(t1,k,k),diag(N)))
    Q2 <- cBind(Q1,Matrix(0, N*k, k))
    Q3 <- rBind(Q2,cBind(Matrix(rnorm(N*k*k), k, N*k), Diagonal(k)))
    tcrossprod(Q3)
}

check_density <- function(CV.sparse, prec) {
    chol.CV <- Cholesky(CV.sparse)
    if (prec) sigma <- solve(CV.dense) else sigma <- CV.dense

    x.sp <- rmvn.sparse(s, mu, chol.CV, prec=prec)
    d.sp <- dmvn.sparse(x.sp, mu, chol.CV, prec=prec)
    d.dens <- dmvnorm(x.sp, mu, sigma, log=TRUE)
    all.equal(d.sp,d.dens)
}

run_bench <- function(D, reps=10) {

    s <- D$s ## number of random samples
    k <- D$k ## heterogeneous variables
    N <- D$N ## number of agents

  ##  cat("s =  ",s,"\tN = ",N,"\tk = ",k,"\n")

    mu <- rep(0,k*N + k)  ## assume mean at origin

    CV.sparse <- build_mat(N, k)
    CV.dense <- as(CV.sparse, "matrix")  ## dense covariance
    chol.CV <- Cholesky(CV.sparse)

    ## check_cov <- check_density(CV.sparse, FALSE)
    ## check_prec <- check_density(CV.sparse, TRUE)
    ## stopifnot(check_cov & check_prec)

    x <- rmvn.sparse(s, mu, chol.CV, prec=FALSE)

    bench <- microbenchmark(
        chol_sparse = Cholesky(CV.sparse),
        chol_dense = chol(CV.dense),
        solve_dense = solve(CV.dense),
        rand_sparse_cov = rmvn.sparse(s, mu, chol.CV, prec=FALSE),
        rand_sparse_prec = rmvn.sparse(s, mu, chol.CV, prec=TRUE),
        density_sparse_cov = dmvn.sparse(x, mu, chol.CV, prec=FALSE),
        density_sparse_prec = dmvn.sparse(x, mu, chol.CV, prec=TRUE),
        rand_dense = rmvnorm(s, mu, CV.dense, method="chol"),
        density_dense = dmvnorm(x, mu, CV.dense, log=TRUE),
        times = reps
    )

    vals <- plyr::ddply(data.frame(bench), "expr",
                        function(x) return(data.frame(expr=x$expr,
                                                      time=x$time,
                                                      rep=1:length(x$expr))))

    data.frame(s=s, N=N, k=k, vals)
}


get_batch <- function(i, cases, reps=10) {
    plyr::ddply(cases, c("s","N","k"), run_bench, reps=reps) %>%
        mutate(batch=i)
}


reps <- 100
times <- ceiling(total/cores)

## times in milliseconds
cases <- expand.grid(s = 1000,
                     N = c(20,100,250,500), ##,200,300,400,500), ##, 500, 1000),
                     k = 3) %>%
    mutate(nvars=(N+1)*k,
           nels = nvars^2,
           nnz = N*k^2 + k^2 + 2*N*k*k,
           nnzLT = (N+1) * k*(k+1)/2 + N*k*k,
           pct.nnz = nnz/nels)

print(system.time(runtimes <- foreach(batch=1:cores, .combine=rbind) %dopar%
                      get_batch(batch, cases, reps=times)))


save(cases, runtimes, file="vignettes/runtimes.Rdata")


