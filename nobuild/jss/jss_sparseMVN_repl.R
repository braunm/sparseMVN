library("sparseMVN")
library("microbenchmark")
library("Matrix")
library("mvtnorm")
library("dplyr")
library("tidyr")
library("stringr")
library("reshape2")
set.seed(123)


build_mat <- function(N, k) {
    t1 <- exp(rnorm(k*k))
    Q1 <- tril(kronecker(diag(N),Matrix(t1,k,k)))
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
        rand_dense_cov = rmvnorm(s, mu, CV.dense, method="chol"),
        density_dense_cov = dmvnorm(x, mu, CV.dense, log=TRUE),
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


cores <- 10
reps <- 200
times <- ceiling(reps/cores)


## times in milliseconds
cases <- expand.grid(s = 100,#0,
                     N = c(10, 20), # 50, 100, 200, 300, 400, 500),
                     k = c(2,4)) %>%
    mutate(nvars=(N+1)*k,
           nels = nvars^2,
           nnz = N*k^2 + k^2 + 2*N*k*k,
           nnzLT = (N+1) * k*(k+1)/2 + N*k*k,
           pct.nnz = nnz/nels)

RT <- plyr::ldply(1:cores, get_batch, cases, reps=times)

labs <- str_split_fixed(RT[['expr']],"_",3)
colnames(labs) <- c("stat","pattern","type")
runtimes <- cbind(RT, labs)

