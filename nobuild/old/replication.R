library("sparseMVN")
library("microbenchmark")
library("Matrix")
library("mvtnorm")
library("dplyr")

##library(tidyr)
library(doParallel)


registerDoParallel(cores=12)

get_times <- function(D, reps=100) {

    s <- D$s ## number of random samples
    k <- D$k ## heterogeneous variables
    N <- D$N ## number of agents
    prec <- D$prec ##

    cat("s =  ",s,"\tN = ",N,"\tk = ",k,"\tprec = ",prec,"\n")

    mu <- rep(0,k*N + k)  ## assume mean at origin
    Q1 <- tril(kronecker(Matrix(seq(0.1,k,length=k*k),k,k),diag(N)))
    Q2 <- cBind(Q1,Matrix(0, N*k, k))
    Q3 <- rBind(Q2,cBind(Matrix(rnorm(N*k*k), k, N*k), Diagonal(k)))
    CV.sparse <- tcrossprod(Q3)
    CV.dense <- as(CV.sparse, "matrix")  ## dense covariance

    if (prec) sigma <- solve(CV.dense) else sigma <- CV.dense

    chol.CV <- Cholesky(CV.sparse)
    x.sp <- rmvn.sparse(s, mu, chol.CV, prec=prec)
    d.sp <- dmvn.sparse(x.sp, mu, chol.CV, prec=prec)
    d.dens <- dmvnorm(x.sp, mu, sigma, log=TRUE)
    stopifnot(all.equal(d.sp,d.dens))

    bench <- microbenchmark(
        chol_sparse = Cholesky(CV.sparse),
        chol_dense = chol(CV.dense),
        solve_sparse = solve(CV.sparse),
        solve_dense = solve(CV.dense),
        rand_sparse = rmvn.sparse(s, mu, chol.CV, prec=prec),
        density_sparse = dmvn.sparse(x.sp, mu, chol.CV, prec=prec),
        rand_dense = rmvnorm(s, mu, sigma, method="chol"),
        density_dense = dmvnorm(x.sp, mu, sigma, log=TRUE),
        times = reps
    )

    vals <- plyr::ddply(data.frame(bench), "expr",
                        function(x) return(data.frame(expr=x$expr,
                                                      time=x$time,
                                                      rep=1:length(x$expr))))

    data_frame(s=s, N=N, k=k, prec=prec, bench=vals)
}


reps <- 20
## times in milliseconds
cases <- expand.grid(s = 250,
                     N = c(25, 100, 250, 500, 1000),
                     k = c(2,4),
                     prec = c(FALSE, TRUE)) %>%
    mutate(nvars=(N+1)*k,
           nels = nvars^2,
           nnz = N*k^2 + k^2 + 2*N*k*k,
           nnzLT = (N+1) * k*(k+1)/2 + N*k*k,
           pct.nnz = nnz/nels)

runtimes <- plyr::ddply(cases, c("s","N","k","prec"), get_times, reps=reps,
                  .parallel=TRUE)




##save(cases, runtimes, file="vignettes/runtimes.Rdata")


