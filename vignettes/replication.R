library(plyr)
library(tidyr)
library(stringr)
library(dplyr)
library(Matrix)
library(mvtnorm)
library(microbenchmark)
library(doParallel)
library(sparseMVN)

registerDoParallel(cores=8)

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
    
    data.frame(s=s, N=N, k=k, prec=prec, bench=vals)
}


reps <- 100
## times in milliseconds
cases <- expand.grid(s = 1000,
                     N = c(25, 250, 1000),
                     k = 2,
                     prec = c(FALSE, TRUE)) %>%
    mutate(nvars=(N+1)*k,
           nels = nvars^2,
           nnz = N*k^2 + k^2 + 2*N*k*k,
           nnzLT = (N+1) * k*(k+1)/2 + N*k*k,
           pct.nnz = nnz/nels) 

runtimes <- ddply(cases, c("s","N","k","prec"), get_times, reps=reps,
                  .parallel=TRUE) %>%
    group_by(s, N, k, prec, bench.expr) %>%
    summarize(mean = mean(bench.time/10^6),
              sd = sd(bench.time/10^6)) %>%
    ungroup() %>%
    tidyr::separate(bench.expr, c("step", "storage"), sep="_")
    

save(cases, runtimes, file="vignettes/runtimes.Rdata")


