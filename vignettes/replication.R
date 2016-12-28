library(plyr)
library(tidyverse)
library(Matrix)
library(mvtnorm)
library(microbenchmark)
library(doParallel)
registerDoParallel(cores=4)

get_times <- function(D, reps=100) {

    s <- D$s ## number of random samples
    p <-  D$p ## population variables
    k <- D$k ## heterogeneous variables
    N <- D$N ## number of agents
    prec <- D$prec ## 

    mu <- rep(0,k*N + p)  ## assume mean at origin
    Q1 <- tril(kronecker(Matrix(seq(0.1,k,length=k*k),k,k),diag(N)))
    Q2 <- cBind(Q1,Matrix(0, N*k, p))
    Q3 <- rBind(Q2,cBind(Matrix(rnorm(N*k*p), p, N*k), Diagonal(k)))
    CV.sparse <- tcrossprod(Q3)
    CV.dense <- as(CV.sparse, "matrix")  ## dense covariance

    if (prec) sigma <- solve(CV.dense) else sigma <- CV.dense

    chol.CV <- Cholesky(CV.sparse)
    x.sp <- rmvn.sparse(s, mu, chol.CV, prec=prec)
    d.sp <- dmvn.sparse(x.sp, mu, chol.CV, prec=prec)
    d.dens <- dmvnorm(x.sp, mu, sigma, log=TRUE)
    stopifnot(all.equal(d.sp,d.dens))
    
    bench <- microbenchmark(
        chol = Cholesky(CV.sparse),
        solve = solve(CV.dense),
        r_sparse = rmvn.sparse(s, mu, chol.CV, prec=prec),
        d_sparse = dmvn.sparse(x.sp, mu, chol.CV, prec=prec),
        r_dense = rmvnorm(s, mu, sigma, method="chol"),
        d_dense = dmvnorm(x.sp, mu, sigma, log=TRUE),
        times = reps
    )

    vals <- plyr::ddply(data.frame(bench), "expr",
                        function(x) return(data.frame(expr=x$expr,
                                                      time=x$time,
                                                      rep=1:length(x$expr))))
    
    data.frame(s=s, N=N, k=k, p=p, prec=prec, bench=vals)
}


reps <- 10
## times in milliseconds
runtimes <- expand.grid(s = c(20, 100),
                 N = c(5, 10, 15),
                 k = 3,
                 p = 3,
                 prec = c(FALSE, TRUE)) %>%
    ddply(c("s","N","k","p","prec"), get_times, reps=reps,
          .parallel=FALSE) %>%
    group_by(s, N, k, p, prec, bench.expr) %>%
    summarize(mean = mean(bench.time/10^6),
              sd = sd(bench.time/10^6)) %>%
    rename(step=bench.expr) %>%
    mutate(nvars=N*k+p,
           nels = nvars^2,
           nnz = N*k^2 + p^2 + 2*N*k*p,
           nnzLT = N * k*(k+1)/2 + p*(p+1)/2 + p*N*k,
           pct.nnz = nnz/nels) %>%
    ungroup()

save(runtimes, file="vignettes/runtimes.Rdata")


