## Part of the sparseMVN package
## Copyright (C) 2013-2017 Michael Braun

## Functions to compute objective function, gradient, and Hessian for
## binary choice example, and some unit tests.


#' @name binary
#' @title Binary choice example
#' @description Functions for binary choice example in the vignette.
#' @param P Numeric vector of length \eqn{(N+1)k}.  First \eqn{Nk}
#' elements are heterogeneous coefficients. The remaining k elements are population parameters.
#' @param data Named list of data matrices Y and X, and choice count integer T
#' @param priors Named list of matrices inv.Omega and inv.A.
#' @param order.row Determines order of heterogeneous coefficients in
#' parameter vector. If TRUE, heterogeneous coefficients are ordered by unit.  If FALSE, they are ordered by covariate.
#' @return For binary.f, binary.df and binary.hess, the log posterior density, gradient and Hessian, respectively. The Hessian is a dgCMatrix object. binary.sim returns a list with simulated Y and X, and the input T.
#' @details These functions are used by the heterogeneous binary
#' choice example in the vignette. There are N heterogeneous units, each making T binary choices.  The choice probabilities depend on k covariates.  binary.sim simulates a dataset suitable for running the example.
#' @rdname binary
#' @export
binary.f <- function(P, data, priors, order.row=FALSE) {

    N <- length(data$Y)
    k <- NROW(data$X)
    beta <- matrix(P[1:(N*k)], k, N, byrow=order.row)
    mu <- P[(N*k+1):(N*k+k)]

    bx <- colSums(data$X * beta)

    ## can't use log1p if P is complex
    log.p <- bx - log(1+exp(bx))
    log.p1 <- -log(1+exp(bx))

    data.LL <- sum(data$Y*log.p + (data$T-data$Y)*log.p1)

    Bmu <- apply(beta, 2, "-", mu)

    prior <- -0.5 * sum(diag(tcrossprod(Bmu) %*% priors$inv.A))
    hyp <- -0.5 * t(mu) %*% priors$inv.Omega %*% mu
    res <- data.LL + prior + hyp
    if(is.complex(P)) return(as.complex(res)) else return(as.numeric(res))

}

#' @rdname binary
#' @export
binary.grad <- function(P, data, priors, order.row=FALSE) {

    Y <- data$Y
    X <- data$X
    T <- data$T
    inv.A <- priors$inv.A
    inv.Omega <- priors$inv.Omega

  q1 <- .dlog.f.db(P, Y, X, T, inv.Omega, inv.A, order.row=order.row)
  q2 <- .dlog.f.dmu(P, Y, X, T, inv.Omega, inv.A, order.row=order.row)
  res <- c(q1, q2)
  return(res)
}

#' @rdname binary
#' @export
binary.hess <- function(P, data, priors, order.row=FALSE) {

    Y <- data$Y
    X <- data$X
    T <- data$T
    inv.A <- priors$inv.A
    inv.Omega <- priors$inv.Omega
    N <- length(Y)
    k <- NROW(X)

    SX <- Matrix(inv.A)
    XO <- Matrix(inv.Omega)
    B1 <- .d2.db(P, Y, X, T, SX, order.row)

    if(order.row) {
        pvec <- NULL
        for (i in 1:N) {
            pvec <- c(pvec, seq(i, i+(k-1)*N, length=k))
        }

        pmat <- as(pvec, "pMatrix")
        B2 <- Matrix::t(pmat) %*% B1 %*% pmat
        cross <- .d2.cross(N, SX) %*% pmat
    } else {
        B2 <- B1
        cross <- .d2.cross(N, SX)
    }

    Bmu <- .d2.dmu(N,SX, XO)
    res <- rbind(cbind(B2, Matrix::t(cross)),cbind(cross, Bmu))

    return(drop0(res))
}

#' @param N Number of heterogeneous units
#' @param k Number of heterogeneous parameters
#' @param T Observations per household
#' @rdname binary
#' @export
binary.sim <- function(N, k, T) {

    x.mean <- rep(0,k)
    x.cov <- diag(k)
    x.cov[1,1] <- .02
    x.cov[k,k] <- x.cov[1,1]
    mu <- seq(-2,2,length=k)
    Omega <- diag(k)

    X <- t(mvtnorm::rmvnorm(N, mean=x.mean, sigma=x.cov)) ## k x N
    B <- t(mvtnorm::rmvnorm(N, mean=mu, sigma=Omega)) ## k x N
    XB <- colSums(X * B)
    log.p <- XB - log1p(exp(XB))
    Y <- sapply(log.p, function(q) return(stats::rbinom(1,T,exp(q))))

    binary <- list(Y=Y, X=X, T=T)
    return(binary)
}




.dlog.f.db <- function(P, Y, X, T, inv.Omega, inv.A, order.row) {

    N <- length(Y)
    k <- NROW(X)


    beta <- matrix(P[1:(N*k)], k, N, byrow=order.row)
    mu <- P[(N*k+1):length(P)]
    bx <- colSums(X * beta)

    p <- exp(bx)/(1+exp(bx))

    tmp <- Y - T*p

    dLL.db <- apply(X,1,"*",tmp)

    Bmu <- apply(beta, 2, "-", mu)
    dprior <- -inv.A %*% Bmu

    res <- t(dLL.db) + dprior
    if (order.row) res <- t(res)

    return(as.vector(res))

}

.dlog.f.dmu <- function(P, Y, X, T, inv.Omega, inv.A, order.row) {

    N <- length(Y)
    k <- NROW(X)

    beta <- matrix(P[1:(N*k)], k, N, byrow=order.row)
    mu <- P[(N*k+1):length(P)]
    Bmu <- apply(beta, 2, "-", mu)

    res <- inv.A %*% (rowSums(Bmu)) -  inv.Omega %*% mu
    return(res)
}

.d2.db <- function(P, Y, X, T, inv.A, order.row) {

    N <- length(Y)
    k <- NROW(X)

    beta <- matrix(P[1:(N*k)], k, N, byrow=order.row)
    mu <- P[(N*k+1):length(P)]
    ebx <- exp(colSums(X * beta))

    p <- ebx/(1+ebx)

    q <- vector("list",length=N)
    for (i in 1:N) {
        q[[i]] <- -T*p[i]*(1-p[i])*tcrossprod(X[,i]) - inv.A
    }
    B <- bdiag(q)
    return(B)
}

.d2.dmu <- function(N, inv.A, inv.Omega) {
  return(-N*inv.A-inv.Omega)
}

.d2.cross <- function(N, inv.A) {
  res <- kronecker(matrix(rep(1,N),nrow=1),inv.A)
  return(res)
}




