#' @rdname dmvn.sparse
#' @inherit rmvn.sparse
#' @export
dmvn.sparse <- function(x, mu, CH, prec=TRUE, log=TRUE) {

    if (is.vector(x) | (is.atomic(x) & NCOL(x)==1)) {
        x <- matrix(x,nrow=1)
    }

    k <- length(mu)
    n <- NROW(x)
    if (!(k>0)) {
        stop("mu must have positive length")
    }

    if (!(k==dim(CH)[1])) {
        stop("dimensions of mu and CH do not conform")
    }
    if (k!=NCOL(x)) {
        stop("x must have same number of columns as the length of mu")
    }
    if (!is.logical(prec)) {
        stop("prec must be either TRUE or FALSE")
    }
    A <- expand(CH)
    detL <- sum(log(Matrix::diag(A$L)))
    C <- -0.918938533204672669541*k ## -k*log(2*pi)/2
    xmu <- t(x)-mu
    z <- as.matrix(A$P %*% xmu)

    if (prec) {
        y <- Matrix::crossprod(A$L,z)  ## L' %*% x
        log.dens <- C + detL - Matrix::colSums(y*y)/2
    } else {
        y <- solve(A$L, z) ## Ly = x
        log.dens <- C - detL - Matrix::colSums(y*y)/2
    }

    if (log) return (log.dens) else return (exp(log.dens))
}
