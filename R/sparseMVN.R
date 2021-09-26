#' @name sparseMVN-package
#' @aliases sparseMVN-package sparseMVN
#' @docType package
#' @title Multivariate Normal Functions for Sparse Covariate and Precision Matrices
#' @description MVN functions for sparse covariance and precision matrices.
#' @details Computes multivariate normal (MVN) densities, and samples
#' from MVN distributions, when either the covariance or precision
#' matrix is stored as a sparse Matrix (a dsCMatrix object, as defined
#' in the Matrix package.  The user can provide the precision matrix
#' directly, rather than convert it to a covariance via matrix
#' inversion.
#' @import Matrix
#' @import methods
#' @importFrom stats rnorm
#' @encoding UTF-8
#' @keywords package
"_PACKAGE"
