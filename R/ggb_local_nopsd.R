#' Evaluate the GGB-Local Penalty Proximal Operator Along Grid of Lambda Values
#'
#'  Returns \code{Sig}, a (length \code{nlam}) list of \code{Matrix} objects,
#'  with \code{Sig[[l]]} being the solution for the l-th lambda value.
#'
#' @param S a symmetric, p by p matrix
#' @param g seed graph
#' @param lambda a vector of positive tuning parameters. If NULL, uses
#'        \code{nlam}, \code{flmin} to produce a grid
#' @param w weights (currently not implemented)
#' @param flmin ratio of smallest to largest lambda value (ignored if
#'        \code{lambda} is nonnull)
#' @param nlam number of lambda values (ignored if \code{lambda} is nonnull)
#' @param max_depths maximal bandwidths considered.  Should be a p-vector of
#'        non-negative integers. Default is NULL, which is equivalent to taking
#'        max_depths to be >= diameter(g).  (These are referred to as M_j in
#'        the paper.)
#' @param maxiter maximum number of passes through coordinates
#' @param tol convergence threshold
#' @param verbose level of verbosity in printed output (2, 1, 0)
ggb_local_nopsd <- function(S, g, lambda = NULL, w = NULL, flmin = 0.01,
                            nlam = 20, max_depths = NULL, maxiter = 500,
                            tol = 1e-4, verbose = 1) {
  if (!is.null(w)) stop("general w not implemented in c function yet.")
  D <- igraph::shortest.paths(g)
  D[D == Inf] <- -1
  M <- apply(D, 1, max)
  if (!is.null(max_depths)) {
    stopifnot(length(max_depths) == vcount(g),
              max_depths >= 0,
              max_depths == round(max_depths))
    M <- pmin(M, max_depths)
  }
  D[D == -1] <- max(M) + 1 # instead of Inf, we use 1 larger than all the M's.
  out <- prox_bcd(R = S, D = D, M = M, lambda = lambda, nlam = nlam,
                  flmin = flmin, maxiter = maxiter, tol = tol,
                  verbose = verbose)
  list(Sig = out$sumV, lambda = out$lambda, obj = out$obj)
}

#' R Wrapper to C Function for GGB-Local Penalty Proximal Operator
#'
#'  Returns sumvv which can be thought of as a (p choose 2)-by-nlam matrix,
#'  each column being the lower triangle of sum_jb V_jb.
#'
#'  argmin 0.5 * || R - sum_jb V_jb ||_F^2 + lam * sum_jb w_jb || V_jb ||_F
#'  where V_jb is 0 off of g_jb.
#'
#' @param R a symmetric, p by p matrix
#' @param D matrix of graph distances
#' @param M p-vector of maximal bandwidths
#' @param lambda a vector of positive tuning parameters. If NULL, uses
#'        \code{nlam}, \code{flmin} to produce a grid (formation of lambda grid
#'        in this case is done within C)
#' @param nlam number of lambda values (ignored if \code{lambda} is nonnull)
#' @param flmin ratio of smallest to largest lambda value (ignored if
#'        \code{lambda} is nonnull)
#' @param maxiter maximum number of passes through coordinates
#' @param tol convergence threshold
#' @param verbose level of verbosity in printed output (2, 1, 0)
#' @useDynLib ggb
prox_bcd <- function(R, D, M, lambda = NULL, nlam = 20, flmin = 0.01,
                     maxiter = 100, tol = 1e-4, verbose = 1) {
  # calls C function
  #
  # void pathwiseprox_bcd3(double *rr, int *dd, int *M, double *lamlist,
  # int *nlam, double *flmin, int *p, double *sumvv, double *obj, int *maxiter,
  # double *tol, int *verbose)
  #
  p <- nrow(R)
  stopifnot(is.matrix(D))
  if (is.null(lambda)) {
    lambda <- rep(0, nlam) # this will be ignored within C
    stopifnot(nlam > 0, flmin > 0, flmin < 1)
  } else {
    stopifnot(lambda > 0)
    nlam = length(lambda)
    flmin = -1 # this is used in C to determine that lambda is being passed
  }
  out <- .C("pathwiseprox_bcd3",
            rr = as.double(R[lower.tri(R)]),
            as.integer(D[lower.tri(D)]),
            as.integer(M),
            lambda = as.double(lambda),
            as.integer(nlam),
            as.double(flmin),
            as.integer(p),
            sumvv = as.double(rep(0, nlam * choose(p, 2))),
            obj = as.double(rep(0, nlam)),
            as.integer(maxiter),
            as.double(tol),
            as.integer(verbose),
            PACKAGE = "ggb")
  sumvv <- matrix(out$sumvv, ncol = nlam)
  sumV <- list()
  for (l in seq(nlam)) {
    sumV[[l]] <- lowertri2mat(sumvv[, l], p)
    diag(sumV[[l]]) <- diag(R)
  }
  list(sumV = sumV, obj = out$obj, lambda = out$lambda)
}

lowertri2mat <- function(rr, p) {
  # lower tri part to symmetric Matrix (diagonal remains zero)
  R <- Matrix::Matrix(0, p, p, sparse = TRUE)
  R[lower.tri(R)] <- rr
  return(Matrix::forceSymmetric(R, uplo = "L"))
}
