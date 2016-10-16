#' ggb: Graph-Guided Banding for Covariance Estimation
#'
#' This package implements the two methods introduced in
#' Bien, J. (2016) "Graph-Guided Banding of the Covariance Matrix"
#' \url{http://arxiv.org/abs/1606.00451}
#'
#' The main function is \code{\link{ggb}}.
#'
#' @docType package
#' @import hsm
#' @name ggb
NULL

#' Graph-Guided Banding Estimators
#'
#' Computes the global or local GGB estimator (depending on \code{type}
#' argument) along a grid of lambda values.
#'
#' To evaluate the proximal operator, it uses the "closed form" of Yan and Bien
#' (2015) for the global GGB estimator as implemented in the R package \pkg{hsm}
#' and uses blockwise coordinate descent for the local GGB estimator.
#'
#' If \code{delta} is non-NULL, then alternates between evaluating proximal
#' operator and projecting eigenvalues.
#'
#' The weights used in this penalty are the square root of the group size.
#'
#' @param S p-by-p sample covariance matrix
#' @param g seed graph (not necessarily connected)
#' @param type either "global" or "local"
#' @param lambda positive tuning parameter(s)
#' @param delta lower bound on eigenvalues (default, NULL, means no bound)
#' @param flmin ratio of smallest to largest lambda value (ignored if
#'        \code{lambda} is nonnull)
#' @param nlam number of lambda values (ignored if \code{lambda} is nonnull)
#' @param max_depths maximal bandwidth(s) considered.  If type is "global",
#'        this must be a single non-negative integer.  If type is "local", this
#'        must be a p-vector of non-negative integers. Default is NULL, which
#'        is equivalent to taking max_depths to be >= diameter(g).  (These are
#'        referred to as M and M_j in the paper.)
#' @param out_iter number of iterations of outer loop (i.e., # of eigenvalue
#'        decompositions). Ignored when \code{delta} is NULL.
#' @param in_iter number of iterations on inner loop (i.e., # cycles over
#'        row/cols) when \code{type}
#'        is "local".
#' @param out_tol convergence threshold for outer loop BCD (ignored when
#'        \code{delta} is NULL)
#' @param in_tol convergence threshold for inner loop BCD when \code{type}
#'        is "local".
#' @param verbose level of verbosity in printed output (3, 2, 1, 0)
#' @export
#' @seealso \code{\link{cv_ggb}}
#' @examples
#' set.seed(123)
#' n <- 100
#' p <- 20
#' g <- igraph::graph.lattice(c(5, 4))
#' x <- matrix(rnorm(n * p), n, p)
#' S <- cov(x)
#' fit_global <- ggb(S, g, type = "global", nlam = 10)
#' fit_local <- ggb(S, g, type = "local", nlam = 10)
ggb <- function(S, g, type = c("global", "local"), lambda = NULL,
                delta = NULL, flmin = 0.01, nlam = 20, max_depths = NULL,
                out_iter = 100, in_iter = 500, out_tol = 1e-4, in_tol = 1e-4,
                verbose = 0) {
  p <- nrow(S)
  stopifnot(igraph::vcount(g) == p)
  type <- type[1]
  if (!(type %in% c("global", "local")))
    stop("type must be either 'global' or 'local'.")
  if (!is.null(max_depths)) {
    stopifnot(max_depths >= 0, max_depths == round(max_depths))
    if (type == "global")
      if(length(max_depths) != 1)
        stop("max_depths must be a scalar when type is 'global'.")
    else if (type == "local")
      if (length(max_depths) != p)
        stop("max_depths must be of length p when type is 'local'.")
  }
  if (is.null(delta)) {
    if (type == "global")
      out <- ggb_global_nopsd(S, g, lambda = lambda, flmin = flmin,
                              nlam = nlam, max_depths = max_depths)
    else if (type == "local")
      out <- ggb_local_nopsd(S, g, lambda = lambda, flmin = flmin, nlam = nlam,
                             max_depths = max_depths,
                             maxiter = in_iter, tol = in_tol,
                             verbose = verbose - 1)
    else stop("invalid type")
    out <- list(Sig = out$Sig, lambda = out$lambda, type = type, delta = delta,
                max_depths = max_depths, obj = out$obj)
    class(out) <- "ggb"
    return(out)
  }
  # if delta is non-NULL, then we perform BCD as described in paper,
  # alternating proximal and eigenvalue steps.
  if (is.null(lambda)) {
    # find a lambda that is large enough to ensure all V are 0
    # (w/o psd constraint)
    if (type == "global")
      lammax <- compute_lammax_global(S = S, g = g, max_depths = max_depths)
    else if (type == "local")
      lammax <- compute_lammax_local(S = S, g = g, max_depths = max_depths)
    else stop("invalid type")
    lambda <- lammax * exp(seq(0, log(flmin), length = nlam))
  } else {
    nlam <- length(lambda)
  }
  Sig <- list()
  obj <- rep(NA, nlam)
  C <- matrix(0, p, p) # dual var corresponding to eigen constraint
  for (l in seq(nlam)) {
    if (verbose > 0) cat(sprintf("lam[%s] = %s", l, lambda[l]), fill = TRUE)
    for (i in seq(out_iter)) {
      # update over B: solve Prox on S + C
      R <- S + C
      if (type == "global") {
        out <- ggb_global_nopsd(R, g, lambda = lambda[l],
                                max_depths = max_depths)
      } else if (type == "local") {
        out <- ggb_local_nopsd(R, g, lambda = lambda[l], maxiter = in_iter,
                               tol = in_tol, verbose = verbose - 1,
                               max_depths = max_depths)
      } else stop("invalid type.")
      B <- R - out$Sig[[1]]
      if (verbose > 0)
        cat("  (min eval is", # out$Sig[[1]] is same as S + C - B
            min(eigen(out$Sig[[1]], only.values = TRUE)$val), ")",
            fill = TRUE)
      if (i > 1) if (max(abs(out$Sig[[1]] - oldSig)) < out_tol) {
        if (verbose > 0)
          cat("BCD converged in", i, "iterations.", fill = TRUE)
        break
      }
      oldSig <- out$Sig[[1]]
      # update over C: adjust eigenvalues
      eig <- eigen(B - S)
      C <- eig$vec %*% diag(pmax(eig$val + delta, 0)) %*% t(eig$vec)
      # keep using same w without recomputing
      w <- out$w
    }
    Sig[[l]] <- out$Sig[[1]] #S + C - B
    obj[l] <- out$obj + (sum((S - Sig[[l]])^2) - sum((R - Sig[[l]])^2))/2
  }
  out <- list(Sig = Sig, lambda = lambda, obj = obj)
  class(out) <- "ggb"
  out
}

compute_lammax_global <- function(S, g, max_depths = NULL) {
  if (!is.null(max_depths)) stopifnot(length(max_depths) == 1)
  D <- igraph::shortest.paths(g)
  D[is.infinite(D)] <- 0
  depth <- max(D) # diameter of largest connected component
  if (!is.null(max_depths)) depth <- min(depth, max_depths)
  m <- rep(NA, depth)
  for (d in seq(depth)) {
    m[d] <- mean(S[D <= d & D > 0]^2)
  }
  sqrt(max(m))
}


compute_lammax_local <- function(S, g, max_depths = NULL) {
  D <- igraph::shortest.paths(g)
  D[is.infinite(D)] <- 0
  depths <- apply(D, 1, max)
  if (!is.null(max_depths)) depths <- pmin(depths, max_depths)
  tot <- sum(depths)
  m <- rep(NA, tot)
  i <- 1
  for (j in seq(nrow(S))) {
    if (depths[j] == 0) next
    for (d in seq(depths[j])) {
      m[i] <- mean(S[j, D[j, ] <= d & D[j, ] > 0]^2)
      i <- i + 1
    }
  }
  sqrt(max(m))
}
