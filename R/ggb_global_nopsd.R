#' Compute the Global GGB estimator without psd constraint
#'
#' Uses the path BCD "closed form" of Yan and Bien (2015) as implemented
#' in R package \pkg{hsm}.
#'
#' @param S p-by-p sample covariance matrix
#' @param g seed graph (not necessarily connected)
#' @param lambda positive tuning parameter(s)
#' @param w weights from b = 1 to diam(g). Default is sqrt( group size )
#' @param flmin ratio of smallest to largest lambda value (ignored if
#'        \code{lambda} is nonnull)
#' @param nlam number of lambda values (ignored if \code{lambda} is nonnull)
#' @param max_depths maximal bandwidth considered.  Should be a single
#'        non-negative integer.  Default is NULL, which
#'        is equivalent to taking max_depths to be >= diameter(g).  (This is
#'        referred to as M in the paper.)
ggb_global_nopsd <- function(S, g, lambda = NULL, w = NULL, flmin = 0.01,
                             nlam = 20, max_depths = NULL) {
  p <- nrow(S)
  stopifnot(igraph::vcount(g) == p)

  D <- igraph::shortest.paths(g)
  D[D == Inf] <- 0 # Inf occurs when g is not connected.
  # Edges b/w components should not be assigned to any group
  D[D > max_depths] <- 0 # no pairs that are beyond max_depths from each other
  if (is.null(w)) {
    # default is sqrt of group size
    w <- as.numeric(sqrt(cumsum(tabulate(D[D != 0]))))
  }
  if (is.null(lambda)) {
    # produce smallest lambda that ensures Sig is diagonal
    lammax <- compute_lammax_global(S = S, g = g, max_depths = max_depths)
    lambda <- lammax * exp(seq(0, log(flmin), length = nlam))
  } else {
    nlam <- length(lambda)
  }
  obj <- rep(NA, nlam) # penval[l] is lambda*Omega(Sigma) for lth lambda
  Sig <- list()
  diagS <- diag(S)
  for (l in seq(nlam)) {
    fit <- pathprox(S, lambda[l], D, w)
    Sig[[l]] <- Matrix::Matrix(fit$b, p, p)
    diag(Sig[[l]]) <- diagS # diagonal had been zero
    obj[l] <- 0.5 * sum((S - as.matrix(Sig[[l]]))^2) + fit$pen
  }
  list(Sig = Sig, lambda = lambda, w = w, obj = obj)
}

#' R Wrapper to C Function
#'
#' Calls C function defined in R package \pkg{hsm}.
#'
#' @param y array of p elements
#' @param lam positive tuning parameter
#' @param assign array indexing each element in y with node indices (with
#'        0 meaning element is not in the graph)
#' @param w array of positive weights, one per node
#' @return b, w, and pen (value of lam * Omega(b) )
pathprox <- function(y, lam, assign, w = NULL) {
  # calls pathgraph_prox(double *y, double *lam, double *w, int *assign,
  #                      int *dim, int *depth)
  stopifnot(lam > 0, w > 0)
  p <- length(y)
  depth <- max(assign[is.finite(assign)])
  penal <- 0
  stopifnot(length(assign) == p)
  if (is.null(w)) {
    # wts are sqrt(num in grp):
    w <- sqrt(cumsum(table(assign[!(assign %in% c(0, Inf))])))
  }
  else
    stopifnot(length(w) == depth)
  out <- .C("pathgraph_prox2",
            r = as.double(y),
            as.double(lam),
            as.double(w),
            as.integer(assign),
            as.integer(p),
            as.integer(depth),
            penalval = as.double(penal),
            PACKAGE = "hsm")
  list(b = out$r, w = w, pen = out$penalval)
}
