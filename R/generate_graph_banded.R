#' Generate a Graph-Banded Covariance Matrix
#'
#' Returns a p-by-p matrix covariance matrix that is banded with respect to the
#' graph \code{g} and with local bandwidths given in the p-vector \code{b}.
#'
#' @param g a seed graph (with p vertices)
#' @param b a p-vector of local bandwidths
#' @param sigma sd of independent error to add to each node
#' @param cor whether to scale rows/cols to have variance 1
#' @export
#' @examples
#' g <- igraph::graph.lattice(c(5, 4))
#' b <- rep(1, 20)
#' Sig <- generate_gb_covariance(g, b)
generate_gb_covariance <- function(g, b, sigma = 0.01, cor = TRUE) {
  p <- length(b)
  stopifnot(igraph::vcount(g) == p)
  stopifnot(b >= 0, b == round(b))
  D <- igraph::shortest.paths(g)
  b <- pmin(b, apply(D, 1, max))
  A <- D <= b # A_jk indicates whether k is in N_j
  B <- A | t(A)
  Sig <- B / D
  diag(Sig) <- 0
  # add to diag so min eval=sig^2
  diag(Sig) <- max(-min(eigen(Sig)$val), 0) + sigma^2
  if (cor) {
    sig <- sqrt(diag(Sig))
    D <- diag(1 / sig)
    Sig <- D %*% Sig %*% D
  }
  Sig
}
