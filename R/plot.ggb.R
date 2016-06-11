#' Plot a GGB object
#'
#' Uses image function from the \pkg{Matrix} R package.
#'
#' @param x ggb object as returned by \code{\link{ggb}}
#' @param subset which lambda values to plot, expressed as a subset of indices,
#'        between 1 and \code{length(x$Sig)}
#' @param ... additional arguments to pass to \code{Matrix::image}
#' @export
plot.ggb <- function(x, subset = seq_along(x$Sig), ...) {
  num_sig <- length(x$Sig)
  stopifnot(subset %in% seq_along(x$Sig))
  num_sig <- length(subset)
  nrow <- floor(sqrt(num_sig))
  ncol <- ceiling(num_sig / nrow)
  a <- lapply(x$Sig[subset], image_covariance)
  for (i in seq_along(subset)) {
    a[[i]]$main <- list(label = as.character(round(x$lambda[subset[i]], 5)),
                        cex = 0.75)
    a[[i]]$x.scales$draw <- FALSE
    a[[i]]$y.scales$draw <- FALSE
    a[[i]]$par.settings$layout.widths=list(left.padding = -2,
                                           right.padding = -2)
    a[[i]]$par.settings$layout.heights=list(top.padding = 1,
                                           bottom.padding = 0)
  }
  if (length(a) == 1) return(print(a[[1]]))
  for (i in seq(length(subset) - 1)) {
    print(a[[i]], split = c(icolrow(i, ncol), ncol, nrow), more = TRUE, ...)
  }
  print(a[[i + 1]], split = c(icolrow(i + 1, ncol), ncol, nrow), ...)
}

icolrow <- function(i, ncol) {
  # allows for par(mfrow)-like behavior from lattice's split
  irow <- floor((i  - 1) / ncol)
  icol <- (i - 1) - ncol * irow
  1 + c(icol, irow)
}

#' Show the Image of a Covariance Matrix
#'
#' Uses the \code{image} function defined in \pkg{Matrix} package.
#' @param Sig covariance matrix
#' @param sub subtitle, default NULL
#' @param xlab x label, default NULL
#' @param ylab y label, default NULL
#' @param ... additional arguments to pass to \code{Matrix::image}
#' @export
image_covariance <- function(Sig, sub = NULL, xlab = NULL, ylab = NULL, ...) {
  Matrix::image(Matrix::Matrix(Sig), sub = sub, xlab = xlab, ylab = ylab, ...)
}
