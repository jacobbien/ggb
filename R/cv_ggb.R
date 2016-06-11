#' Cross Validation to Select Lambda for GGB Estimator
#'
#' Performs nfold-fold cross validation to select the tuning parameter for GGB.
#' Takes \code{fit} object created by \code{\link{ggb}} and computes CV error at
#' the sequence of lambda values from that object.
#'
#' @param x n-by-p data matrix
#' @param fit object of class ggb (as outputted by \code{\link{ggb}})
#' @param g seed graph (should match that used to compute fit)
#' @param errfun a user-specified function measuring the loss incurred by
#' estimating \code{est} (first argument) when the true covariance matrix
#' is \code{true} (second argument).  Default: Squared Frobenius norm.
#' @param nfolds number of folds for cross-validation (default: 5)
#' @param ... additional arguments (passed to ggb)
#' @export
#' @return
#' \describe{
#' \item{errs: }{A \code{nlam}-by-\code{nfolds} matrix of errors.
#' \code{errs[l, j]} is error incurred in using \code{lamlist[l]} on fold \code{j}}
#' \item{m: }{CV error error for each value of lambda.}
#' \item{se: }{Standard error (estimated over folds) for each value of lambda}
#' \item{lambda_best: }{Value of \code{lambda} minimizing CV error.}
#' \item{ibest: }{Index of \code{lambda} minimizing CV error.}
#' \item{lambda_1se: }{Selected value of lambda using the one-standard-error
#' rule, a common heuristic that favors a sparser model when there isn't strong
#'  evidence against it.}
#' \item{i1se: }{Index of \code{lambda} of one-standard-error rule.}
#' }
#' @seealso \code{\link{ggb}}
#' @examples
#' set.seed(123)
#' n <- 20
#' p <- 10
#' g <- igraph::graph.lattice(c(5, 2))
#' x <- matrix(rnorm(n * p), n, p)
#' S <- stats::cov(x)
#' fit <- ggb(S, g, type = "local")
#' cv <- cv_ggb(x, fit, g)
#' plot(cv)
#' Sighat <- fit$Sig[[cv$i1se]]
cv_ggb <- function(x, fit, g, errfun = NULL, nfolds = 5, ...) {
  p <- ncol(x)
  n <- nrow(x)
  stopifnot(igraph::vcount(g) == p)
  stopifnot(class(fit) == "ggb")
  if (length(fit$Sig) == 0)
    stop("fit must have at least one lambda value.")
  stopifnot(nrow(fit$Sig[[1]]) == p)
  if (is.null(errfun)) {
    errfun <- function(est, true) sum(as.matrix(est) - as.matrix(true))^2
  }
  folds <- make_folds(n, nfolds)
  nlam <- length(fit$lambda)
  errs <- matrix(NA, nlam, nfolds)
  for (i in seq(nfolds)) {
    # train on all but i-th fold (and use settings from fit):
    fitcv <- ggb(stats::cov(x[-folds[[i]], ]), g = g, type = fit$type,
                 lambda = fit$lambda, delta=fit$delta, w = fit$w, ...)
    # evaluate this on left-out fold:
    Sig.te <- stats::cov(x[folds[[i]], ])
    for (l in seq(nlam)) errs[l, i] <- errfun(fit$Sig[[l]], Sig.te)
  }
  m <- rowMeans(errs)
  se <- apply(errs, 1, stats::sd) / sqrt(nfolds)
  ibest <- which.min(m)
  i1se <- min(which(m < m[ibest] + se[ibest]))
  out <- list(errs = errs, m = m, se = se, lambda_best = fit$lambda[ibest],
              ibest = ibest, lambda_1se = fit$lambda[i1se], i1se = i1se,
              lambda = fit$lambda,
              nonzero = unlist(lapply(fit$Sig, Matrix::nnzero)))
  class(out) <- "cv.ggb"
  out
}

#' Make folds for cross validation
#'
#' Divides the indices \code{1:n} into \code{nfolds} random folds of about the same size.
#'
#' @param n sample size
#' @param nfolds number of folds
make_folds <- function(n, nfolds) {
  nn <- round(n / nfolds)
  sizes <- rep(nn, nfolds)
  sizes[nfolds] <- sizes[nfolds] + n - nn * nfolds
  b <- c(0, cumsum(sizes))
  ii <- sample(n)
  folds <- list()
  for (i in seq(nfolds))
    folds[[i]] <- ii[seq(b[i] + 1, b[i + 1])]
  folds
}
