
#' Print output for a GGB object
#'
#' @param x ggb object as returned by \code{\link{ggb}}
#' @param ... additional arguments (not used)
#'
#' @export
print.ggb <- function(x, ...) {
  cat(sprintf("GGB object of type %s with %s lambda values.", x$type,
              length(x$lambda)), fill = TRUE)
  tab <- data.frame(lambda = x$lambda,
                    "num nonzeros" = unlist(lapply(x$Sig, Matrix::nnzero)))
  print(tab, row.names = FALSE)
}
