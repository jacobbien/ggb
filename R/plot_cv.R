#' Plot a GGB-CV object
#'
#' This function is adapted from the \pkg{hierNet} R package.
#'
#' @param x an object of class \code{cv.ggb}, as returned by
#'        \code{\link{cv_ggb}}
#' @param ... additional arguments (not used)
#' @export
plot.cv.ggb <- function(x, ...) {
  graphics::par(mar = c(5, 5, 5, 1))
  yrang = range(c(x$m - x$se, x$m + x$se))
  graphics::plot(log(x$lambda), x$m, xlab = "log(lambda)",
                 ylab = "Cross-validation Error",
       type = "n", ylim = yrang)
  graphics::axis(3, at = log(x$lambda), labels = paste(x$nonzero), srt = 90,
       adj = 0)
  graphics::mtext("Number of nonzeros", 3, 4, cex = 1.2)
  graphics::axis(2, at = c(0, 0.2, 0.4, 0.6, 0.8))
  error_bars(log(x$lambda), x$m - x$se, x$m + x$se, width = 0.01,
             col = "darkgrey")
  graphics::points(log(x$lambda), x$m, col = 2, pch = 19)
  graphics::abline(v = log(x$lambda_best), lty = 3)
  graphics::abline(v = log(x$lambda_1se), lty = 3)
  invisible()
}

error_bars <- function (x, upper, lower, width = 0.02, ...) {
  # this function came from hierNet package
  xlim <- range(x)
  barw <- diff(xlim) * width
  graphics::segments(x, upper, x, lower, ...)
  graphics::segments(x - barw, upper, x + barw, upper, ...)
  graphics::segments(x - barw, lower, x + barw, lower, ...)
  range(upper, lower)
}
