context("global with no psd")

n <- 20
p <- 10
set.seed(123)
x <- matrix(rnorm(n * p), n, p)
S <- cov(x)
g <- igraph::graph.lattice(c(5, 2))
tol <- 1e-4

test_that("output beats Convex.jl", {
  lam <- 0.01
  #cvx <- ggb_global_convexjl(S, lam, g = g, delta = NULL)
  fit1 <- ggb::ggb(S, g, type = "global", lambda = lam, delta = NULL)
  #expect_true(fit1$obj < cvx$obj)
  expect_true(fit1$obj < 0.151305)

  lam <- 0.1
  #cvx <- ggb_global_convexjl(S, lam, g = g, delta = NULL)
  fit2 <- ggb::ggb(S, g, type = "global", lambda = lam, delta = NULL)
  #expect_true(fit2$obj < cvx$obj)
  expect_true(fit2$obj < 1.10802)

  lam <- 0.2
  #cvx <- ggb_global_convexjl(S, lam, g = g, delta = NULL)
  fit3 <- ggb::ggb(S, g, type = "global", lambda = lam, delta = NULL)
  #expect_true(fit3$obj < cvx$obj + tol)
  expect_true(fit3$obj < 1.36565 + tol)

})

test_that("largest lambda gives diagonal matrix", {
  fit <- ggb::ggb(S, g, type = "global", nlam = 2, delta = NULL)
  #fit$Sig[[1]][abs(fit$Sig[[1]]) < 1e-16] <- 0
  expect_true(Matrix::isDiagonal(fit$Sig[[1]]))

  # example from vignette:
  set.seed(123)
  gg <- igraph::graph.lattice(c(5, 4))
  pp <- 5 * 4
  bb <- rep(1, pp)
  Sigg <- generate_gb_covariance(gg, bb)
  nn <- 30
  eigg <- eigen(Sigg)
  AA <- diag(sqrt(eigg$values)) %*% t(eigg$vectors)
  xx <- matrix(rnorm(nn * pp), nn, pp) %*% AA
  fitt <- ggb::ggb(S = cov(xx), g = gg, type = "global", nlam = 2, delta = NULL)
  expect_true(Matrix::isDiagonal(fitt$Sig[[1]]))
})

test_that("solving pathwise gives same results as one-at-a-time", {
  fit <- ggb::ggb(S, g, type = "global", delta = NULL)
  Sig2 <- list()
  for (l in seq_along(fit$lambda)) {
    Sig2[[l]] <- ggb::ggb(S, g, type = "global", lambda = fit$lambda[l],
                          delta = NULL)$Sig[[1]]
  }
  expect_equal(fit$Sig, Sig2, tolerance = 1e-4)
})
