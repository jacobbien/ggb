context("local without psd")

n <- 20
p <- 10
set.seed(123)
x <- matrix(rnorm(n * p), n, p)
S <- cov(x)
g <- igraph::graph.lattice(c(5, 2))

test_that("output beats Convex.jl", {
  lam <- 0.01
  #cvx <- ggb_local_convexjl(S, lam, g = g, delta = NULL)
  fit1 <- ggb::ggb(S, g, type = "local", lambda = lam, delta = NULL)
  #expect_true(fit1$obj < cvx$obj)
  expect_true(fit1$obj < 0.1367512)

  lam <- 0.1
  #cvx <- ggb_local_convexjl(S, lam, g = g, delta = NULL)
  fit2 <- ggb::ggb(S, g, type = "local", lambda = lam, delta = NULL)
  #expect_true(fit2$obj < cvx$obj)
  expect_true(fit2$obj < 0.9956542)

  lam <- 0.2
  #cvx <- ggb_local_convexjl(S, lam, g = g, delta = NULL)
  fit3 <- ggb::ggb(S, g, type = "local", lambda = lam, delta = NULL)
  #expect_true(fit3$obj < cvx$obj)
  expect_true(fit3$obj < 1.329507)

})

test_that("largest lambda gives diagonal matrix", {
  fit <- ggb::ggb(S, g, type = "local", nlam = 2, delta = NULL)
  expect_true(Matrix::isDiagonal(fit$Sig[[1]]))
})

test_that("solving pathwise gives same results as one-at-a-time", {
  fit <- ggb::ggb(S, g, type = "local", delta = NULL, in_tol = 1e-5)
  Sig2 <- list()
  for (l in seq_along(fit$lambda)) {
    Sig2[[l]] <- ggb::ggb(S, g, type = "local", lambda = fit$lambda[l],
                          delta = NULL, in_tol = 1e-5)$Sig[[1]]
  }
  expect_equal(fit$Sig, Sig2, tolerance = 1e-4)
})
