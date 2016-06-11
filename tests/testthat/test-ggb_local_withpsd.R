context("local with psd")

n <- 20
p <- 10
set.seed(123)
x <- matrix(rnorm(n * p), n, p)
S <- cov(x)
g <- igraph::graph.lattice(c(5, 2))
tol <- 1e-7

test_that("output beats Convex.jl", {
  lam <- 0.01
  del <- 0.2
  #cvx <- ggb_local_convexjl(S, lam, g = g, delta = del)
  fit1 <- ggb::ggb(S, g, type = "local", lambda = lam, delta = del,
                   out_tol = 1e-7)
  #expect_true(fit1$obj < cvx$obj)
  expect_true(fit1$obj < 0.1390355)
  expect_true(min(eigen(fit1$Sig[[1]], only.values = TRUE)$values) > del - tol)

  lam <- 0.1
  del <- 0.4
  #cvx <- ggb_local_convexjl(S, lam, g = g, delta = del)
  fit2 <- ggb::ggb(S, g, type = "local", lambda = lam, delta = del,
                   out_tol = 1e-7)
  #expect_true(fit2$obj < cvx$obj)
  expect_true(fit2$obj < 0.9956679)
  expect_true(min(eigen(fit2$Sig[[1]], only.values = TRUE)$values) > del - tol)

  lam <- 0.2
  del <- 0.6
  #cvx <- ggb_local_convexjl(S, lam, g = g, delta = del)
  fit3 <- ggb::ggb(S, g, type = "local", lambda = lam, delta = del,
                   out_tol = 1e-7)
  #expect_true(fit3$obj < cvx$obj)
  expect_true(fit3$obj < 1.335467)
  expect_true(min(eigen(fit3$Sig[[1]], only.values = TRUE)$values) > del - tol)

})

test_that("largest lambda gives diagonal matrix", {
  fit <- ggb::ggb(S, g, type = "local", nlam = 2, delta = 1)
  expect_true(Matrix::isDiagonal(fit$Sig[[1]]))
})

test_that("solving pathwise gives same results as one-at-a-time", {
  del <- 1
  fit <- ggb::ggb(S, g, type = "local", delta = del, nlam = 10,
                  in_tol = 1e-5, out_tol = 1e-5)
  Sig2 <- list()
  for (l in seq_along(fit$lambda)) {
    Sig2[[l]] <- ggb::ggb(S, g, type = "local", lambda = fit$lambda[l],
                          delta = del, in_tol = 1e-5, out_tol = 1e-5)$Sig[[1]]
  }
  expect_equal(fit$Sig, Sig2, tolerance = 2e-4)
})
