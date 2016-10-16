context("disconnected seed graph")

n <- 20
p <- 10
set.seed(123)
x <- matrix(rnorm(n * p), n, p)
S <- cov(x)
adj <- matrix(0, p, p)
adj[1, 2] <- adj[2, 3] <- adj[3, 4] <- 1
adj[5, 6] <- adj[6, 7] <- 1
g <- igraph::graph.adjacency(adj, "undirected")
# igraph::components(g)$csize # components of size 4 3 1 1 1
tol <- 1e-4

test_that("output beats Convex.jl, global, no psd", {
  lamlist <- c(0.01, 0.1, 0.2)
  #cvx1 <- ggb_global_convexjl(S, lamlist[1], g = g, delta = NULL)
  #cvx2 <- ggb_global_convexjl(S, lamlist[2], g = g, delta = NULL)
  #cvx3 <- ggb_global_convexjl(S, lamlist[3], g = g, delta = NULL)
  fit <- ggb::ggb(S, g, type = "global", lambda = lamlist, delta = NULL)
  #expect_true(fit$obj[1] < cvx1$obj)
  expect_true(fit$obj[1] < 1.205731)
  #expect_true(fit$obj[2] < cvx2$obj)
  expect_true(fit$obj[2] < 1.348231)
  #expect_true(fit$obj[3] < cvx3$obj)
  expect_true(fit$obj[3] < 1.365835)
})

test_that("output beats Convex.jl, global, with psd", {
  lamlist <- c(0.01, 0.1, 0.2)
  del <- 0.5
  #cvx1 <- ggb_global_convexjl(S, lamlist[1], g = g, delta = del)
  #cvx2 <- ggb_global_convexjl(S, lamlist[2], g = g, delta = del)
  #cvx3 <- ggb_global_convexjl(S, lamlist[3], g = g, delta = del)
  fit <- ggb::ggb(S, g, type = "global", lambda = lamlist, delta = del)
  #expect_true(fit$obj[1] < cvx1$obj)
  expect_true(fit$obj[1] < 1.211355)
  expect_true(min(eigen(fit$Sig[[1]], only.values = TRUE)$values) > del - tol)
  #expect_true(fit$obj[2] < cvx2$obj)
  expect_true(fit$obj[2] < 1.34855)
  expect_true(min(eigen(fit$Sig[[2]], only.values = TRUE)$values) > del)
  #expect_true(fit$obj[3] < cvx3$obj)
  expect_true(fit$obj[3] < 1.36587)
  expect_true(min(eigen(fit$Sig[[3]], only.values = TRUE)$values) > del)
})

test_that("output beats Convex.jl, local, no psd", {
  lamlist <- c(0.01, 0.1, 0.2)
  #cvx1 <- ggb_local_convexjl(S, lamlist[1], g = g, delta = NULL)
  #cvx2 <- ggb_local_convexjl(S, lamlist[2], g = g, delta = NULL)
  #cvx3 <- ggb_local_convexjl(S, lamlist[3], g = g, delta = NULL)
  fit <- ggb::ggb(S, g, type = "local", lambda = lamlist, delta = NULL)
  #expect_true(fit$obj[1] < cvx1$obj)
  expect_true(fit$obj[1] < 1.204627)
  #expect_true(fit$obj[2] < cvx2$obj)
  expect_true(fit$obj[2] < 1.337951)
  #expect_true(fit$obj[3] < cvx3$obj)
  expect_true(fit$obj[3] < 1.365294)
})

test_that("output beats Convex.jl, global, with psd", {
  lamlist <- c(0.01, 0.1, 0.2)
  del <- 0.5
  #cvx1 <- ggb_local_convexjl(S, lamlist[1], g = g, delta = del)
  #cvx2 <- ggb_local_convexjl(S, lamlist[2], g = g, delta = del)
  #cvx3 <- ggb_local_convexjl(S, lamlist[3], g = g, delta = del)
  fit <- ggb::ggb(S, g, type = "local", lambda = lamlist, delta = del)
  #expect_true(fit$obj[1] < cvx1$obj + tol)
  expect_true(fit$obj[1] < 1.210306 + tol)
  expect_true(min(eigen(fit$Sig[[1]], only.values = TRUE)$values) > del - tol)
  #expect_true(fit$obj[2] < cvx2$obj)
  expect_true(fit$obj[2] < 1.338516 + tol)
  expect_true(min(eigen(fit$Sig[[2]], only.values = TRUE)$values) > del)
  #expect_true(fit$obj[3] < cvx3$obj)
  expect_true(fit$obj[3] < 1.365049 + 3 * tol)
  expect_true(min(eigen(fit$Sig[[3]], only.values = TRUE)$values) > del)
})

test_that("largest lambda gives diagonal matrix", {
  fit <- ggb::ggb(S, g, type = "local", nlam = 2, delta = 1)
  expect_true(Matrix::isDiagonal(fit$Sig[[1]]))

  fit <- ggb::ggb(S, g, type = "global", nlam = 2, delta = 1)
  expect_true(Matrix::isDiagonal(fit$Sig[[1]]))
})
