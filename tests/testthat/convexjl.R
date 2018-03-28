#julia <- "/Applications/Julia-0.6.app/Contents/Resources/julia/bin/julia"
julia <- "/Applications/Julia-0.5.app/Contents/Resources/julia/bin/julia"

ggb_local_convexjl <- function(S, lam, g, delta = NULL) {
  # Computes the local GGB estimator using Convex.jl
  #
  # Args:
  #  S:     sample covariance matrix
  #  lam:   positive tuning parameter
  #  g:     connected graph
  #  delta: (optional) nonnegative value giving minimum eigenvalue. If NULL,
  #         then no eigenvalue constraint included
  # this first line should be changed to location of julia on machine...
  p <- nrow(S)
  stopifnot(igraph::vcount(g) == p)
  D <- igraph::shortest.paths(g)
  D[D == Inf] <- 0
  depths <- apply(D, 1, max)
  w <- rep(NA, sum(depths))
  i <- 1
  for (j in seq(p)) {
    if (depths[j] == 0) next
    for (d in seq(depths[j])) {
      w[i] <- sqrt(2 * sum(D[j, ] <= d & D[j, ] > 0))
      i <- i + 1
    }
  }
  if (is.null(delta)) {
    psd <- ""; delta <- 0
  } else
    psd <- "pr.constraints += lambdamin(Sig) >= delta"
  pr.def <- paste(
    paste0("pr = minimize(0.5 * sumsquares(S - Sig) ",
           "+ lam * sum([w[i] * vecnorm(V[i],2) for i=1:sum(depths)]))"),
    "pr.constraints += (Sig == sum(V) + diagm(dd))",
    "i = 1",
    "for j = 1:p",
    "for d = 1:depths[j]",
    "if depths[j] == 0",
    "continue",
    "end",
    "notj = setdiff(1:p, j)",
    "pr.constraints += V[i][notj, notj] == 0",
    "njd = find(0 .< D[j, :] .<= d)", # d-neighborhood of j
    "njdc = setdiff(1:p, njd)",
    "pr.constraints += V[i][j, njdc] == 0",
    "pr.constraints += V[i][njdc, j] == 0",
    "i = i + 1",
    "end",
    "end",
    psd,
    sep = ";")
  # the following code changes the array of V to separate matrices V1, V2, ...
  code.after <- "function string_as_varname(s::String,v::Any)
  s=symbol(s)
  @eval (($s) = ($v))
  end;[string_as_varname(string('V',i),V[i]) for i = 1:length(V)]"
  opt.vars <- paste("V = Variable[Variable(p,p) for i=1:sum(depths)]",
                    "dd = Variable(p)",
                    "Sig = Variable(p, p)",
                    sep = ";")
  fit <- convexjulia::callconvex(opt.vars = opt.vars, pr.def = pr.def,
                                 code.after = code.after,
                                 const.vars = list(S = S, D = D,
                                                   lam = lam, p = p, w = w,
                                                   depths = depths,
                                                   delta = delta),
                                 opt.var.names = c(paste("V", 1:sum(depths),
                                                       sep = ""), "dd", "Sig"),
                                 julia.call = julia)
  V <- array(NA, c(p, p, sum(depths)))
  for (i in seq(sum(depths))) {
    V[, , i] <- fit[[paste0("V", i)]]
    fit[[paste0("V", i)]] <- NULL
  }
  fit$V <- V
  penval <- rep(0, p)
  norms <- sqrt(apply(V^2, 3, sum))
  i <- 1
  for (j in seq(p)) {
    if (depths[j] == 0) next
    for (d in seq(depths[j])) {
      penval[j] <- penval[j] + lam * w[i] * norms[i]
      i <- i + 1
    }
  }
  fit$w <- w
  fit$penval <- penval
  fit$obj <- 0.5 * sum((S - fit$Sig)^2) + sum(penval)
  fit
}

ggb_global_convexjl <- function(S, lam, g, delta = NULL) {
  # Computes the global GGB estimator using Convex.jl
  #
  # Args:
  #  S:     sample covariance matrix
  #  lam:   positive tuning parameter
  #  g:     connected graph
  #  delta: (optional) nonnegative value giving minimum eigenvalue. If NULL,
  #         then no eigenvalue constraint included
  # this first line should be changed to location of julia on machine...
  p <- nrow(S)
  stopifnot(igraph::vcount(g) == p)
  D <- igraph::shortest.paths(g)
  D[D == Inf] <- 0
  depth <- max(D)
  w <- rep(NA, depth)
  for (d in seq(depth)) {
    w[d] <- sqrt(sum(D <= d & D > 0))
  }
  if (is.null(delta)) {
    psd <- ""; delta <- 0
  } else
    psd <- "pr.constraints += lambdamin(Sig) >= delta"
  pr.def <- paste(
    paste0("pr = minimize(0.5 * sumsquares(S - Sig) ",
           "+ lam * sum([w[d] * vecnorm(V[d],2) for d=1:depth]))"),
    "pr.constraints += (Sig == sum(V) + diagm(dd))",
    "for d = 1:depth",
    "nd = find(0 .< D .<= d)", # d-neighborhood of j
    "ndc = setdiff(1:(p^2), nd)",
    "pr.constraints += V[d][ndc] == 0",
    "end",
    psd,
    sep = ";")
  # the following code changes the array of V to separate matrices V1, V2, ...
  code.after <- "function string_as_varname(s::String,v::Any)
  s=symbol(s)
  @eval (($s) = ($v))
  end;[string_as_varname(string('V',d),V[d]) for d = 1:length(V)]"
  opt.vars <- paste("V = Variable[Variable(p,p) for d=1:depth]",
                    "dd = Variable(p)",
                    "Sig = Variable(p, p)",
                    sep = ";")
  fit <- convexjulia::callconvex(opt.vars = opt.vars, pr.def = pr.def,
                                 code.after = code.after,
                                 const.vars = list(S = S, D = D,
                                                   lam = lam, p = p, w = w,
                                                   depth = depth,
                                                   delta = delta),
                                 opt.var.names = c(paste("V", 1:depth,
                                                         sep = ""), "dd", "Sig"),
                                 julia.call = julia)
  V <- array(NA, c(p, p, depth))
  for (d in seq(depth)) {
    V[, , d] <- fit[[paste0("V", d)]]
    fit[[paste0("V", d)]] <- NULL
  }
  fit$V <- V
  norms <- sqrt(apply(V^2, 3, sum))
  fit$w <- w
  fit$penval <- lam * sum(w * norms)
  fit$obj <- 0.5 * sum((S - fit$Sig)^2) + fit$penval
  fit
}
