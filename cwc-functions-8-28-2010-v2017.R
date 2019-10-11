## Functions for Clustering with Confidence (CWC) (5-28-2010)
## ==========================================================
## 8-25-2010: Added functions for assessing significance of splits

## We want to construct a cluster tree for which we are confident that
## the leaves actually correspond to modes of the density underlying
## the data and are not just to sampling variability.

## We have two ways of constructing a tree for a given confidence
## level 1-alpha: (1) based on significance tests and (2) based on
## CWC.

## (1) Pick one of excess mass, persistence, or size as a measure of
## prominence. Determine the distribution of the largest runt
## prominence value for uniform data sets of the same dimension and
## size as the sample and find the cutoff point for significance level
## alpha. Declare a node of the sample cluster tree as real if its
## runt prominence if larger than the cutoff.
## Why choose a uniform? How about choosing a unimodal density which
## is in some sense closest to the data? This is what one would do for
## 1-d data. How do we find such a unimodal?

## (2) Pick a type of confidence bounds (empirical, gaussian, or
## poisson).
## Empirical intervals are based on the quantiles of the resample
## weights. Gaussian intervals are of the form mu +- gamma_g * sd.
## Poisson intervals are of the form mu +- gamma _p * sqrt(mu).
## Find bounds for simultaneous coverage prob 1-alpha.  Do
## CWC

## The three variations of (2) may have very different real (as
## compared to nominal) signficance levels when applied to uniform
## data. For each them we can determine the corresponding nominal
## simultaneous coverage prob for a given real significance level.
## This calibration makes sense if we want to compare variations of
## (2) with (1). We can do power comparisons of the variations of (2)
## for the same simultaneous coverage prob.

## Code needed:
## For Gaussian and Poisson bounds find gamma_g and gamma_p for given
## simultaneous coverage prob. For empirical bounds find
## non-simultaneous coverage prob leading to desired simltaneous
## coverage prob.
## For poisson bounds find simultaneous coverage prob for given
## gamma_p.



###################################################################
## Chapter 1. Functions for computing resample vertex and edge weights
##
## Here is a function that takes as its arguments the vertices and
## edges of a test graph, and a density estimate. It computes the edge
## weights (minimum density along each of the edges) and the vertex
## weights (density at each vertex) using grid search.

cwc.ve.weights <- function(tg.vertices, tg.edges, density, ngrid) {
  nedges <- nrow(tg.edges)
  m <- ncol(tg.vertices)
  edge.weights <- rep(0, nedges)
  vertex.weights <- density(tg.vertices)
  t <- seq(0, 1, length = ngrid)[2:(ngrid-1)]
  Tmat <- matrix(t, nrow = ngrid-2, ncol = m)
  for (i in 1:nedges) {
    v1.index <- tg.edges[i, 1]
    v2.index <- tg.edges[i, 2]
    v1 <- tg.vertices[v1.index,]
    v2 <- tg.vertices[v2.index,]
    V1 <- matrix(v1, nrow = ngrid-2, ncol = m, byrow = T)
    V2 <- matrix(v2, nrow = ngrid-2, ncol = m, byrow = T)
    E <- (1-Tmat) * V1 + Tmat * V2
    phat <- c(vertex.weights[v1.index], density(E), vertex.weights[v2.index])
    edge.weights[i] <- min(phat)
  }
  return(list(vertex.weights = vertex.weights, edge.weights = edge.weights))
}

##-----------------------------------------------------------------
## Compute vertex and edge weights for original sample and resamples
## using kernel density estimate with LS CV.
##
## Version of 5-19-2010 below: Allow for bootstrap or half-sampling
## Cross validation fails for bootstrap resamples, presumably because
## of ties. Need to work out weighted version of LS CV??
## For the moment implement simple and not entirely kosher version:
## Use estimated bandwidth for original sample on all the bootstrap
## samples
##
## Choices for bandwidth:
## "cv.mean"   average CV bandwidth over resamples
## "cv"        individual CV bandwidth for each resamples
## a number    use this bandwidth for all resamples
##
## Revision 6-30-2010

cwc.resample.kernel.ve.weights <- function(X, tg.vertices = NULL, tg.edges = NULL,
                                           kmax = 20, include.mst = T,
                                           nres = 0, 
                                           resamples = NULL, trial.h = NULL,
                                           ngrid = 10, bandwidth = "cv",
                                           resample.type = "half.sampling",
                                           debug.filename = NULL) {
  debug <- !is.null(debug.filename)
  if (is.null(tg.vertices)){
    tg.edges <- cwc.determine.edges(X, kmax = kmax, include.mst = include.mst)
    tg.vertices <- X
  }
  if (is.null(trial.h))
    trial.h <- c(0.01, 0.05, 0.10, 0.15, 0.20, 0.30, 0.40, 0.60, 0.80,
                 1.00, 1.50, 2.00)
  nvertices <- nrow(tg.vertices)
  nedges <- nrow(tg.edges)
  n <- nrow(X)
  ##
  ## Compute vertex and edge weights for original sample
  ##
  if (debug) cat("\n\n\nComputing edge weights for original sample",
                 file = debug.filename, append = T)
  cvs.orig <- NULL
  if (is.numeric(bandwidth)) h.orig <- bandwidth
  if (!is.numeric(bandwidth)) {
    orig.cvs <- cv.search(X, trial.h)
    orig.h <- orig.cvs$opt.smopar
  }
  density <- make.gaussian.kernel.density.estimate(X, orig.h)
  cwcvw.out <- cwc.ve.weights(tg.vertices, tg.edges, density, ngrid)
  vertex.weights <- cwcvw.out$vertex.weights
  edge.weights <- cwcvw.out$edge.weights

  resample.vertex.weights <- NULL
  resample.edge.weights <- NULL
  resample.h <- NULL
  resample.cvs <- NULL

  if ((nres > 0) | (!is.null(resamples))) {
    if (is.null(resamples)) {
      if (resample.type == "half.sampling") rs.size <- round(nrow(X) / 2)
      if (resample.type == "bootstrap") rs.size <- nrow(X)
      resamples <- matrix(0, nrow = rs.size, ncol = nres)
      for (i in 1:nres) {
        if (resample.type == "half.sampling") rs <- sample(1:n, rs.size, replace = F)
        if (resample.type == "bootstrap") rs <- sample(1:n, rs.size, replace = T)
        resamples[,i] <- rs
      }
    }
    nres <- ncol(resamples)
    resample.type <- "half.sampling"
    if (nrow(resamples) == n) resample.type <- "bootstrap"
    resample.vertex.weights <- matrix(0, nrow = nvertices, ncol = nres)
    resample.edge.weights <- matrix(0, nrow = nedges, ncol = nres)
    resample.cvs <- vector("list", nres)
    resample.h <- rep(0, nres)

    if (debug) cat("\n\nFinding bandwidths for resamples",
                   file = debug.filename, append = T)
    if ((resample.type == "bootstrap") & !is.numeric(bandwidth)) {
      cvso <- cv.search(X, trial.h)
      resample.cvs[[1]] <- cvso
      resample.h[1:nres] <- cvso$opt.smopar
      if (debug) cat("\n", resample.h[1], file = debug.filename, append = T)
    } 
    if ((resample.type == "half.sampling") & !is.numeric(bandwidth)) {
      for (i in 1:nres) {
        rs <- resamples[,i]
        X.rs <- as.matrix(X[rs, ], nrow = rs.size)
        cvso <- cv.search(X.rs, trial.h)
        resample.cvs[[i]] <- cvso
        resample.h[i] <- cvso$opt.smopar
        if (debug) cat("\n", i, resample.h[i], file = debug.filename, append = T)
      }
      mean.h <- mean(resample.h[!is.na(resample.h)])
    }
    if (debug) cat("\n\nComputing resample vertex and edge weights\n",
                   file = debug.filename, append = T)
    for (i in 1:nres) {
      if (is.numeric(bandwidth)) h <- bandwidth
      if (!is.numeric(bandwidth)) {
        if (resample.type == "bootstrap") h <- resample.h[1]
        if (resample.type == "half.sampling") {
          h <- mean.h
          if ((bandwidth == "cv") & (!is.na(resample.h[i]))) h <- resample.h[i]
        }
      }
      rs <- resamples[,i]
      X.rs <- as.matrix(X[rs, ], nrow = rs.size)
      density <- make.gaussian.kernel.density.estimate(X.rs, h)
      cwcvw.out <- cwc.ve.weights(tg.vertices, tg.edges, density, ngrid)
      ##  browser()
      resample.vertex.weights[,i] <- cwcvw.out$vertex.weights
      resample.edge.weights[,i] <- cwcvw.out$edge.weights
      if(debug) cat(" ", i, file = debug.filename, append = T)
    }
  }
  return(list(tg.vertices = tg.vertices,
              tg.edges = tg.edges,
              vertex.weights = vertex.weights,
              edge.weights = edge.weights,
              orig.cvs = orig.cvs,
              resample.vertex.weights = resample.vertex.weights,
              resample.edge.weights = resample.edge.weights,
              resamples = resamples,
              resample.h = resample.h,
              resample.cvs = resample.cvs))
}

##=================================================================
## Now do same thing for nn density estimate

cwc.make.nn.density.estimate <- function(X.train) {
  density.estimate <- function(X.eval) {
     phat <- cwc.nn.density.estimate(X.train, X.eval)
     return(phat)
   }
  return(density.estimate)
}

cwc.nn.density.estimate <- function(X.obs, X.eval) {
  if (is.vector(X.obs)) {
    X.obs <- matrix(X.obs, ncol = 1)
    X.eval <- matrix(X.eval, ncol = 1)
  }
  m <- ncol(X.obs)
  nobs <- nrow(X.obs)
  neval <- nrow(X.eval)
  dens <- rep(0, neval)
  dist <- gsl.interpoint.distance.matrix(X.eval, X.obs)
  dist[dist < 0] <- 0
  dist <- sqrt(dist)
  for (i in 1:neval) {
    r <- min(dist[i,])
    ##if (r == 0) dens[i] <- my.Inf else dens[i] <- 1 / (n * r^m)
    if (r == 0) dens[i] <- my.Inf else dens[i] <- 1 / r
  }
  return(dens)
}

##-----------------------------------------------------------------

cwc.resample.nn.ve.weights <- function(X, tg.vertices = NULL, tg.edges = NULL,
                                       kmax = 20, include.mst = T,
                                       nres = 0, resamples = NULL, 
                                       ngrid = 10, resample.type = "half.sampling",
                                       debug.filename = NULL) {
  debug <- !is.null(debug.filename)
  if (is.null(tg.vertices)){
    tg.edges <- cwc.determine.edges(X, kmax = kmax, include.mst = include.mst)
    tg.vertices <- X
  }
  nvertices <- nrow(tg.vertices)
  nedges <- nrow(tg.edges)
  n <- nrow(X)
  ##
  ## Compute vertex and edge weights for original sample
  ##
  if (debug) cat("\n\n\nComputing edge weights for original sample",
                 file = debug.filename, append = T)
  density <- cwc.make.nn.density.estimate(X)
  cwcvw.out <- cwc.ve.weights(tg.vertices, tg.edges, density, ngrid)
  vertex.weights <- cwcvw.out$vertex.weights
  edge.weights <- cwcvw.out$edge.weights

  resample.vertex.weights <- NULL
  resample.edge.weights <- NULL

  if ((nres > 0) | (!is.null(resamples))) {
    if (is.null(resamples)) {
      if (resample.type == "half.sampling") rs.size <- round(nrow(X) / 2)
      if (resample.type == "bootstrap") rs.size <- nrow(X)
      resamples <- matrix(0, nrow = rs.size, ncol = nres)
      for (i in 1:nres) {
        if (resample.type == "half.sampling") rs <- sample(1:n, rs.size, replace = F)
        if (resample.type == "bootstrap") rs <- sample(1:n, rs.size, replace = T)
        resamples[,i] <- rs
      }
    }
    nres <- ncol(resamples)
    resample.type <- "half.sampling"
    if (nrow(resamples) == n) resample.type <- "bootstrap"
    resample.vertex.weights <- matrix(0, nrow = nvertices, ncol = nres)
    resample.edge.weights <- matrix(0, nrow = nedges, ncol = nres)
    if (debug) cat("\n\nComputing resample vertex and edge weights\n",
                   file = debug.filename, append = T)
    for (i in 1:nres) {
      rs <- resamples[,i]
      X.rs <- as.matrix(X[rs, ], nrow = rs.size)
      density <- cwc.make.nn.density.estimate(X.rs)
      cwcvw.out <- cwc.ve.weights(tg.vertices, tg.edges, density, ngrid)
      ##  browser()
      resample.vertex.weights[,i] <- cwcvw.out$vertex.weights
      resample.edge.weights[,i] <- cwcvw.out$edge.weights
      if(debug) cat(" ", i, file = debug.filename, append = T)
    }
  }
  return(list(tg.vertices = tg.vertices,
              tg.edges = tg.edges,
              vertex.weights = vertex.weights,
              edge.weights = edge.weights,
              resample.vertex.weights = resample.vertex.weights,
              resample.edge.weights = resample.edge.weights,
              resamples = resamples))
}


##-----------------------------------------------------------------

cwc.make.resample.mdsfun <- function(resample.weights) {
  tg.edges <- resample.weights$tg.edges
  rvw <- resample.weights$resample.vertex.weights
  rew <- resample.weights$resample.edge.weights
  nres <- ncol(rvw)
  nedges <- nrow(tg.edges)
  n <- nrow(rvw)
  resample.mdsfun <- function(i, j) {
    if (i == j) return(rvw[i,])
    for (k in 1:nedges) {
      if ((tg.edges[k, 1] == i) & (tg.edges[k, 2] == j))
        return(rew[k,])
      if ((tg.edges[k, 1] == j) & (tg.edges[k, 2] == i))
        return(rew[k,])
    }
    return(rep(0, nres))
  }
  attr(resample.mdsfun, "nres") <- nres
  attr(resample.mdsfun, "n") <- n
  return(resample.mdsfun)
}

##-----------------------------------------------------------------
## kmax = # of nearest neighbors, not including the point itself.

cwc.determine.edges <- function(X, kmax = 20, include.mst = T,
                                debug = F) {
  if (is.vector(X)) X <- matrix(X, ncol = 1)
  n <- nrow(X)
  p <- ncol(X)
  already.computed <- matrix(F, nrow = n, ncol = n)
  edges <- matrix(0, nrow = n * kmax + (n-1), ncol = 2)
  nedges <- 0
  Dist <- gsl.interpoint.distance.matrix(X, X)
  if (include.mst) {
    mst.mat.out <- gsl.mst.mat(Dist)
    mst.edges <- mst.mat.out$edges
    for (i in 1:(n-1)) {
      nedges <- nedges + 1
      edges[nedges,] <- mst.edges[i,]
      v1 <- mst.edges[i, 1]
      v2 <- mst.edges[i, 2]
      already.computed[v1, v2] <- T
      already.computed[v2, v1] <- T
    }
  }
  if(kmax >= 1) {
    for (i in 1:n) {
      v1 <- i
      NN <- order(Dist[i,])
      for (j in 2:(kmax + 1)) {
        v2 <- NN[j]
        if (!already.computed[v1, v2]) {
          nedges <- nedges + 1
          edges[nedges, 1] <- v1
          edges[nedges, 2] <- v2
          already.computed[v1, v2] <- T
          already.computed[v2, v1] <- T
        }
      }
    }
  }
  return(edges[1:nedges, ])
}

##################################################################
## Chapter 2: Original "Clustering with confidence" using simultaneous
## Confidence bands

## Current version of gsl.pruning.criteria can be used for CWC with
## arbitrary form of u(e) (upper confidence bound for edge weight) and
## l(v) (lower confidence bound for vertex weights).
## We choose a significance level ahead of time. The natural choice
## the pruning parameter is 0.
##
## We now add a version with poisson like confidence
## bounds for edge and vertex weights:
## Upper confidence bound u(e) for weight of edge e
## u(e) = wbar(e) + gamma * sqrt(wbar(e))
## Lower confidence bound l(v) for weight of vertex v
## l(v) = wbar(v) - gamma * sqrt(wbar(v))
##
## For this version we can use the resample means as the edge and
## vertex weights for  the test graph. Logically, we should use the
## upper confidence bounds, but they are monotonically increasing in
## the means, and therefore the mast will be the same.
##
## We can therefore compute, for each node of the mast dendogram, the
## the largest value of gamma for which the left resample rise would
## be bigger that 0, and ditto for the right resample rise. Call those
## left.gamma and right.gamma. Then we can compute
## runt.gamma = min(left.gamma, right.gamma).
##
## We can then assign a significance to each node which is the
## simultaneous coverage prob corresponding to runt.gamma.
##
## Ordinarily we prune for small values of the criterion (like runt
## size). Same here: large value of gamma means wide confidence band,
## which means large confidence level.
##

binary.root.search <- function(fun, xmin, xmax, maxit = 8) {
  nit = 0
  f.xmin <- fun(xmin)
  f.xmax <- fun(xmax)
  if (sign(f.xmin * f.xmax) > 0)
    stop("binary.root.search: initial bounds not a bracket")
  if (sign(f.xmin * f.xmax) == 0) {
    if (f.xmin == 0) return(xmin) else return(xmax)
  }
  for (i in 1:maxit) { 
    xmiddle <- (xmin + xmax) / 2
    f.xmiddle <- fun(xmiddle)
    if (f.xmiddle == 0) return(xmiddle)
    if (sign(f.xmiddle * f.xmin) < 0) {
      xmax <- xmiddle
      f.xmax <- f.xmiddle
      next
    }
    if (sign(f.xmiddle * f.xmax) < 0) {
      xmin <- xmiddle
      f.xmin <- f.xmiddle
      next
    }
  }
  return(xmiddle)
}

## test.fun <- function(x) return(round(10 * sin(x)))
## binary.root.search(test.fun, -3, 3)
  
##-----------------------------------------------------------------    
## For Gaussian and poisson interval, cov.parameter = gamma.
## For empirical intervals, cov.parameter = nout

cwc.simul.cov.prob <- function(resample.weights, cov.parameter,
                               band.type = "empirical") {
  rvw <- resample.weights$resample.vertex.weights
  rew <- resample.weights$resample.edge.weights
  nres <- ncol(rew)
  ne <- nrow(rew)
  nv <- nrow(rvw)

  if (band.type == "empirical") {
    nout <- cov.parameter
    lower <- nout + 1
    upper <- nres - nout
    if (lower >= upper) {
      lower <- lower - 1
      upper <- upper + 1
    }
    upper.ew <- rep(0, ne)
    upper.vw <- rep(0, nv)
    lower.vw <- rep(0, nv)
    for (i in 1:ne) {
      sew <- sort(rew[i,])
      upper.ew[i] <- sew[upper]
    }
    for (i in 1:nv) {
      svw <- sort(rvw[i,])
      upper.vw[i] <- svw[upper]
      lower.vw[i] <- svw[lower]
    }
  }
  if (band.type != "empirical") {
    gamma <- cov.parameter
    ew.mean <- apply(rew, 1, mean)
    ew.sd <- apply(rew, 1, sd)
    vw.mean <- apply(rvw, 1, mean)
    vw.sd <- apply(rvw, 1, sd)
  }
  if (band.type == "gaussian") {
    upper.ew <- ew.mean + gamma * ew.sd
    upper.vw <- vw.mean + gamma * vw.sd
    lower.vw <- vw.mean - gamma * vw.sd
  }
  if (band.type == "poisson") {
    vw.sd <- sqrt(vw.mean)
    upper.ew <- ew.mean + gamma * sqrt(ew.mean)
    upper.vw <- vw.mean + gamma * sqrt(vw.mean)
    lower.vw <- vw.mean - gamma * sqrt(vw.mean)
  }
  nout <- 0
  out.resamples <- rep(0, nres)
  for (i in 1:nres) {
    out <- sum(rew[,i] > upper.ew) + sum(rvw[,i] > upper.vw) +
      sum(rvw[,i] < lower.vw)
    if (out > 0) {
      nout <- nout + 1
      out.resamples[nout] <- i
    }
  }
  return(list(simul.cov.prob = 1 - nout/nres, out.resamples =
    out.resamples[1:nout]))

}


##-----------------------------------------------------------------

cwc.simul.confidence.bounds <- function(resample.weights, target.simul.cov.prob,
                                  band.type = "empirical") {
  make.fun <- function(rvw, rew, band.type, target.sim.cov.prob) {
    fun <- function(cov.parameter) {
      return(cwc.simul.cov.prob(resample.weights, cov.parameter,
                                band.type)$simul.cov.prob
             - target.sim.cov.prob)
    }
    return(fun)
  }
  rvw <- resample.weights$resample.vertex.weights
  rew <- resample.weights$resample.edge.weights
  ne <- nrow(rew)
  nv <- nrow(rvw)
  nres <- ncol(rvw)

  if ((band.type == "gaussian") | (band.type == "poisson")) {
    xmin <- 0
    ## xmax <- 4
    xmax <- 10  ## changed 11-02-2018
    fun <- make.fun(rvw, rew,  band.type, target.simul.cov.prob)
    gamma <- binary.root.search(fun, xmin, xmax, maxit = 12)
    ew.mean <- apply(rew, 1, mean)
    ew.sd <- apply(rew, 1, sd)
    vw.mean <- apply(rvw, 1, mean)
    vw.sd <- apply(rvw, 1, sd)
  }
  if (band.type == "gaussian") {
    upper.ew <- ew.mean + gamma * ew.sd
    upper.vw <- vw.mean + gamma * vw.sd
    lower.vw <- vw.mean - gamma * vw.sd
  }
  if (band.type == "poisson") {
    vw.sd <- sqrt(vw.mean)
    upper.ew <- ew.mean + gamma * sqrt(ew.mean)
    upper.vw <- vw.mean + gamma * sqrt(vw.mean)
    lower.vw <- vw.mean - gamma * sqrt(vw.mean)
  }
  if (band.type == "empirical") {
    xmin <- 0
    xmax <- nres / 2
    fun <- make.fun(rvw, rew, band.type, target.simul.cov.prob)
    nout <- round(binary.root.search(fun, xmin, xmax, maxit = 8))
    gamma <- nout
    lower <- nout + 1
    upper <- nres - nout
    if (lower >= upper) {
      lower <- lower - 1
      upper <- upper + 1
    }
    upper.ew <- rep(0, ne)
    upper.vw <- rep(0, nv)
    lower.vw <- rep(0, nv)
    for (i in 1:ne) {
      sew <- sort(rew[i,])
      upper.ew[i] <- sew[upper]
    }
    for (i in 1:nv) {
      svw <- sort(rvw[i,])
      upper.vw[i] <- svw[upper]
      lower.vw[i] <- svw[lower]
    }
  }
  return(list(upper.ew = upper.ew, upper.vw = upper.vw, lower.vw =
         lower.vw, cov.parameter = gamma))
}

##-----------------------------------------------------------------

cwc.nonsimul.confidence.bounds <- function(resample.weights, target.cov.prob,
                                           band.type = "empirical") {
  rvw <- resample.weights$resample.vertex.weights
  rew <- resample.weights$resample.edge.weights
  ne <- nrow(rew)
  nv <- nrow(rvw)
  nres <- ncol(rvw)

  if (band.type == "empirical") {
    nout <- ceiling(nres * (1-target.cov.prob)/2)
    lower <- nout + 1
    upper <- nres - nout
    if (lower >= upper) {
      lower <- lower - 1
      upper <- upper + 1
    }
    upper.ew <- rep(0, ne)
    lower.vw <- rep(0, nv)
    upper.vw <- rep(0, nv)
    for (i in 1:ne) {
      sew <- sort(rew[i,])
      upper.ew[i] <- sew[upper]
    }
    for (i in 1:nv) {
      svw <- sort(rvw[i,])
      lower.vw[i] <- svw[lower]
      upper.vw[i] <- svw[upper]
    }
    cov.parameter <- nout
  }
  if ((band.type == "gaussian") | (band.type == "poisson")) {
    z.cov <- qnorm(1 - (1-target.cov.prob)/2)
    cov.parameter <- z.cov
    ew.mean <- apply(rew, 1, mean)
    ew.sd <- apply(rew, 1, sd)
    vw.mean <- apply(rvw, 1, mean)
    vw.sd <- apply(rvw, 1, sd)
    if(band.type == "gaussian") {
      upper.ew <- ew.mean + z.cov * ew.sd
      lower.vw <- vw.mean - z.cov * vw.sd
      upper.vw <- vw.mean + z.cov * vw.sd
    }
    if(band.type == "poisson") {
      upper.ew <- ew.mean + z.cov * sqrt(ew.mean)
      lower.vw <- vw.mean - z.cov * sqrt(vw.mean)
      upper.vw <- vw.mean + z.cov * sqrt(vw.mean)
    }
  }
  return(list(upper.ew = upper.ew, upper.vw = upper.vw, lower.vw =
         lower.vw, cov.parameter = cov.parameter))
}

##-----------------------------------------------------------------

cwc.confidence.bounds <- function(resample.weights, target.cov.prob,
                                  band.type = "empirical", simultaneous = F){
  if(simultaneous)
    confidence.bounds <- cwc.simul.confidence.bounds(resample.weights, target.cov.prob,
                                  band.type = band.type)
  if(!simultaneous)
    confidence.bounds <- cwc.nonsimul.confidence.bounds(resample.weights, target.cov.prob,
                                  band.type = band.type)
  return(confidence.bounds)
}

##-----------------------------------------------------------------

cwc.nonsimul.upper.poisson.bounds <- function(rew, rvw, alpha) {
  xmin <- 0
  xmax <- 4
  maxit <- 12
  mean.rew <- apply(rew, 1, mean)
  mean.rvw <- apply(rvw, 1, mean)
  nedges <- nrow(rew)
  n <- nrow(rvw)
  nres <- ncol(rew)
  fun <- function(cov.parameter) {
    max.out <- 0
    for (i in 1:nedges) {
      upper <- mean.rew[i] + cov.parameter * sqrt(mean.rew[i])
      out <- sum(rew[i,] > upper)
      if (out > max.out) max.out <- out
    }
    for (i in 1:n) {
      upper <- mean.rvw[i] + cov.parameter * sqrt(mean.rvw[i])
      out <- sum(rvw[i,] > upper)
      if (out > max.out) max.out <- out
    }
    return(max.out/nres - alpha)
  }
  gamma <- binary.root.search(fun, xmin, xmax, maxit = maxit)
  upper.ew <- mean.rew + gamma * sqrt(mean.rew)
  upper.vw <- mean.rvw + gamma * sqrt(mean.rvw)
  return(list(upper.vw = upper.vw, upper.ew = upper.ew))
}

##-----------------------------------------------------------------

cwc.make.upper.confidence.mdsfun <- function(resample.weights, confidence.bounds) {
  tg.vertices <- resample.weights$tg.vertices
  tg.edges <- resample.weights$tg.edges
  n <- nrow(tg.vertices)
  nedges <- nrow(tg.edges)
  uc.mdsmat <- matrix(0, nrow = n, ncol = n)
  diag(uc.mdsmat) <- confidence.bounds$upper.vw
  for (i in 1:nedges) {
    v1 <- tg.edges[i, 1]
    v2 <- tg.edges[i, 2]
    uc.mdsmat[v1, v2] <- confidence.bounds$upper.ew[i]
    uc.mdsmat[v2, v1] <- confidence.bounds$upper.ew[i]
  }
  up <- as.vector(uc.mdsmat[upper.tri(uc.mdsmat, diag = T)])
  upper.confidence.mdsfun <- gsl.make.mdsfun.from.mdsmat(up)
  return(upper.confidence.mdsfun)
}

##-----------------------------------------------------------------

cwc.make.resample.mean.mdsfun <- function(resample.weights) {
  tg.vertices <- resample.weights$tg.vertices
  tg.edges <- resample.weights$tg.edges
  n <- nrow(tg.vertices)
  nedges <- nrow(tg.edges)
  uc.mdsmat <- matrix(0, nrow = n, ncol = n)
  resample.edge.weights <- resample.weights$resample.edge.weights
  resample.vertex.weights <- resample.weights$resample.vertex.weights
  mean.edge.weights <- apply(resample.edge.weights, 1, mean)
  mean.vertex.weights <- apply(resample.vertex.weights, 1, mean)
  diag(uc.mdsmat) <- mean.vertex.weights
  for (i in 1:nedges) {
    v1 <- tg.edges[i, 1]
    v2 <- tg.edges[i, 2]
    uc.mdsmat[v1, v2] <- mean.edge.weights[i]
    uc.mdsmat[v2, v1] <- mean.edge.weights[i]
  }
  up <- as.vector(uc.mdsmat[upper.tri(uc.mdsmat, diag = T)])
  upper.confidence.mdsfun <- gsl.make.mdsfun.from.mdsmat(up)
  return(upper.confidence.mdsfun)
}

##-----------------------------------------------------------------

cwc.make.mdsfun.from.weights <- function(tg.edges, vertex.weights,
                                         edge.weights) {
  n <- length(vertex.weights)
  nedges <- nrow(tg.edges)
  uc.mdsmat <- matrix(0, nrow = n, ncol = n)
  diag(uc.mdsmat) <- vertex.weights
  for (i in 1:nedges) {
    v1 <- tg.edges[i, 1]
    v2 <- tg.edges[i, 2]
    uc.mdsmat[v1, v2] <- edge.weights[i]
    uc.mdsmat[v2, v1] <- edge.weights[i]
  }
  up <- as.vector(uc.mdsmat[upper.tri(uc.mdsmat, diag = T)])
  mdsfun <- gsl.make.mdsfun.from.mdsmat(up)
  return(mdsfun)
}



##=================================================================
## Utilities 

  
dump.to.string <- function(object.names) {
  tc <- textConnection("dump.out", "w")
  sink(tc)
  dump(object.names, file = "")
  sink()
  close(tc)
  return(dump.out)
}

source.from.string <- function(string) {
  tc <- textConnection(string, "r")
  source(tc)
  close(tc)
}

## resample.weights.filename <- function(data.name, resample.type,
##                                       nres, bandwidth,
##                                       kmax, ngrid) {
##   rt.abb <- "hs"
##   if (resample.type == "bootstrap") rt.abb <- "bt"
##   bw.abb <- "cvm"
##   if (bandwidth == "cv") bw.abb <- "cv"
##   if (bandwidth == "nn") bw.abb <- "nn"
##   if (is.numeric(bandwidth)) bw.abb <- paste("n", bandwidth, sep = "")
##   resample.weights.name <- paste(data.name, "rw", rt.abb, nres,
##                                  bw.abb, kmax, ngrid, sep = "-")
##   resample.weights.filename <- paste(resample.weights.name, ".R",
##                                      sep = "")
##   return(resample.weights.filename)
## }

##-----------------------------------------------------------------
## Version below does not add extension
resample.weights.filename <- function(data.name, resample.type,
                                      nres, bandwidth,
                                      kmax, ngrid) {
  rt.abb <- "hs"
  if (resample.type == "bootstrap") rt.abb <- "bt"
  bw.abb <- "cvm"
  if (bandwidth == "cv") bw.abb <- "cv"
  if (bandwidth == "nn") bw.abb <- "nn"
  if (is.numeric(bandwidth)) bw.abb <- paste("n", bandwidth, sep = "")
  resample.weights.filename <- paste(data.name, "rw", rt.abb, nres,
                                 bw.abb, kmax, ngrid, sep = "-")
  return(resample.weights.filename)
}


#################################################################
## Chapter 3: Assessing confidence in splits.
## The idea of using simultaneous confidence bands seems to fail for
## high-dimensional data where bandwidths are small.
##
## Alternative idea: Assess confidence in individual splits In a non-
## simultaneous way.

gsl.draw.tree.pick.node <- function(gsl.cluster.out, tree.type = "cluster.tree", ...) {
  cluster.tree <- gsl.cluster.out$cluster.tree
  obs.leaf.code <- gsl.cluster.out$leaf.code
  obs.cluster.core <- gsl.cluster.out$cluster.core
  if (tree.type == "cluster.tree") {
    node.xy <- gsl.draw.cluster.tree(gsl.cluster.out, ...)
    title("Cluster tree")
  }
  else {
    node.xy <- gsl.draw.cluster.dendogram(gsl.cluster.out, ...)
    title("Cluster dendogram")
  }
  selected.node <- identify(node.xy$x, node.xy$y, labels = rep(" ", length(node.xy$x)),
                            n = 1)
  points(node.xy$x[selected.node], node.xy$y[selected.node], pch = 18, cex = 2, 
         col = "blue1")
  return(selected.node)
}


##-----------------------------------------------------------------
## There are many ways of deciding which resample samples support
## a split. Here is one:
## 
## A resample sample supports the split if its edge weight for the
## split edge is <= some threshold opt.level and there are
## observations in the left and right cluster cores for which the
## resample density is > opt.level. Here opt.level is the same for
## all resample samples
##


cwc.find.bs.samples.supporting.split <- function(broken.edge.weights,
                                                 min.core.maxes,
                                                 uniform.level = T) {
  nres <- length(broken.edge.weights)
  if (uniform.level) {
    obew <- order(broken.edge.weights)
    sbew <- sort(broken.edge.weights)
    smcw <- min.core.maxes[obew]
    supporting.bs <- NULL
    for (i in 1:nres) {
      n.supporting.bs <- sum(smcw[1:i] > sbew[i])
      if (n.supporting.bs > length(supporting.bs)) {
        opt.level <- sbew[i]
        supporting.bs <- (1:nres)[(broken.edge.weights <= opt.level) &
                                   (min.core.maxes > opt.level)]
       }
    }
  }
  else {
    opt.level <- NULL
    supporting.bs <- (1:nres)[min.core.maxes > broken.edge.weights]
  }
  return(list(opt.level = opt.level, supporting.bs = supporting.bs))
}

##=================================================================
## Here is another version:
##
## Local version of CWC: Grow tree using non-simultaneous upper
## 1-alpha/2 confidence bounds for edge weights. Get left and right
## cluster cores C_l, C_r.
## Option 1: Let e be the split edge and w(e) the split edge weight.
## Let S be the set of resamples for which w(e,i) <= w(e).
## For each resample sample i in S and each point in C_l or C_r
## determine whether w(v, i) > w(e). Question: Are there points
## j_l, j_r such that for 1-alpha/2 of the resample samples in S
## w(j_l, i) > w(e) and w(j_r, i) > w(e). If such points exist,
## declare the split significant at level alpha
## 

cwc.assess.split.significance <- function(gsl.cluster.out, inode,
                                          resample.mdsfun) {
  ##browser()
  nres <- attr(resample.mdsfun, "nres")
  n <- attr(resample.mdsfun, "n")                                    
  cluster.tree <- gsl.cluster.out$cluster.tree
  mast.dendogram <- gsl.cluster.out$mast.dendogram
  obs.leaf.code <- gsl.cluster.out$leaf.code
  ## nobs <- length(obs.leaf.code)
  obs.in.core <- gsl.cluster.out$cluster.core
  leftson <- cluster.tree[inode, "leftson"]
  if (leftson == 0)
    stop("assess.split.significance: you should select an internal node")
  rightson <- cluster.tree[inode, "rightson"]
  split.edge.index <- cluster.tree[inode, "dendogram.node"]
  split.edge <- mast.dendogram$edges[split.edge.index,]
  split.edge.weight <- mast.dendogram$height[split.edge.index]
  left.leaf.code <- cluster.tree[leftson, "leaf.code"]
  right.leaf.code <- cluster.tree[rightson, "leaf.code"]
  split.edge.bt.weights <- resample.mdsfun(split.edge[1], split.edge[2])
  in.left.core <- gsl.select.obs.in.node(obs.leaf.code, obs.in.core,
                                         left.leaf.code)$in.core
  in.right.core <- gsl.select.obs.in.node(obs.leaf.code, obs.in.core,
                                          right.leaf.code)$in.core
  obs.in.left.core <- (1:n)[in.left.core]
  obs.in.right.core <- (1:n)[in.right.core]
  nleft <- length(obs.in.left.core)
  nright <- length(obs.in.right.core)
  left.core.bt.weights <- matrix(0, nrow = nleft, ncol = nres)
  right.core.bt.weights <- matrix(0, nrow = nright, ncol = nres)
  for (i in 1:nleft) {
    obs.ind <- obs.in.left.core[i]
    left.core.bt.weights[i, ] <- resample.mdsfun(obs.ind, obs.ind)
  }
  for (i in 1:nright) {
    obs.ind <- obs.in.right.core[i]
    right.core.bt.weights[i, ] <- resample.mdsfun(obs.ind, obs.ind)
  }
  split.edge.indicators <- split.edge.bt.weights <= split.edge.weight
  left.core.indicators <- left.core.bt.weights > split.edge.weight
  right.core.indicators <- right.core.bt.weights > split.edge.weight
  max.simul.count <- 0
  best.left.ind <- 0
  best.right.ind <- 0
  for (i in 1:nleft) {
    for (j in 1:nright) {
      simul.count <- sum(left.core.indicators[i,] *
                         right.core.indicators[j,]*
                         split.edge.indicators)
      if (simul.count > max.simul.count) {
        best.left.ind <- i
        best.right.ind <- j
        max.simul.count <- simul.count
      }
    }
  }
  best.left.obs <- obs.in.left.core[best.left.ind]
  best.right.obs <- obs.in.right.core[best.right.ind]
  supporting.resamples <- (1:nres)[as.logical(split.edge.indicators *
                                    left.core.indicators[best.left.ind,] *
                                    right.core.indicators[best.right.ind,])]
  return(list(split.edge = split.edge,
              best.left.obs = best.left.obs,
              best.right.obs = best.right.obs,
              supporting.resamples = supporting.resamples))
}

##-----------------------------------------------------------------

cwc.assess.all.nodes <- function(gsl.cluster.out, resample.mdsfun) {
  cluster.tree <- gsl.cluster.out$cluster.tree
  nnodes <- nrow(cluster.tree)
  nres <- attr(resample.mdsfun, "nres")
  node.ass.out <- vector("list", nnodes)
  node.confidence.level <- rep(0,nnodes)
  for (i in 1:nnodes){
    if (cluster.tree[i, "leftson"] == 0) next
    ass.out <- cwc.assess.split.significance(gsl.cluster.out, i,
                                             resample.mdsfun)
    n.supporting.resamples <- length(ass.out$supporting.resamples)
    node.confidence.level[i] <- n.supporting.resamples / nres
    node.ass.out[[i]] <- ass.out
  }
  return(list(node.confidence.level = node.confidence.level,
              node.ass.out = node.ass.out))
}

##=================================================================
## And yet another one
##
## Grow tree using original or bagged edge weights.
## For chosen edge mast edge e get left and right cluster cores C_l, C_r.
## Let v_l, v_r be the points in C_l, C_r with maximum vertex weights.
## Find minimum alpha_l for which upper (1-alpha) cb for edge weights
## of e is <= lower that lower (1-alpha) cb for vertex weights of v_l
## Find alpha_r the same way. Set alpha = max(alpha_l, alpha_r).

cwc.assess.split.significance.v2 <- function(gsl.cluster.out, inode,
                                             resample.mdsfun) {
  ##browser()
  nres <- attr(resample.mdsfun, "nres")
  n <- attr(resample.mdsfun, "n")                                    
  cluster.tree <- gsl.cluster.out$cluster.tree
  mast.dendogram <- gsl.cluster.out$mast.dendogram
  obs.leaf.code <- gsl.cluster.out$leaf.code
  obs.in.core <- gsl.cluster.out$cluster.core
  obs.density <- gsl.cluster.out$obs.density
  ##
  leftson <- cluster.tree[inode, "leftson"]
  if (leftson == 0)
    stop("assess.split.significance: you should select an internal node")
  rightson <- cluster.tree[inode, "rightson"]
  split.edge.index <- cluster.tree[inode, "dendogram.node"]
  split.edge <- mast.dendogram$edges[split.edge.index,]
  left.leaf.code <- cluster.tree[leftson, "leaf.code"]
  right.leaf.code <- cluster.tree[rightson, "leaf.code"]
  in.left.core <- gsl.select.obs.in.node(obs.leaf.code, obs.in.core,
                                         left.leaf.code)$in.core
  in.right.core <- gsl.select.obs.in.node(obs.leaf.code, obs.in.core,
                                          right.leaf.code)$in.core
  obs.in.left.core <- (1:n)[in.left.core]
  obs.in.right.core <- (1:n)[in.right.core]
  left.densities <- obs.density[obs.in.left.core]
  right.densities <- obs.density[obs.in.right.core]
  best.left.obs <- obs.in.left.core[which.max(left.densities)]
  best.right.obs <- obs.in.right.core[which.max(right.densities)]
  left.max.resample.weights <- resample.mdsfun(best.left.obs, best.left.obs)
  right.max.resample.weights <- resample.mdsfun(best.right.obs, best.right.obs)
  split.edge.resample.weights <- resample.mdsfun(split.edge[1], split.edge[2])
  ##
  rev.sorted.serw <- rev(sort(split.edge.resample.weights))
  sorted.lmrw <- sort(left.max.resample.weights)
  sorted.rmrw <- sort(right.max.resample.weights)
  for (il in 1:nres) if(rev.sorted.serw[il] <= sorted.lmrw[il]) break
  for (ir in 1:nres) if(rev.sorted.serw[il] <= sorted.lmrw[ir]) break
  sig <- 2 * max(il, ir) / nres
  confidence.level <- 1 - sig
  
  return(list(split.edge = split.edge,
              best.left.obs = best.left.obs,
              best.right.obs = best.right.obs,
              confidence.level = confidence.level))
}

##-----------------------------------------------------------------

cwc.assess.all.nodes.v2 <- function(gsl.cluster.out, resample.mdsfun) {
  cluster.tree <- gsl.cluster.out$cluster.tree
  nnodes <- nrow(cluster.tree)
  nres <- attr(resample.mdsfun, "nres")
  node.ass.out <- vector("list", nnodes)
  node.confidence.level <- rep(0,nnodes)
  for (i in 1:nnodes){
    if (cluster.tree[i, "leftson"] == 0) next
    ass.out <- cwc.assess.split.significance.v2(gsl.cluster.out, i,
                                                resample.mdsfun)
    node.confidence.level[i] <- ass.out$confidence.level
    node.ass.out[[i]] <- ass.out
  }
  return(list(node.confidence.level = node.confidence.level,
              node.ass.out = node.ass.out))
}

##-----------------------------------------------------------------

##=================================================================
## Here is version 3
##
## Grow tree using original or bagged edge weights.
## For chosen edge mast edge e get left and right cluster cores C_l, C_r.
## Let v_l, v_r be the points in C_l, C_r with maximum vertex weights.
## Compute resampling distribution for
##   T = log(min(phat(v_l), phat(v_r)) - log(phat(e))
## and P(T < 0)


cwc.assess.split.significance.v3 <- function(gsl.cluster.out, inode,
                                             resample.mdsfun) {
  ##browser()
  nres <- attr(resample.mdsfun, "nres")
  n <- attr(resample.mdsfun, "n")                                    
  cluster.tree <- gsl.cluster.out$cluster.tree
  mast.dendogram <- gsl.cluster.out$mast.dendogram
  obs.leaf.code <- gsl.cluster.out$leaf.code
  obs.in.core <- gsl.cluster.out$cluster.core
  obs.density <- gsl.cluster.out$obs.density
  log.dip <- rep(NA, nres)
  ##
  leftson <- cluster.tree[inode, "leftson"]
  if (leftson == 0)
    stop("assess.split.significance: you should select an internal node")
  rightson <- cluster.tree[inode, "rightson"]
  split.edge.index <- cluster.tree[inode, "dendogram.node"]
  split.edge <- mast.dendogram$edges[split.edge.index,]
  left.leaf.code <- cluster.tree[leftson, "leaf.code"]
  right.leaf.code <- cluster.tree[rightson, "leaf.code"]
  in.left.core <- gsl.select.obs.in.node(obs.leaf.code, obs.in.core,
                                         left.leaf.code)$in.core
  in.right.core <- gsl.select.obs.in.node(obs.leaf.code, obs.in.core,
                                          right.leaf.code)$in.core
  obs.in.left.core <- (1:n)[in.left.core]
  obs.in.right.core <- (1:n)[in.right.core]
  left.densities <- obs.density[obs.in.left.core]
  right.densities <- obs.density[obs.in.right.core]
  best.left.obs <- obs.in.left.core[which.max(left.densities)]
  best.right.obs <- obs.in.right.core[which.max(right.densities)]
  left.max.resample.weights <- resample.mdsfun(best.left.obs, best.left.obs)
  right.max.resample.weights <- resample.mdsfun(best.right.obs, best.right.obs)
  split.edge.resample.weights <- resample.mdsfun(split.edge[1], split.edge[2])
  ##
  for (i in 1:nres) {
    log.dip[i] <- log(min(left.max.resample.weights[i],
                          right.max.resample.weights[i])) -
                  log(split.edge.resample.weights[i])
  }
  significance <- sum(log.dip < 0) / nres
  return(list(split.edge = split.edge,
              best.left.obs = best.left.obs,
              best.right.obs = best.right.obs,
              log.dip = log.dip,
              significance = significance))
}

##-----------------------------------------------------------------

cwc.assess.all.nodes.v3 <- function(gsl.cluster.out, resample.mdsfun) {
  cluster.tree <- gsl.cluster.out$cluster.tree
  nnodes <- nrow(cluster.tree)
  nres <- attr(resample.mdsfun, "nres")
  node.ass.out <- vector("list", nnodes)
  node.significance <- rep(NA, nnodes)
  for (i in 1:nnodes){
    if (cluster.tree[i, "leftson"] == 0) next
    ass.out <- cwc.assess.split.significance.v3(gsl.cluster.out, i,
                                                resample.mdsfun)
    node.significance[i] <- ass.out$significance
    node.ass.out[[i]] <- ass.out
  }
  return(list(node.significance = node.significance,
              node.ass.out = node.ass.out))
}

##-----------------------------------------------------------------

cwc.make.mdsfun.for.tree <- function(mdsfun, tree.edges) {
  n <- attr(mdsfun, "n")
  mdsmat <- matrix(0, nrow = n, ncol = n)
  for (i in 1:n) mdsmat[i, i] <- mdsfun(i, i)
  for (i in 1:(n-1)) {
    v1 <- tree.edges[i, 1]
    v2 <- tree.edges[i, 2]
    mdsmat[v1, v2] <- mdsfun(v1, v2)
    mdsmat[v2, v1] <- mdsmat[v1, v2]
  }
  up <- as.vector(mdsmat[upper.tri(mdsmat, diag = T)])
  mdsfun <- gsl.make.mdsfun.from.mdsmat(up)
  return(mdsfun)
}

##-----------------------------------------------------------------

cwc.runt.statistics.for.edges <- function(gsl.cluster.out) {
  mast.dendogram.edges <- (gsl.cluster.out$mast.dendogram)$edges
  pc.values <- gsl.cluster.out$pc.values
  nedges <- nrow(mast.dendogram.edges)
  rse <- matrix(0, nrow = nedges, ncol = 5)
  colnames(rse) <- c("runt.size", "runt.excess.mass", "runt.rise",
                     "runt.bootstrap.rise", "runt.log.ratio" )
  rse[,"runt.size"] <- pmin(pc.values[,"left.size"], pc.values[,"right.size"])
  rse[,"runt.excess.mass"] <- pmin(pc.values[,"left.excess.mass"],
                                   pc.values[,"right.excess.mass"])
  rse[,"runt.rise"] <- pmin(pc.values[,"left.rise"], pc.values[,"right.rise"])
  rse[,"runt.bootstrap.rise"] <- pmin(pc.values[,"left.bootstrap.rise"],
                                      pc.values[,"right.bootstrap.rise"])
  rse[,"runt.log.ratio"] <- pmin(pc.values[,"left.log.ratio"],
                                 pc.values[,"right.log.ratio"])
  return(list(mast.edges = mast.dendogram.edges, runt.statistics = rse))
}



#################################################################
## Chapter 4: Prediction strength
##
## Function produces a matrix of prediction strengths with maxclust
## rows and 2 columns. The column pred.strength[, "test"] has the
## prediction strengths for the original assignment of the attributes
## "train" and "test". The column pred.strength[, "train"] has the
## prediction strengths for the reverse assignment.

prediction.strength <- function(X.train, mdsfun.train, X.test, mdsfun.test,
                                pruning.criterion, maxclust) {
  test.plot <- function() {
    colors <- c("black", "red", "blue", "green", "yellow", "orange") ###
    mfrow.orig <- par("mfrow")
    par(mfrow = c(2,2))
    plot(X1, type = "n")
    for (j in 1:i) points(X1[lab1 == j,], col = colors[j], pch = 20)
    title("lab1")
    plot(X1, type = "n")
    for (j in 1:i) points(X1[lab1.pred == j,], col = colors[j], pch = 20)
    title("lab1.pred")
    ##
    plot(X2, type = "n")
    for (j in 1:i) points(X2[lab2 == j,], col = colors[j], pch = 20)
    title("lab2")
    plot(X2, type = "n")
    for (j in 1:i) points(X2[lab2.pred == j,], col = colors[j], pch = 20)
    title("lab2.pred")
    par(mfrow = mfrow.orig)
  }
  ##
  ntrain <- nrow(X.train)
  ntest <- nrow(X.test)
  cluster.train <- gsl.cluster(X.train, mdsfun.train, assign.fluff = F,
                               pruning.criterion = pruning.criterion)
  cluster.test <- gsl.cluster(X.test, mdsfun.test, assign.fluff = F,
                              pruning.criterion = pruning.criterion)
  pred.strength <- matrix(1, nrow = maxclust, ncol = 2)
  colnames(pred.strength) <- c("train", "test")
  for (i in 2:maxclust) {
    rpc.train <- gsl.runt.pruning.criterion(cluster.train)[i-1]
    rpc.test <- gsl.runt.pruning.criterion(cluster.test)[i-1]
    cluster.train.pruned <- gsl.cluster(X.train, mdsfun.train,
                                        pruning.threshold = rpc.train,
                                        pruning.criterion = pruning.criterion)
    cluster.test.pruned <- gsl.cluster(X.test, mdsfun.test,
                                       pruning.threshold = rpc.test,
                                       pruning.criterion = pruning.criterion)
    lab.train.orig <- gsl.observation.labels(cluster.train.pruned)
    lab.test.orig <- gsl.observation.labels(cluster.test.pruned)
    ## Recode labels so that they are between 1 and i
    lab.train <- as.numeric(factor(lab.train.orig, labels = 1:i))
    lab.test <-  as.numeric(factor(lab.test.orig, labels = 1:i))
    lab.train.pred <- as.vector(gsl.knn.classifier(X.test, lab.test, X.train,
                                                   kk = 1)$ypred)
    lab.test.pred <- as.vector(gsl.knn.classifier(X.train, lab.train, X.test,
                                                  kk = 1)$ypred)
    ## browser()
    ## For each test cluster, we compute the proportion of
    ## observation pairs in that cluster that are also assigned to
    ## the same cluster by the training set centroids. The prediction
    ## strength is the minimum of this quantity over the k test
    ## clusters.
    clust.pred.strength <- rep(1, i)
    for (iclust in 1:i) {
      iclust.obs <- (1:ntest)[lab.test == iclust]
      ni <- length(iclust.obs)
      n.same <- 0
      for (k1 in iclust.obs) {
        for(k2 in iclust.obs) {
          if ((k1 != k2) & (lab.test.pred[k1] == lab.test.pred[k2]))
            n.same <- n.same + 1
        }
      }
      clust.pred.strength[iclust] <- n.same / (ni * (ni - 1))
    }
    pred.strength[i, "test"] <- min(clust.pred.strength)
    ## Now exchange of "train" and "test", while we are at it.
    for (iclust in 1:i) {
      iclust.obs <- (1:ntrain)[lab.train == iclust]
      ni <- length(iclust.obs)
      n.same <- 0
      for (k1 in iclust.obs) {
        for(k2 in iclust.obs) {
          if ((k1 != k2) & (lab.train.pred[k1] == lab.train.pred[k2]))
            n.same <- n.same + 1
        }
      }
      clust.pred.strength[iclust] <- n.same / (ni * (ni - 1))
    }
    pred.strength[i, "train"] <- min(clust.pred.strength)
  }
  return(pred.strength)
}

##-----------------------------------------------------------------
## This version assumes that we are using the test graph for the
## entire sample and only use two different density estimates for the
## two half-samples. So we obtain clusterings of the entire sample for
## each of the two half-samples

prediction.strength.cwc <- function(X, mdsfun.train, mdsfun.test,
                                    pruning.criterion, maxclust) {
  test.plot <- function() {
    colors <- c("black", "red", "blue", "green", "yellow", "orange") ###
    mfrow.orig <- par("mfrow")
    par(mfrow = c(2,2))
    plot(X1, type = "n")
    for (j in 1:i) points(X1[lab1 == j,], col = colors[j], pch = 20)
    title("lab1")
    plot(X1, type = "n")
    for (j in 1:i) points(X1[lab1.pred == j,], col = colors[j], pch = 20)
    title("lab1.pred")
    ##
    plot(X2, type = "n")
    for (j in 1:i) points(X2[lab2 == j,], col = colors[j], pch = 20)
    title("lab2")
    plot(X2, type = "n")
    for (j in 1:i) points(X2[lab2.pred == j,], col = colors[j], pch = 20)
    title("lab2.pred")
    par(mfrow = mfrow.orig)
  }
  ##
  n <- nrow(X)
  cluster.train <- gsl.cluster(X, mdsfun.train, assign.fluff = F,
                               pruning.criterion = pruning.criterion)
  cluster.test <- gsl.cluster(X, mdsfun.test, assign.fluff = F,
                              pruning.criterion = pruning.criterion)
  pred.strength <- matrix(1, nrow = maxclust, ncol = 2)
  colnames(pred.strength) <- c("train", "test")
  for (i in 2:maxclust) {
    rpc.train <- gsl.runt.pruning.criterion(cluster.train)[i-1]
    rpc.test <- gsl.runt.pruning.criterion(cluster.test)[i-1]
    cluster.train.pruned <- gsl.cluster(X, mdsfun.train,
                                        pruning.threshold = rpc.train,
                                        pruning.criterion = pruning.criterion)
    cluster.test.pruned <- gsl.cluster(X, mdsfun.test,
                                       pruning.threshold = rpc.test,
                                       pruning.criterion = pruning.criterion)
    in.core.train <- gsl.observation.in.core(cluster.train.pruned)
    in.core.test <- gsl.observation.in.core(cluster.test.pruned)
    in.core.one <- as.logical(in.core.train + in.core.test)
    ## lab.train.orig <- gsl.observation.labels(cluster.train.pruned)
    ## lab.test.orig <- gsl.observation.labels(cluster.test.pruned)
    lab.train.orig <- gsl.observation.labels(cluster.train.pruned)[in.core.one]
    lab.test.orig <- gsl.observation.labels(cluster.test.pruned)[in.core.one]
    ## Recode labels so that they are between 1 and i
    lab.train <- as.numeric(factor(lab.train.orig, labels = 1:i))
    lab.test <-  as.numeric(factor(lab.test.orig, labels = 1:i))
    ## browser()
    ## For each test cluster, we compute the proportion of
    ## observation pairs in that cluster that are also assigned to
    ## the same cluster by the training set centroids. The prediction
    ## strength is the minimum of this quantity over the k test
    ## clusters.
    clust.pred.strength <- rep(1, i)
    for (iclust in 1:i) {
      ## iclust.obs <- (1:n)[lab.test == iclust]
      iclust.obs <- (1:sum(in.core.one))[lab.test == iclust]
      ni <- length(iclust.obs)
      n.same <- 0
      for (k1 in iclust.obs) {
        for(k2 in iclust.obs) {
          if ((k1 != k2) & (lab.train[k1] == lab.train[k2]))
            n.same <- n.same + 1
        }
      }
      clust.pred.strength[iclust] <- n.same / (ni * (ni - 1))
    }
    pred.strength[i, "test"] <- min(clust.pred.strength)
    ## Now exchange of "train" and "test", while we are at it.
    for (iclust in 1:i) {
      ## iclust.obs <- (1:n)[lab.train == iclust]
      iclust.obs <- (1:sum(in.core.one))[lab.train == iclust]
      ni <- length(iclust.obs)
      n.same <- 0
      for (k1 in iclust.obs) {
        for(k2 in iclust.obs) {
          if ((k1 != k2) & (lab.test[k1] == lab.test[k2]))
            n.same <- n.same + 1
        }
      }
      clust.pred.strength[iclust] <- n.same / (ni * (ni - 1))
    }
    pred.strength[i, "train"] <- min(clust.pred.strength)
  }
  return(pred.strength)
}


#################################################################
## Chapter 5. Functions for simulations to assess level and power
##


#################################################################
## Appendix 1: Utilities

###################################################################
## Appendix 2: Modified GSL functions

gsl.pruning.criteria <- function(mast.dendogram, obs.density,
                                 obs.density.lower.confidence.bounds = NULL) {
  MY.INF <- 1.e30
  merge <- mast.dendogram$merge
  height <- mast.dendogram$height
  left.size <- rep(0, nrow(merge))
  right.size <- rep(0, nrow(merge))
  left.excess.mass <- rep(0, nrow(merge))
  right.excess.mass <- rep(0, nrow(merge))
  left.rise <- rep(0, nrow(merge))
  right.rise <- rep(0, nrow(merge))
  left.log.ratio <- rep(0, nrow(merge))
  right.log.ratio <- rep(0, nrow(merge))
  left.bootstrap.rise <- rep(- MY.INF, nrow(merge))
  right.bootstrap.rise <- rep(- MY.INF, nrow(merge))
  left.poisson.gamma <- rep(0, nrow(merge))
  right.poisson.gamma <- rep(0, nrow(merge))
  ##
  ## 7-16-2010
  find.gamma <- function(height.inode, obs.densities) {
    fun <- function(gamma) {
      obs.densities.lower.confidence.bounds <-
        obs.densities - gamma * sqrt(obs.densities)
      height.inode.upper.confidence.bound <-
        height.inode + gamma * sqrt(height.inode)
      return(max(obs.densities.lower.confidence.bounds) -
             height.inode.upper.confidence.bound)
    }
    xmin <- 0
    xmax <- 5
    maxit <- 15
    gamma <- binary.root.search(fun, xmin, xmax, maxit = maxit)
    return(gamma)
  }
  pc.int <- function(inode) {
    ileft <- merge[inode, 1]
    iright <- merge[inode, 2]
    if (ileft < 0) left.leaves <- -ileft
    if (ileft > 0) left.leaves <- pc.int(ileft)
    if (iright < 0) right.leaves <- -iright
    if (iright > 0) right.leaves <- pc.int(iright)
    ## Should this be >=
    left.high.density.cluster <- left.leaves[obs.density[left.leaves] >
                                             height[inode]]
    right.high.density.cluster <- right.leaves[obs.density[right.leaves] >
                                               height[inode]]
    if (length(left.high.density.cluster) > 0) {
      left.size[inode] <<- length(left.high.density.cluster)
      left.excess.mass[inode] <<- sum(1 - height[inode] /
                                      obs.density[left.high.density.cluster])/
                                        (nrow(merge) + 1)
      ##browser()
      left.max <- max(obs.density[left.high.density.cluster])
      if (left.max == Inf) {
        left.rise[inode] <<- Inf
        left.log.ratio[inode] <<- Inf
      }
      else {
        left.rise[inode] <<- left.max - height[inode]
        left.log.ratio[inode] <<- log(left.max) - log(height[inode])
      }
      ##left.poisson.gamma[inode] <<- find.gamma(height[inode],
      ##                                 obs.density[left.high.density.cluster])
      left.poisson.gamma[inode] <<- 0 ## 9-1-2010
      if (!is.null(obs.density.lower.confidence.bounds)) 
        left.bootstrap.rise[inode] <<-
          max(obs.density.lower.confidence.bounds[left.high.density.cluster]) -
            height[inode]
    }
    if (length(right.high.density.cluster) > 0) {
      right.size[inode] <<- length(right.high.density.cluster)
      right.excess.mass[inode] <<- sum(1 - height[inode] /
                                       obs.density[right.high.density.cluster])/
                                         (nrow(merge) + 1)
      right.max <- max(obs.density[right.high.density.cluster])
      if (right.max == Inf) {
        right.rise[inode] <<- Inf
        right.log.ratio[inode] <<- Inf
      }
      else {
        right.rise[inode] <<- right.max - height[inode]
        right.log.ratio[inode] <<- log(right.max) - log(height[inode])
      }
      ##right.poisson.gamma[inode] <<- find.gamma(height[inode],
      ##                                  obs.density[right.high.density.cluster])
      ##right.poisson.gamma[inode] <<- 0
      if (!is.null(obs.density.lower.confidence.bounds))
        right.bootstrap.rise[inode] <<-
          max(obs.density.lower.confidence.bounds[right.high.density.cluster]) -
            height[inode]
    }
    return(c(left.leaves, right.leaves))
  }
  pc.int(nrow(merge))
  res <- cbind(left.size, right.size,
               left.excess.mass, right.excess.mass,
               left.rise, right.rise,
               left.bootstrap.rise, right.bootstrap.rise,
               left.poisson.gamma, right.poisson.gamma,
               left.log.ratio, right.log.ratio
               )
  colnames(res) <- c("left.size", "right.size", "left.excess.mass", "right.excess.mass",
                     "left.rise", "right.rise", "left.bootstrap.rise",
                     "right.bootstrap.rise", "left.poisson.gamma",
                     "right.poisson.gamma", "left.log.ratio", "right.log.ratio")
  return(res)
}

##-----------------------------------------------------------------
## Modifcation (2-26-2011)
## Introduce new pruning criterions "size.with.tiebreak" so that we can
## always obtain any number of clusters between 1 and n when we use
## the nn density estimate and "size" as the pruning criterion.
##
## Modification (6-30-2011)
## Introduce runt.log.ratio as a pruning criterion

gsl.compute.cluster.tree <- function(dendogram, gsl.pc.out, pruning.criterion = "size",
                                     pruning.threshold = 0) {
  runt.size <- pmin(gsl.pc.out[,"left.size"], gsl.pc.out[,"right.size"])
  runt.excess.mass <- pmin(gsl.pc.out[,"left.excess.mass"],
                           gsl.pc.out[,"right.excess.mass"])
  runt.rise <- pmin(gsl.pc.out[,"left.rise"], gsl.pc.out[,"right.rise"])
  runt.bootstrap.rise <- pmin(gsl.pc.out[,"left.bootstrap.rise"],
                              gsl.pc.out[,"right.bootstrap.rise"])
  runt.poisson.gamma <- pmin(gsl.pc.out[,"left.poisson.gamma"],
                             gsl.pc.out[,"right.poisson.gamma"])
  runt.log.ratio <- pmin(gsl.pc.out[,"left.log.ratio"], gsl.pc.out[,"right.log.ratio"])
  if (pruning.criterion == "size.with.tiebreak") 
    runt.crit <- runt.size + seq(0, 0.1, length.out = length(runt.size))
  if (pruning.criterion == "size") runt.crit <- runt.size
  if (pruning.criterion == "excess.mass") runt.crit <- runt.excess.mass
  if (pruning.criterion == "rise") runt.crit <- runt.rise
  if (pruning.criterion == "bootstrap.rise") runt.crit <-
    runt.bootstrap.rise
  if (pruning.criterion == "poisson.gamma") runt.crit <- runt.poisson.gamma
  if (pruning.criterion == "log.ratio") runt.crit <- runt.log.ratio
  ##browser()
  merge <- dendogram$merge
  level <- dendogram$height
  nobs <- nrow(merge) + 1

  cluster.tree <- matrix(0, nrow = 2 * nobs - 1, ncol = 18)
  cn <- c("leftson", "rightson", "runt.size", "runt.excess.mass",
          "leaf.code", "level", "dendogram.node", "size",
          "excess.mass", "rise", "bootstrap.rise",
          "poisson.gamma", "log.ratio", "runt.rise", "runt.bootstrap.rise",
          "runt.poisson.gamma", "runt.log.ratio", "runt.pruning.crit")
  colnames(cluster.tree) <- cn
  nnodes <- 1
  cluster.tree[1, "size"] <- nobs
  cluster.tree[1, "excess.mass"] <- 1
  cluster.tree[1, "level"] <- 0
  cluster.tree[1, "leaf.code"] <- 1
  cluster.tree[1, "dendogram.node"] <- nrow(merge)

  large.runt <- (runt.crit >= pruning.threshold) & (runt.rise > 0)
  ## large.runt <- (runt.crit >= pruning.threshold)

  find.large.runt.node.in.subdend <- function(root) {
    ## Note: Because of the monotonicity of the pruning criteria either
    ## the root has large runt or at most one of the two subdendograms.
    if (root < 0) return(0)
    if (large.runt[root]) return(root)
    large.runt.node <- 0
    if (merge[root, 1] > 0)
      large.runt.node <- find.large.runt.node.in.subdend(merge[root, 1])
    if ((large.runt.node == 0) & (merge[root, 2] > 0))
      large.runt.node <-
        find.large.runt.node.in.subdend(merge[root, 2])
    return(large.runt.node)
  }
  make.cluster.tree.sons <- function(parent) {
    ## browser()
    if (cluster.tree[parent, "dendogram.node"] < 0) return()
    large.runt.node <-
      find.large.runt.node.in.subdend(cluster.tree[parent, "dendogram.node"])
    if (large.runt.node > 0) {
      cluster.tree[parent, "dendogram.node"] <<- large.runt.node
      cluster.tree[parent, "runt.size"] <<- runt.size[large.runt.node]
      cluster.tree[parent, "runt.excess.mass"] <<- runt.excess.mass[large.runt.node]
      cluster.tree[parent, "runt.rise"] <<- runt.rise[large.runt.node]
      cluster.tree[parent, "runt.bootstrap.rise"] <<- runt.bootstrap.rise[large.runt.node]
      cluster.tree[parent, "runt.poisson.gamma"] <<- runt.poisson.gamma[large.runt.node]
      cluster.tree[parent, "runt.pruning.crit"] <<- runt.crit[large.runt.node]
      cluster.tree[parent, "runt.log.ratio"] <<- runt.log.ratio[large.runt.node]
      nnodes <<- nnodes + 1
      cluster.tree[parent, "leftson"] <<- nnodes
      cluster.tree[nnodes, "leaf.code"] <<- 2 * cluster.tree[parent, "leaf.code"]
      cluster.tree[nnodes, "size"] <<- gsl.pc.out[large.runt.node, "left.size"]
      cluster.tree[nnodes, "excess.mass"] <<-
        gsl.pc.out[large.runt.node, "left.excess.mass"]
      cluster.tree[nnodes, "rise"] <<-
        gsl.pc.out[large.runt.node, "left.rise"]
      cluster.tree[nnodes, "bootstrap.rise"] <<-
        gsl.pc.out[large.runt.node, "left.bootstrap.rise"]
      cluster.tree[nnodes, "log.ratio"] <<-
        gsl.pc.out[large.runt.node, "left.log.ratio"]
      cluster.tree[nnodes, "poisson.gamma"] <<-
        gsl.pc.out[large.runt.node, "left.poisson.gamma"]
      cluster.tree[nnodes, "level"] <<- level[large.runt.node]
      cluster.tree[nnodes, "dendogram.node"] <<- merge[large.runt.node, 1]
      left.parent <- nnodes
      make.cluster.tree.sons(left.parent)
      nnodes <<- nnodes + 1
      cluster.tree[parent, "rightson"] <<- nnodes
      cluster.tree[nnodes, "leaf.code"] <<- 2 * cluster.tree[parent, "leaf.code"] + 1
      cluster.tree[nnodes, "size"] <<- gsl.pc.out[large.runt.node, "right.size"]
      cluster.tree[nnodes, "excess.mass"] <<-
        gsl.pc.out[large.runt.node, "right.excess.mass"]
      cluster.tree[nnodes, "rise"] <<-
        gsl.pc.out[large.runt.node, "right.rise"]
      cluster.tree[nnodes, "bootstrap.rise"] <<-
        gsl.pc.out[large.runt.node, "right.bootstrap.rise"]
      cluster.tree[nnodes, "poisson.gamma"] <<-
        gsl.pc.out[large.runt.node, "right.poisson.gamma"]
      cluster.tree[nnodes, "log.ratio"] <<-
        gsl.pc.out[large.runt.node, "right.log.ratio"]
      cluster.tree[nnodes, "level"] <<- level[large.runt.node]
      cluster.tree[nnodes, "dendogram.node"] <<- merge[large.runt.node, 2]
      right.parent <- nnodes
      make.cluster.tree.sons(right.parent)
    }
    return()
  }
  make.cluster.tree.sons(1)
  if (nnodes == 1) {
    ct.out <- matrix(cluster.tree[1,], nrow = 1, byrow = T)
    colnames(ct.out) <- cn
    return(ct.out)
  }
  return(cluster.tree[1:nnodes,])
}

##-----------------------------------------------------------------
## Change runt functions so that they only return stats for internal
## nodes

gsl.runt.rise <- function(gsl.cluster.out) {
  runt.size <- gsl.cluster.out$cluster.tree[, "runt.size"]
  return(rev(sort(runt.size[runt.size > 0])))
}

gsl.runt.excess.mass <- function(gsl.cluster.out) {
  runt.size <- gsl.cluster.out$cluster.tree[, "runt.size"]
  runt.excess.mass <- gsl.cluster.out$cluster.tree[, "runt.excess.mass"]
  return(rev(sort(runt.excess.mass[runt.size > 0])))
}

gsl.runt.rise <- function(gsl.cluster.out) {
  runt.rise <- gsl.cluster.out$cluster.tree[, "runt.rise"]
  runt.size <- gsl.cluster.out$cluster.tree[, "runt.size"]
  return(rev(sort(runt.rise[runt.size > 0])))
}

gsl.runt.bootstrap.rise <- function(gsl.cluster.out) {
  runt.bootstrap.rise <- gsl.cluster.out$cluster.tree[, "runt.bootstrap.rise"]
  runt.size <- gsl.cluster.out$cluster.tree[, "runt.size"]
  return(rev(sort(runt.bootstrap.rise[runt.size > 0])))
}


##-----------------------------------------------------------------

gsl.runt.stats <- function(gsl.cluster.out, order.by = "nothing") {
  ct <- gsl.cluster.out$cluster.tree
  nnodes <- nrow(ct)
  rs.order <- 1:nnodes
  if (order.by == "size") rs.order <- rev(order(ct[, "runt.size"]))
  if (order.by == "excess.mass") rs.order <- rev(order(ct[, "runt.excess.mass"]))
  if (order.by == "rise") rs.order <- rev(order(ct[, "runt.rise"]))
  if (order.by == "bootstrap.rise") rs.order <- rev(order(ct[, "runt.bootstrap.rise"]))
  runt.stats <- matrix(0, nrow = nrow(ct), ncol = 4)
  runt.stats[,1] <- ct[rs.order, "runt.size"]
  runt.stats[,2] <- ct[rs.order, "runt.excess.mass"]
  runt.stats[,3] <- ct[rs.order, "runt.rise"]
  runt.stats[,4] <- ct[rs.order, "runt.bootstrap.rise"]
  colnames(runt.stats) <- c("runt.size", "runt.excess.mass", "runt.rise",
                            "runt.bootstrap.rise")
  ## interior <- runt.stats[, "runt.size"] > 0
  ## return(runt.stats[interior,])
  return(runt.stats)
}

##-----------------------------------------------------------------
## Draw possibly annotated cluster tree

gsl.draw.tree <- function(cluster.tree, draw.cluster.numbers = F,
                          draw.cluster.leaf.codes = F,
                          draw.runt.excess.mass = F,
                          draw.runt.size = F, interior.annotations = NULL) {
  spacing <- 1
  nnodes <- nrow(cluster.tree)
  left.width <- rep(0, nnodes)
  right.width <- rep(0, nnodes)
  node.x <- rep(0, nnodes)
  node.y <- cluster.tree[, "level"]
  runt.excess.mass <- cluster.tree[, "runt.excess.mass"]
  runt.size <- cluster.tree[, "runt.size"]
  n <- cluster.tree[1, "size"]
  plot.symbol <- rep(0, nnodes)
  current <- 1
  for (i in 1:nnodes) {
    if (cluster.tree[i, 1] == 0) {
      if (draw.cluster.numbers) {
        plot.symbol[i] <- toString(current)
        current <- current + 1
      }
      if (draw.cluster.leaf.codes) 
        plot.symbol[i] <- toString(cluster.tree[i, "leaf.code"])
    }
  }
  is.leaf <- cluster.tree[, 1] == 0
  compute.width <- function(i, spacing) {
    left.width[i] <<- 0
    right.width[i] <<- 0
    if(cluster.tree[i, 1] > 0) {
      left.width[i] <<- compute.width(cluster.tree[i, 1], spacing)
      right.width[i] <<- compute.width(cluster.tree[i, 2], spacing)
    }
    return(left.width[i] + right.width[i] + spacing)
  }
  compute.node.xy <- function(i, lr, mother.x, mother.y, spacing) {
    if(lr == "left") {
      node.x[i] <<- mother.x - (right.width[i] + 0.5 * spacing)
    }
    else {
      node.x[i] <<- mother.x + (left.width[i] + 0.5 * spacing)
    }
    if (cluster.tree[i, 1] > 0) {
      compute.node.xy(cluster.tree[i, 1], "left", node.x[i], node.y[i], spacing)
      compute.node.xy(cluster.tree[i, 2], "right", node.x[i], node.y[i], spacing)
    }
  }
  draw.edges <- function(i) {
    if(cluster.tree[i, 1] > 0) {
      leftson <- cluster.tree[i, 1]
      rightson <- cluster.tree[i, 2]
      lines(c(node.x[i], node.x[leftson]), c(node.y[i], node.y[leftson]))
      lines(c(node.x[i], node.x[rightson]), c(node.y[i], node.y[rightson]))
      draw.edges(leftson)
      draw.edges(rightson)
    }
  }
  compute.width(1, spacing)
  compute.node.xy(1, "left", 0, 0, spacing)
  ylim = c(-0.05 * max(node.y), 1.05 * max(node.y))
  if (draw.cluster.numbers | draw.cluster.leaf.codes) {
    plot(node.x, node.y, pch = 19, xaxt = "n", xlab = "", ylab = "level",
         ylim = ylim, cex.axis = 1.5, cex.lab = 1.5)
    text(node.x[is.leaf], node.y[is.leaf] + 0.05 * max(node.y),
         labels = plot.symbol[is.leaf], cex = 1.5)
    points(node.x[!is.leaf], node.y[!is.leaf], pch = 19)
  }
  else {
    plot(node.x, node.y, pch = 19, xaxt = "n", xlab = "", ylab = "level", ylim = ylim,
         cex.axis = 1.5, cex.lab = 1.5)
  }
  if(draw.runt.excess.mass) {
    ##    text(node.x[!is.leaf] + 0.2, node.y[!is.leaf],
    ##       round(n * runt.excess.mass[!is.leaf], 0))
    text(node.x[!is.leaf], node.y[!is.leaf] - 0.05 * max(node.y),
         round(n * runt.excess.mass[!is.leaf], 0), cex = 1.5)
  }
  if(draw.runt.size) {
    text(node.x[!is.leaf], node.y[!is.leaf] - 0.05 * max(node.y),
         runt.size[!is.leaf], cex = 1.5)
  }
  if(!is.null(interior.annotations)) {
    text(node.x[!is.leaf], node.y[!is.leaf] - 0.05 * max(node.y),
         interior.annotations[!is.leaf], cex = 1.5)
  }
  draw.edges(1)
  return(list(x = node.x, y = node.y))
}

##-----------------------------------------------------------------

rand.index <- function(true.labels, estimated.labels) {
  tab <- table(true.labels, estimated.labels)
  nidot <- apply(tab, 1, sum)
  ndotj <- apply(tab, 2, sum)
  n <- sum(tab)
  num <- sum(choose(tab, 2)) - 
    sum(choose(nidot, 2)) * sum(choose(ndotj, 2)) / choose(n, 2)
  denom <- 0.5 * (sum(choose(nidot, 2)) + sum(choose(ndotj, 2))) -
    sum(choose(nidot, 2)) * sum(choose(ndotj, 2)) / choose(n, 2)
  return(num/denom)
}

##-----------------------------------------------------------------


#################################################################
#################################################################
#################################################################

##-----------------------------------------------------------------
## Compute bootstrap vertex and edge weights using nearest neighbor density
## estimate

## cwc.make.nn.density.estimate <- function(X.train) {
##   density.estimate <- function(X.eval) {
##      phat <- cwc.nn.density.estimate(X.train, X.eval)
##      return(phat)
##    }
##   return(density.estimate)
## }

## cwc.nn.density.estimate <- function(X.obs, X.eval) {
##   m <- ncol(X.obs)
##   n <- nrow(X.obs)
##   if (is.vector(X.obs)) {
##     X.obs <- matrix(X.obs, ncol = 1)
##     X.eval <- matrix(X.eval, ncol = 1)
##   }
##   n.eval <- nrow(X.eval)
##   n.obs <- nrow(X.obs)
##   p <- ncol(X.obs)
##   dens <- rep(0, n.eval)
##   dist <- sqrt(gsl.interpoint.distance.matrix(X.eval, X.obs)) 
##   for (i in 1:n.eval) {
##     r <- min(dist[i,])
##     if (r == 0) dens[i] <- Inf else dens[i] <- 1 / (n * r^m)
##   }
##   return(dens)
## }

## ## -----------------------------------------------------------------

## cwc.boot.nn.ve.weights <- function(tg.vertices, tg.edges, X, nboot,
##                                    ngrid = 10, interactive = T) {
##   n.vertices <- nrow(tg.vertices)
##   n.edges <- nrow(tg.edges)
##   n <- nrow(X)
##   bs.size <- round(nrow(X) / 2)
##   bootstrap.vertex.weights <- matrix(0, nrow = n.vertices, ncol = nboot)
##   bootstrap.edge.weights <- matrix(0, nrow = n.edges, ncol = nboot)
##   bootstrap.samples <- matrix(0, nrow = bs.size, ncol = nboot)
##   if (interactive) cat("\n\nComputing bootstrap vertex and edge weights\n")
##   for (i in 1:nboot) {
##     bs <- sample(1:n, bs.size)
##     bootstrap.samples[,i] <- bs
##     X.boot <- as.matrix(X[bs, ], nrow = bs.size)
##     density <- cwc.make.nn.density.estimate(X.boot)
##     cwcvw.out <- cwc.ve.weights(tg.vertices, tg.edges, density, ngrid)
##     bootstrap.vertex.weights[,i] <- cwcvw.out$vertex.weights
##     bootstrap.edge.weights[,i] <- cwcvw.out$edge.weights
##     if(interactive) cat(" ", i)
##   }
##   return(list(bootstrap.vertex.weights = bootstrap.vertex.weights,
##               bootstrap.edge.weights = bootstrap.edge.weights,
##               bootstrap.samples = bootstrap.samples,
##               bootstrap.h = NULL,
##               cvs.out = NULL))
## }

## ## -----------------------------------------------------------------
## ## Make mdsfun from edges, edge weights, end vertex weights. First
## ## a simple version that uses a list to store the weights. Convert the
## ## vertex indices of a list into a number

## cwc.make.mdsfun.simple <- function(edges, edge.weights, vertex.weights) {
##   n <- length(vertex.weights)
##   ne <- length(edge.weights)
##   mult <- 10^(ceiling(log10(n + 1)))
##   weight.list <- list()
##   convert <- function(i, j) {
##     return(mult * i + j)
##   }
##   for (i in 1:n) weight.list[[convert(i, i)]] <- vertex.weights[i]
##   for (i in 1:ne)
##     weight.list[[convert(edges[i, 1], edges[i, 2])]] <- edge.weights[i]
##   mdsfun <- function(i, js) {
##     nj <- length(js)
##     val <- rep(0, nj)
##     for (k in 1:nj) {
##       key <- convert(i, js[k])
##       val[k] <- weight.list[[key]]
##     }
##     return(val)
##   }
##   attr(mdsfun, "n") <- n
##   return(mdsfun)
## }
  
## ## =================================================================
## ## Modified versions of gsl functions to deal with additional argument
## ## obs.density.lower.confidence.bounds

## gsl.cluster <- function(X, mdsfun, pruning.criterion = "size",
##                         pruning.threshold = 0, gsl.cluster.out = NULL,
##                         assign.fluff = T,
##                         obs.density.lower.confidence.bounds = NULL) {
##   n <- attr(mdsfun, "n")
##   if (is.null(gsl.cluster.out)) {
##     mast <- gsl.mast(mdsfun)
##     mast.dendogram <- gsl.mast.dendogram(mast)
##     obs.density <- rep(0, n)
##     for (i in 1:n) obs.density[i] <- mdsfun(i, i)
##     pc.values <- gsl.pruning.criteria(mast.dendogram, obs.density,
##                                       obs.density.lower.confidence.bounds =
##                                       obs.density.lower.confidence.bounds)
##   }
##   else {
##     mast.dendogram <- gsl.cluster.out$mast.dendogram
##     pc.values <- gsl.cluster.out$pc.values
##     obs.density <- gsl.cluster.out$obs.density
##   }
##   cluster.tree <- gsl.compute.cluster.tree(mast.dendogram, pc.values,
##                                    pruning.criterion = pruning.criterion,
##                                    pruning.threshold = pruning.threshold)
##   afo.out <- NULL
##   if (assign.fluff) afo.out <- gsl.assign.fluff.oneshot(X, cluster.tree,
##                                                         mast.dendogram, obs.density)
##   return(list(cluster.tree = cluster.tree, leaf.code = afo.out$leaf.code,
##               cluster.core = afo.out$cluster.core, mast.dendogram = mast.dendogram,
##               pc.values = pc.values, obs.density = obs.density))
## }

##=================================================================
## Use versions in cwc-functions-7-17-2010.R

## cwc.ve.weights <- function(tg.vertices, tg.edges, density, ngrid) {
##   n.edges <- nrow(tg.edges)
##   m <- ncol(tg.vertices)
##   edge.weights <- rep(0, n.edges)
##   vertex.weights <- density(tg.vertices)
##   t <- seq(0, 1, length = ngrid)[2:(ngrid-1)]
##   Tmat <- matrix(t, nrow = ngrid-2, ncol = m)
##   for (i in 1:n.edges) {
##     v1.index <- tg.edges[i, 1]
##     v2.index <- tg.edges[i, 2]
##     v1 <- tg.vertices[v1.index,]
##     v2 <- tg.vertices[v2.index,]
##     V1 <- matrix(v1, nrow = ngrid-2, ncol = m, byrow = T)
##     V2 <- matrix(v2, nrow = ngrid-2, ncol = m, byrow = T)
##     E <- (1-Tmat) * V1 + Tmat * V2
##     phat <- c(vertex.weights[v1.index], density(E), vertex.weights[v2.index])
##     edge.weights[i] <- min(phat)
##   }
##   return(list(vertex.weights = vertex.weights, edge.weights = edge.weights))
## }

## ##-----------------------------------------------------------------
## ## Compute bootstrap vertex and edge weights using kernel density estimate
## ## with LS CV.
## ## 
## ## In current version, bootstrap samples are generated by half-sampling
## ## May want to include bootstrap sampling from kernel density estimate
## ## for original sample.
## ##
## ## Choices for bandwidth:
## ## "cv.mean"   average CV bandwidths over bootstrap samples
## ## "cv"        individual CV bandwidth for each bootstrap sample
## ## a number    use this bandwidth for all bootstrap samples


## ## Version of 5-19-2010 below: Allow for bootstrap or half-sampling
## ## Cross validation fails for bootstrap resamples, presumably because
## ## of ties. Need to work out weighted version of LS CV??
## ## For the moment implement simple and not entirely kosher version:
## ## Use estimated bandwidth for original sample on all the bootstrap
## ## samples

## cwc.boot.kernel.ve.weights <- function(tg.vertices, tg.edges, X, nboot = NULL, 
##                                        bootstrap.samples = NULL, trial.h = NULL,
##                                        ngrid = 10, interactive = F,
##                                        bandwidth = "cv.mean",
##                                        resampling.type = "half.sampling") {
##   if (is.null(trial.h))
##     trial.h <- c(0.01, 0.05, 0.10, 0.15, 0.20, 0.30, 0.40, 0.60, 0.80,
##                  1.00, 1.50, 2.00)
##   n.vertices <- nrow(tg.vertices)
##   n.edges <- nrow(tg.edges)
##   n <- nrow(X)
##   if (is.null(bootstrap.samples)) {
##     if (resampling.type == "half.sampling") bs.size <- round(nrow(X) / 2)
##     if (resampling.type == "bootstrap") bs.size <- nrow(X)
##     bootstrap.samples <- matrix(0, nrow = bs.size, ncol = nboot)
##     for (i in 1:nboot) {
##       if (resampling.type == "half.sampling") bs <- sample(1:n, bs.size, replace = F)
##       if (resampling.type == "bootstrap") bs <- sample(1:n, bs.size, replace = T)
##       bootstrap.samples[,i] <- bs
##     }
##   }
##   nboot <- ncol(bootstrap.samples)
##   resampling.type <- "half.sampling"
##   if (nrow(bootstrap.samples) == n) resampling.type <- "bootstrap"
##   bootstrap.vertex.weights <- matrix(0, nrow = n.vertices, ncol = nboot)
##   bootstrap.edge.weights <- matrix(0, nrow = n.edges, ncol = nboot)
##   cvs.out <- vector("list", nboot)
##   bootstrap.h <- rep(0, nboot)

##   if ((resampling.type == "bootstrap") & !is.numeric(bandwidth)) {
##     cvso <- cv.search(X, trial.h)
##     cvs.out[[1]] <- cvso
##     bootstrap.h[1] <- cvso$opt.smopar
##     if (interactive) cat("\n", bootstrap.h[1])
##   } 
##   if ((resampling.type == "half.sampling") & !is.numeric(bandwidth)) {
##     if (interactive) cat("\n\nFinding bandwidths")
##     for (i in 1:nboot) {
##       bs <- bootstrap.samples[,i]
##       X.boot <- as.matrix(X[bs, ], nrow = bs.size)
##       cvso <- cv.search(X.boot, trial.h)
##       cvs.out[[i]] <- cvso
##       bootstrap.h[i] <- cvso$opt.smopar
##       if (interactive) cat("\n", i, bootstrap.h[i])
##     }
##     mean.h <- mean(bootstrap.h[!is.na(bootstrap.h)])
##   }
##   if (interactive) cat("\n\nComputing bootstrap vertex and edge weights\n")
##   for (i in 1:nboot) {
##     if (resampling.type == "bootstrap") h <- bootstrap.h[1]
##     if (resampling.type == "half.sampling") {
##       h <- mean.h
##       if ((bandwidth == "cv") & (!is.na(bootstrap.h[i]))) h <- bootstrap.h[i]
##     }
##     if (is.numeric(bandwidth)) h <- bandwidth
##     bs <- bootstrap.samples[,i]
##     X.boot <- as.matrix(X[bs, ], nrow = bs.size)
##     density <- make.gaussian.kernel.density.estimate(X.boot, h)
##     cwcvw.out <- cwc.ve.weights(tg.vertices, tg.edges, density, ngrid)
##     ##  browser()
##     bootstrap.vertex.weights[,i] <- cwcvw.out$vertex.weights
##     bootstrap.edge.weights[,i] <- cwcvw.out$edge.weights
##     if(interactive) cat(" ", i)
##   }
##   return(list(tg.vertices = tg.vertices,
##               tg.edges = tg.edges,
##               bootstrap.vertex.weights = bootstrap.vertex.weights,
##               bootstrap.edge.weights = bootstrap.edge.weights,
##               bootstrap.samples = bootstrap.samples,
##               bootstrap.h = bootstrap.h,
##               cvs.out = cvs.out))
             
## }

##-----------------------------------------------------------------

## cwc.make.bootstrap.kernel.mdsfun <- function(tg.edges, bootstrap.vertex.weights,
##                                          bootstrap.edge.weights) {
##   n.boot <- ncol(bootstrap.vertex.weights)
##   n.edges <- nrow(tg.edges)
##   n <- nrow(bootstrap.vertex.weights)
##   bootstrap.kernel.mdsfun <- function(i, j) {
##     if (i == j) return(bootstrap.vertex.weights[i,])
##     for (k in 1:n.edges) {
##       if ((tg.edges[k, 1] == i) & (tg.edges[k, 2] == j))
##         return(bootstrap.edge.weights[k,])
##       if ((tg.edges[k, 1] == j) & (tg.edges[k, 2] == i))
##         return(bootstrap.edge.weights[k,])
##     }
##     return(rep(0, n.boot))
##   }
##   attr(bootstrap.kernel.mdsfun, "nboot") <- n.boot
##   attr(bootstrap.kernel.mdsfun, "n") <- n
##   return(bootstrap.kernel.mdsfun)
## }

##-----------------------------------------------------------------
## kmax = # of nearest neighbors, not including the point itself.

## cwc.determine.edges <- function(X, kmax = 20, include.mst = T,
##                                 interactive = F) {
##   if (is.vector(X)) X <- matrix(X, ncol = 1)
##   n <- nrow(X)
##   p <- ncol(X)
##   already.computed <- matrix(F, nrow = n, ncol = n)
##   edges <- matrix(0, nrow = n * kmax + (n-1), ncol = 2)
##   nedges <- 0
##   Dist <- gsl.interpoint.distance.matrix(X, X)
##   if (include.mst) {
##     mst.mat.out <- gsl.mst.mat(Dist)
##     mst.edges <- mst.mat.out$edges
##     for (i in 1:(n-1)) {
##       nedges <- nedges + 1
##       edges[nedges,] <- mst.edges[i,]
##       v1 <- mst.edges[i, 1]
##       v2 <- mst.edges[i, 2]
##       already.computed[v1, v2] <- T
##       already.computed[v2, v1] <- T
##     }
##   }
##   if(kmax >= 1) {
##     for (i in 1:n) {
##       v1 <- i
##       NN <- order(Dist[i,])
##       for (j in 2:(kmax + 1)) {
##         v2 <- NN[j]
##         if (!already.computed[v1, v2]) {
##           nedges <- nedges + 1
##           edges[nedges, 1] <- v1
##           edges[nedges, 2] <- v2
##           already.computed[v1, v2] <- T
##           already.computed[v2, v1] <- T
##         }
##       }
##     }
##   }
##   return(edges[1:nedges, ])
## }




##=================================================================
## A sparse matrix hack. Conclusion: It's too slow

## cwc.make.sparse.matrix <- function(n.row) {
##   mat <- vector(mode = "list", length = n.row)
##   access.mat <- function(i, j, value = NULL) {
##     row.entries <- mat[[i]]
##     if (is.null(value)) {
##        for (k in 1:length(row.entries)) {
##         if (row.entries[[k]][[1]] == j) return(row.entries[[k]][[2]])
##       }
##       return(NULL)
##     }
##     if (!is.null(row.entries)) {
##       for (k in 1:length(row.entries)) {
##         if(row.entries[[k]][[1]] == j) {
##           row.entries[[k]][[2]] <- value
##           return(value)
##         }
##       }
##     }
##     row.entries[[length(row.entries) + 1]] <- vector("list", 2)
##     row.entries[[length(row.entries)]][[1]] <- j
##     row.entries[[length(row.entries)]][[2]] <- value
##     mat[[i]] <<- row.entries
##     return(value)
##   }
##   return(access.mat)
## }

## ## N <- 1000
## ## bla <- cwc.make.sparse.matrix(N)
## ## for (i in 1:N) {
## ##   for (j in 1:20) bla(i, j, i+j)
## ## }

## ## for (i in 1:N) {
## ##   for (j in 1:N) bla(i, j)
## ## }

                            



#################################################################
#################################################################
#################################################################
## Old stuff
##-----------------------------------------------------------------
## Compare three ways of finding upper confidence bounds for edge
## weights and lower confidence bounds for vertex weights"
##
## 1. Empirical upper resp.lower bounds
## 2. Bounds based on Gaussian approximation
## 3. Bounds constructed based on Poisson approximation
##    [phat - gamma.minus * sqrt(phat), phat + gamma.plus * sqrt(phat)]
##
## Option (3)
## If y_i are the confidence bounds and phat_i are the weights, then
## the LS estimate for gamma is
##
## gamma = [sum(y_i sqrt(phat_i)) - sum(phat_i^{3/2)] / sum(phat_i)


## cwc.cowc.confimpare.confidence.bounds <- function(vw, ew, rvw, rew,
##                                           confidence.level = 0.95) {
##   par(mfrow = c(2,1))
##   n.resamples <- ncol(rew)
##   ne <- length(ew)
##   nv <- length(vw)
##   z.cov <- qnorm(1 - (1-confidence.level)/2)

##   nout <- ceiling(n.resamples * (1-confidence.level)/2)
##   lower <- nout + 1
##   upper <- n.resamples - nout
##   if (lower >= upper) {
##     lower <- lower - 1
##     upper <- upper + 1
##   }
##   upper.empirical.ew <- rep(0, ne)
##   lower.empirical.vw <- rep(0, nv)
##   for (i in 1:ne) {
##     sew <- sort(rew[i,])
##     upper.empirical.ew[i] <- sew[upper]
##   }
##   for (i in 1:nv) {
##     svw <- sort(rvw[i,])
##     lower.empirical.vw[i] <- svw[lower]
##   }
##   ew.mean <- apply(rew, 1, mean)
##   ew.sd <- apply(rew, 1, sd)
##   vw.mean <- apply(rvw, 1, mean)
##   vw.sd <- apply(rvw, 1, sd)
##   upper.gaussian.ew <- ew.mean + z.cov * ew.sd
##   lower.gaussian.vw <- vw.mean - z.cov * vw.sd

##   gamma.plus <- (sum(upper.empirical.ew * sqrt(ew.mean)) -
##                  sum(ew.mean^(3/2))) / sum(ew.mean)
##   gamma.minus <- (sum(lower.empirical.vw * sqrt(vw.mean)) -
##                   sum(vw.mean^(3/2))) / sum(vw.mean)
##   upper.poisson.ew <- ew.mean + gamma.plus * sqrt(ew.mean)
##   lower.poisson.vw <- vw.mean + gamma.minus * sqrt(vw.mean)

##   plot(ew.mean, upper.empirical.ew - ew.mean, pch = 20, cex = 0.2)
##   points(ew.mean, upper.gaussian.ew - ew.mean, pch = 20, cex = 0.2, col = "orange")
##   oew <- order(ew.mean)
##   lines(ew.mean[oew], upper.poisson.ew[oew] - ew.mean[oew], col = "red", lwd = 3)
##   title("Upper confidence bounds for edge weights")
##   plot(vw.mean, lower.empirical.vw - vw.mean, pch = 20, cex = 0.2)
##   points(vw.mean, lower.gaussian.vw - vw.mean, pch = 20, cex = 0.2, col = "orange")
##   ovw <- order(vw.mean)
##   lines(vw.mean[ovw], lower.poisson.vw[ovw] - vw.mean[ovw], col = "red", lwd = 3)
##   title("Lower confidence bounds for vertex weights")
##   par(mfrow = c(1,1))

##   nin.empirical <- n.resamples
##   nin.gaussian <- n.resamples
##   nin.poisson <- n.resamples
##   for (i in 1:n.resamples) {
##     out.empirical <- sum(rew[,i] > upper.empirical.ew) +
##       sum(rvw[,i] < lower.empirical.vw)
##     out.gaussian <- sum(rew[,i] > upper.gaussian.ew) +
##       sum(rvw[,i] < lower.gaussian.vw)
##     out.poisson <- sum(rew[,i] > upper.poisson.ew) +
##       sum(rvw[,i] < lower.poisson.vw)
##     if (out.empirical > 0) nin.empirical <- nin.empirical - 1
##     if (out.gaussian > 0) nin.gaussian <- nin.gaussian - 1
##     if (out.poisson > 0) nin.poisson <- nin.poisson - 1
##   }
##   return(list(simul.empirical.cov.prob = nin.empirical / n.resamples,
##               simul.gaussian.cov.prob = nin.gaussian / n.resamples,
##               simul.poisson.cov.prob = nin.poisson / n.resamples))
## }
  

