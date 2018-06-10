## Functions for computing Walther similarity matrix (5-17-07)
## ===========================================================
##
## Revised 11-18-2016}
## ===================

## There are two versions of granulometric level set estimates.
## in both versions, the estimated level set at level lambda is the
## union of spheres around obs x_i with phat(x_i) >= lambda.
##
## In the version described in Walther's 1997 Annals paper, the radius
## of the spheres is fixed, and the sphere around x_i is included if
## is does not contain any obs with phat(xi) < lambda.
##
## In the variable radius version, all obs with phat(x_i) >= lambda
## are included, and the radius of the sphere around x_i is the
## distance to the closest obs with phat < lambda.

## The code below is for the variable radius version which does not
## require choosing a radius r.


compute.radii <- function(X, phat, lambda) {
  n <- nrow(X)
  D2 <- gsl.interpoint.distance.matrix(X, X)
  D2[D2 < 0] <- 0
  D <- sqrt(D2)
  D.aug <- cbind(D, rep(max(D), n))
  phat.aug <- c(phat, 0)
  radii  <- rep(0, n)
  for (i in 1:n) {
    di <- D.aug[i, ]
    radii[i] <- min(di[phat.aug <= lambda])
  }
  return(radii)
}

## This function has a bug
##
## compute.radii.interp <- function(X, phat, lambda, interpolation = F) {
##   argmin <- function(x) {
##     return(((1:length(x))[x == min(x)])[1])
##   }
##   n <- nrow(X)
##   D2 <- gsl.interpoint.distance.matrix(X, X)
##   D2[D2 < 0] <- 0
##   D <- sqrt(D2)
##   D.aug <- cbind(D, rep(max(D), n))
##   phat.aug <- c(phat, 0)
##   radii  <- rep(0, n)
##   for (i in 1:n) {
##     di <- D.aug[i,]
##     j <- argmin(di[phat.aug <= lambda])
##     r <- di[j]
##     radii[i] <- r
##     rstar <- min(di[phat.aug <= lambda])
##     if (r != rstar) browser()
##     if (interpolation) {
##       delta.phat <- phat[i] - phat[j]
##       if (abs(delta.phat) > 0) {
##         a <- (phat[i] - lambda) / delta.phat
##         radii[i] <- a * r
##       }
##     }
##   }
##   return(radii)
## }

## Revised version (6-6-2018)
##
##
compute.radii.interp <- function(X, phat, lambda, interpolation = F) {
  n <- nrow(X)
  D2 <- gsl.interpoint.distance.matrix(X, X)
  D2[D2 < 0] <- 0
  D <- sqrt(D2)
  D.aug <- cbind(D, rep(max(D), n))
  phat.aug <- c(phat, 0)
  radii  <- rep(0, n)
  for (i in 1:n) {
    di <- D.aug[i,]
    dist.to.low.dens <- di[phat.aug <= lambda]
    index.of.low.dens <- (1:n)[phat.aug <= lambda]
    argmin <- which.min(dist.to.low.dens)
    j <- index.of.low.dens[argmin]
    r <- dist.to.low.dens[argmin]
    radii[i] <- r
    if (interpolation) {
      delta.phat <- phat[i] - phat[j]
      if (abs(delta.phat) > 0) {
        a <- (phat[i] - lambda) / delta.phat
        radii[i] <- a * r
      }
    }
  }
  return(radii)
}

compute.weights <- function(X, R, lambdas, phat) {
  n <- nrow(X)
  D2 <- gsl.interpoint.distance.matrix(X, X)
  D2[D2 < 0] <- 0
  D <- sqrt(D2)
  W <- matrix(0, nrow = n, ncol = n)
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      ##browser()
      r <- R[i,] + R[j,]
      cand.lambdas <- lambdas[(lambdas < min(phat[i], phat[j]))]
      cand.r <- r[(lambdas < min(phat[i], phat[j]))]
      W[i, j] <- max(cand.lambdas[cand.r >= D[i, j]])
      W[j, i] <- W[i, j]
    }
    W[i, i] <- max(lambdas[lambdas < phat[i]])
    cat(" ", i)
  }
  W[n,n] <- max(lambdas[lambdas < phat[n]])
  return(W)
}


##-----------------------------------------------------------------

walther.similarity.matrix <- function(X, phat, n.grid = nrow(X),
                                      interpolation = F) {
  compute.radii <- function(lambda) {
    radii  <- rep(0, n)
    for (i in 1:n) {
      di <- D.aug[i,]
      dist.to.low.dens <- di[phat.aug <= lambda]
      index.of.low.dens <- (1:n)[phat.aug <= lambda]
      argmin <- which.min(dist.to.low.dens)
      j <- index.of.low.dens[argmin]
      r <- dist.to.low.dens[argmin]
      radii[i] <- r
      if (interpolation) {
        delta.phat <- phat[i] - phat[j]
        if (abs(delta.phat) > 0) {
          a <- (phat[i] - lambda) / delta.phat
          radii[i] <- a * r
        }
      }
    }
    return(radii)
  }
##-----------------------------------------------------------------
  compute.weights <- function() {
    W <- matrix(0, nrow = n, ncol = n)
    for (i in 1:(n-1)) {
      for (j in (i+1):n) {
        ##browser()
        r <- R[i,] + R[j,]
        cand.lambdas <- lambdas[(lambdas < min(phat[i], phat[j]))]
        cand.r <- r[(lambdas < min(phat[i], phat[j]))]
        W[i, j] <- max(cand.lambdas[cand.r >= D[i, j]])
        W[j, i] <- W[i, j]
      }
      W[i, i] <- max(lambdas[lambdas < phat[i]])
      cat(" ", i)
    }
    W[n,n] <- max(lambdas[lambdas < phat[n]])
    cat(" ", n)
    return(W)
  }
  ## -----------------------------------------------------------------
  n <- nrow(X)
  D2 <- gsl.interpoint.distance.matrix(X, X)
  D2[D2 < 0] <- 0
  D <- sqrt(D2)
  D.aug <- cbind(D, rep(max(D), n))
  phat.aug <- c(phat, 0)
  grid.ind <- floor(seq(1, n, length = n.grid))
  lambdas <- sort(phat)[grid.ind]
  if (n.grid == n) lambdas <- sort(phat)
  R <- matrix(0, nrow = n, ncol = n.grid)
  cat("\nComputing radii ")
  for (i in 1:n.grid) {
    R[,i] <- compute.radii(lambdas[i])
    cat(" ", i)
  }
  lambdas <- c(0, lambdas)
  n.lambdas <- n.grid + 1
  R <- cbind(max(D), R)
  cat("\nComputing weights ")
  W <- compute.weights()
  sim <- W[upper.tri(W, diag = T)]
  return(sim)
}

##=================================================================
## Experiment with nearest neigbor density estimates

nearest.neighbor.density.estimate <- function(X.obs, X.eval, k) {
  n.eval <- nrow(X.eval)
  n.obs <- nrow(X.obs)
  dens <- rep(0, n.eval)
  D2 <- gsl.interpoint.distance.matrix(X.eval, X.obs)
  D2[D2 < 0] <- 0
  D <- sqrt(D2)
  min.pos <- min(D[D > 0])
  D[D == 0] <- 0.1 * min.pos
  for (i in 1:n.eval) {
    dens[i] <- 1.0 / sort(D[i,])[k]
  }
  return(dens)
}

