make.knn.density.estimate <- function(X.train, k = 1) {
  density.estimate <- function(X.eval) {
     phat <- nn.density.estimate(X.train, X.eval, k)
     return(phat)
   }
  return(density.estimate)
}

nn.density.estimate <- function(X.train, X.eval, k = 1) {
  if (is.vector(X.train)) {
    X.train <- matrix(X.train, ncol = 1)
    X.eval <- matrix(X.eval, ncol = 1)
  }
  m <- ncol(X.train)
  nobs <- nrow(X.train)
  neval <- nrow(X.eval)
  dens <- rep(0, neval)
  dist <- gsl.interpoint.distance.matrix(X.eval, X.train)
  dist[dist < 0] <- 0
  dist <- sqrt(dist)
  for (i in 1:neval) {
    r <- sort(dist[i,])[k]
    if (r == 0) dens[i] <- Inf else dens[i] <- 1 / r
  }
  return(dens)
}

## n <- 100
## X <- sort(c(rnorm(n), rnorm(n) + 5))
## hist(X)

## k <- 100
## phat <- nn.density.estimate(X, X, k)
## plot(X, phat, type = "l")

## X.obs <- c(1, 2, 3)
## X.eval <- c(4, 5)
## cwc.nn.density.estimate(X.obs, X.eval, k = 2)
