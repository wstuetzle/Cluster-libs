#################################################################
## Functions for Gaussian mixture models

## Fit Gaussian mixture model

fit.gmm <- function(X, group.id) {
  ngroup <- max(group.id)
  n <- nrow(X)
  mixture.par <- NULL
  if (is.vector(X)) X <- matrix(X, ncol = 1)
  m <- ncol(X)
  for (i in 1:ngroup){
    Xi <- matrix(X[group.id == i,], ncol = m)
    mui <- apply(Xi, 2, mean)
    Sigmai <- var(Xi)
    ni <- sum(group.id == i)
    mixture.par[[i]] <- list(p = ni/n, mu = mui, Sigma = Sigmai)
  }
  return(mixture.par)
}

## nut.filename <- "nut-1.R"
## nut.pathname <- paste(cwc.dir, "JSM-2013", "Data", nut.filename, sep = dir.sep)
## source(nut.pathname, echo = T)

## mixture.par <- fit.gmm(X, group.id)

##-----------------------------------------------------------------
## Generate sample from Gaussian mixture model
## Multinomial = T => sample sizes from mixture components are
## multinomial
## Multinomial = F => sample sizes from mixture components are the
## expected values

rgmm <- function(n, mixture.par, multinomial = T) {
  ncomp <- length(mixture.par)
  m <- length(mixture.par[[1]]$mu)
  prob <- rep(0, ncomp)
  for (i in 1:ncomp) prob[i] <- mixture.par[[i]]$p
  if (multinomial) {
    comp.n <- rmultinom(1, n, prob)
  }
  if (!multinomial) comp.n <- round(n * prob, 0)
  n <- sum(comp.n)
  X <- matrix(0, nrow = n, ncol = m)
  group.id <- rep(0, n)
  istart <- 1
  ## browser()
  for (i in 1:ncomp) {
    ## cat(" ", i)
    mpi <- mixture.par[[i]]
    ni <- comp.n[i]
    ## Xi <- mvrnorm(ni, mpi$mu, mpi$Sigma, empirical = F)
    Xi <- rmvnorm(ni, mpi$mu, mpi$Sigma)
    X[istart:(istart + ni - 1),] <- Xi
    group.id[istart:(istart + ni - 1)] <- i
    istart <- istart + ni
  }
  return(list(X = X, group.id = group.id))
}

## X <- rgmm(1000, mixture.par)$X
## hist(X, nclass = 20)

##-----------------------------------------------------------------
## Evaluate mixture density

dgmm <- function(X, mixture.par) {
  ncomp <- length(mixture.par)
  if (is.vector(X)) X <- matrix(X, nrow = 1)
  n <- nrow(X)
  dens <- rep(0, n)
  for (i in 1:ncomp) {
    mpi <- mixture.par[[i]]
    prob <- mpi$p
    mu <- mpi$mu
    Sigma <- mpi$Sigma
    dens <- dens + prob * dmvnorm(X, mu, Sigma, log=FALSE)
  }
  return(dens)
}

## dens <- dgmm(X, mixture.par)

##-----------------------------------------------------------------

make.dgmm <- function(mixture.par) {
  dens <- function(X) {
    return(dgmm(X, mixture.par))
  }
  return(dens)
}

##-----------------------------------------------------------------

make.dgmm.gradient <- function(mixture.par ){
  ncomp <- length(mixture.par)
  m <- length(mixture.par[[1]]$mu)
  Sigma.inv <- NULL
  for (i in 1:ncomp) {
    Sigma <- mixture.par[[i]]$Sigma
    Sigma.eigen <- eigen(Sigma)
    U <- Sigma.eigen$vectors
    lambda <- Sigma.eigen$values
    for (j in 1:m) U[,j] <- U[,j] / sqrt(lambda[j])
    Sigma.inv[[i]] <- U %*% t(U)
  }
  dgmm.gradient <- function(x) {
    gradient <- rep(0, m)
    for (i in 1:ncomp) {
      mpi <- mixture.par[[i]]
      prob <- mpi$p
      mu <- mpi$mu
      Sigma <- mpi$Sigma
      dens <- prob * dmvnorm(x, mu, Sigma, log=FALSE)
      gradient <- gradient - dens * Sigma.inv[[i]] %*% (x - mu)
    }
    return(as.vector(gradient))
  }
  return(dgmm.gradient)
}


## dens <- make.dgmm(mixture.par)
## x <- rep(0, m)

## eps <- 1.e-8
## fd.gradient <- rep(0, m)
## for (i in 1:m) {
##   y <- x
##   y[i] <- y[i] + eps
##   fd.gradient[i] <- (dens(y) - dens(x)) / eps
## }

## fd.gradient

## dgmm.gradient <- make.dgmm.gradient(mixture.par)
## dgmm.gradient(x)
