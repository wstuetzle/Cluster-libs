## Make data for experiments to evaluate clustering procedures
## ===========================================================
##
## 3-17-2014: Added functions to fit Gaussian mixture models
##
## Write data using "dump" in files with extension ".R". Each file has
## a matrix "X" with the features and a vector "group.id" with the group ids.
## Don't write standardized and sphered versions. Those can be generated
## on the fly

## List of data sets
##------------------
## "three-gaussians-in-2d.R"
## "simplex-3.R"
## "simplex-5.R"
## "simplex-10.R"
## "simplex-3-small.R"
## "simplex-5-small.R"
## "simplex-10-small.R"
## "simplex-3-1999.R"
## "simplex-7-1999.R"
## "bullseye-2.R"
## "bullseye-2-2000.R"
## "bullseye-10-2000.R"
## "bullseye-2-smooth-015.R"
## "olive-5-lda-2d.R"
## "olive-5-lda-5d.R"
## "olive-5-lda-10d.R"
## "olive.R"
## "crab.R"
## "zip-train-256.R"
## "zip-train-64.R"
## "zip-train-16.R"
## "zip-train-500-256.R"
## "zip-train-500-64.R"
## "zip-train-500-16.R"
## "zip-train-500-256-lda-9d.R"
## "tdt-pca-log-ldf-500.R"
## "tdt-pca-log-ldf-16.R"
## "tdt-pca-log-ldf-64.R"
## "tdt-pca-log-ldf-256.R"
## "all7-215-12558.R"
## "all7-215-500.R"
## "all7-215-1000.R"
## "all7-215-2000.R"
## "all7-215-1000-rn.R"
## "all7-215-1000-rn-10PC.R"
## "banana-clump.R"
## "circle-clump.R"
## "baudry-ex41.R"    From Baudry, Raftery et al, "Combining mixture
## "baudry-ex42.R"    components for clustering"
## "baudry-ex43.R"
## "baudry-ex441.R"
## "baudry-ex442.R"
## "baudry-gvhd-minus.R"
## "baudry-gvhd-plus.R"
## "nice-univariate-trimodal.R"
## "swiss-banknote.R"
## "iris.R"

library("MASS")
library("stats")


## source("z:\\Werner\\Clustering\\Progs\\clusfun-11-22-04.R")
## data.directory <- "Z:\\Werner\\Clustering\\Data"
## output.directory <- "Z:\\Werner\\Clustering\\Data\\Test-collection"

dir.sep <- "/"
clustering.dir <- "/Users/wxs/Dropbox/Clustering"
data.directory <- paste(clustering.dir, "Data", sep = dir.sep)
output.directory <- paste(clustering.dir, "Data", "Test-collection", sep = dir.sep)

##=================================================================
## A simple test data set, just to check out the program: 2d, 3
## Gaussian clusters, 20 observations each.

means <- matrix(0, nrow = 3, ncol = 2)
means[1,] <- c(0, 0)
means[2,] <- c(0, sqrt(2))
means[3,] <- c(sqrt(2), 0)

clust.sizes <- c(20, 20, 20)
clust.sds <- rep(0.25, 3)

n <- sum(clust.sizes)
m <- 2
X <- matrix(0, nrow=n, ncol=m)
group.id <- rep(0,n)

low <- 1
for(i in 1:length(clust.sizes)) {
  high <- low + clust.sizes[i] - 1
  X[low:high,] <- matrix(rep(means[i,], clust.sizes[i]), ncol=m, byrow=T) +
    matrix(rnorm(clust.sizes[i]*m, sd=clust.sds[i]), ncol=m)
  group.id[low:high] <- i
  low <- low + clust.sizes[i]
}

output.filename <- "three-gaussians-in-2d.R"
output.pathname <- paste(output.directory, output.filename, sep = "\\")
creation.date <- date()
if (!file.exists(output.pathname)) {
  dump(c("X", "group.id", "creation.date"), output.pathname)
  cat("\nMade three-gaussians-in-2d")
}

##=================================================================
## Make Gaussian and smeared uniform groups for mixture model figure
## showing that # of mixture components is not necessarily =
## number of groups

nobs <- 800
p <- 2
covmat <- matrix(c(1, 0.8, 0.8, 1), nrow = 2, byrow = T)
root.covmat <- chol(covmat)
Z <- matrix(rnorm(nobs*p), nrow = nobs)
unimodal <- Z %*% root.covmat
offset.vector <- c(1, 0)
offset.amount <- 3
parallel.gaussian <- unimodal
for (i in 1:(nobs/2)) {parallel.gaussian[i,] <- unimodal[i,] + offset.vector * offset.amount}
## plot(parallel.gaussian, pch = 20, xlab = "", ylab = "", xaxt = "n", yaxt = "n")
X <- parallel.gausssian
group.id <- rep(2, nobs)
group.id[1:nobs/2] <- 1
output.filename <- "mixmod-demo-gaussian.R"
output.pathname <- paste(output.directory, output.filename, sep = dir.sep)
creation.date <- date()
if (!file.exists(output.pathname)) {
  dump(c("X", "group.id", "creation.date"), output.pathname)
  cat("\nMade mixmod-demo-gaussian")
}

## Make parallel uniformn clusters %------------------------------------
nobs <- 800
p <- 2
kernel.sd <- 0.1
covmat <- matrix(c(1, 0.8, 0.8, 1), nrow = 2, byrow = T)
root.covmat <- chol(covmat)
## t(root.covmat) %*% root.covmat
U <- matrix(runif(nobs*p), nrow = nobs)
G <- matrix(rnorm(nobs*p, sd = kernel.sd), nrow = nobs)
Z <- U + G
unimodal <- Z %*% root.covmat
offset.vector <- c(1, 0)
offset.amount <- 1.1
parallel.uniform <- unimodal
for (i in 1:(nobs/2)) {parallel.uniform[i,] <- unimodal[i,] +
                           offset.vector * offset.amount}
## plot(parallel.uniform, pch = 20, xlab = "", ylab = "", xaxt = "n", yaxt = "n")
X <- parallel.uniform
group.id <- rep(1, nobs)
group.id[1:nobs/2] <- 2
output.filename <- "mixmod-demo-uniform.R"
output.pathname <- paste(output.directory, output.filename, sep = dir.sep)
creation.date <- date()
if (!file.exists(output.pathname)) {
  dump(c("X", "group.id", "creation.date"), output.pathname)
  cat("\nMade mixmod-demo-uniform")
}


##=================================================================
## Make 2 elongated clusters in 2d

make.elongated.groups <- function(n1, n2, cov, offset) {
    n <- n1 + n2
    covmat <- matrix(c(1, cov, cov, 1), nrow = 2, byrow = T)
    root.covmat <- chol(covmat)
    Z <- matrix(rnorm(2 *n), nrow = n)
    X <- Z %*% root.covmat
    group.id <- rep(1, n)
    for (i in (n1 + 1):n) {
        X[i, ] <- X[i,] + c(offset, 0)
        group.id[i] <- 2
    }
    return(list(X = X, group.id = group.id))
}

bla <- make.elongated.groups(100, 100, 0.9, 0)
X <- bla$X
group.id <- rep(1, 200)
output.filename <- "one-elongated-group-in-2d.R"
output.pathname <- paste(output.directory, output.filename, sep = dir.sep)
creation.date <- date()
if (!file.exists(output.pathname)) {
  dump(c("X", "group.id", "creation.date"), output.pathname)
  cat("\nMade one-elongated-group-in-2d")
}

bla <- make.elongated.groups(100, 100, 0.9, 2.5)
X <- bla$X
## plot(X, pch = 20)
group.id <- rep(1, 200)
output.filename <- "two-elongated-groups-in-2d.R"
output.pathname <- paste(output.directory, output.filename, sep = dir.sep)
creation.date <- date()
if (!file.exists(output.pathname)) {
  dump(c("X", "group.id", "creation.date"), output.pathname)
  cat("\nMade two-elongated-groups-in-2d")
}

##=================================================================
## Make 2 intersecting Gaussian clusters in 2d

make.crossing.groups <- function(n, cov) {
  covmat1 <- matrix(c(1, cov, cov, 1), nrow = 2, byrow = T)
  covmat2 <- matrix(c(1, -cov, -cov, 1), nrow = 2, byrow = T)
  root.covmat1 <- chol(covmat1)
  root.covmat2 <- chol(covmat2)
  Z <- matrix(rnorm(2*n), nrow = n)
  group1 <- Z %*% root.covmat1
  group2 <- Z %*% root.covmat2
  cross <- rbind(group1, group2)
  return(cross)
}

n <- 400
cov <- 0.95
X <- make.crossing.groups(n, cov)
## plot(X, pch = 20)
group.id <- rep(2, 2*n)
group.id[1:n] <- 1
output.filename <- "two-crossing-groups-in-2d.R"
output.pathname <- paste(output.directory, output.filename, sep = dir.sep)
creation.date <- date()
if (!file.exists(output.pathname)) {
  dump(c("X", "group.id", "creation.date"), output.pathname)
  cat("\nMade two-elongated-groups-in-2d")
}

##=================================================================
## The simplex data

make.simplex.data <- function(clust.sizes, clust.sds) {
  p <- length(clust.sizes)
  means <- matrix(0, nrow=p, ncol=p)
  means[row(means)==col(means)] <- 1
  n <- sum(clust.sizes)
  X <- matrix(0, nrow=n, ncol=p)
  group.id <- rep(0,n)
  low <- 1
  for(i in 1:p) {
    high <- low + clust.sizes[i] - 1
    X[low:high,] <- matrix(rep(means[i,], clust.sizes[i]), ncol=p, byrow=T) +
                    matrix(rnorm(clust.sizes[i]*p, sd=clust.sds[i]), ncol=p)
    group.id[low:high] <- i
    low <- low + clust.sizes[i]
  }
  return(list(X = X, group.id = group.id))
}

## Three dimensions
p <- 3
clust.sizes <- seq(from = 50, by = 10, length = p)
clust.sds <- rep(0.25, p)
msd.out <- make.simplex.data(clust.sizes, clust.sds)
X <- msd.out$X
group.id <- msd.out$group.id
output.filename <- "simplex-3.R"
output.pathname <- paste(output.directory, output.filename, sep = "\\")
creation.date <- date()
if (!file.exists(output.pathname)) {
  dump(c("X", "group.id", "creation.date"), output.pathname)
  cat("\nMade simplex-3")
}

## Five dimensions
p <- 5
clust.sizes <- seq(from = 50, by = 10, length = p)
clust.sds <- rep(0.25, p)
msd.out <- make.simplex.data(clust.sizes, clust.sds)
X <- msd.out$X
group.id <- msd.out$group.id
output.filename <- "simplex-5.R"
output.pathname <- paste(output.directory, output.filename, sep = "\\")
creation.date <- date()
if (!file.exists(output.pathname)) {
  dump(c("X", "group.id", "creation.date"), output.pathname)
  cat("\nMade simplex-5")
}

## Ten dimensions
p <- 10
clust.sizes <- seq(from = 50, by = 10, length = p)
clust.sds <- rep(0.25, p)
msd.out <- make.simplex.data(clust.sizes, clust.sds)
X <- msd.out$X
group.id <- msd.out$group.id
output.filename <- "simplex-10.R"
output.pathname <- paste(output.directory, output.filename, sep = "\\")
creation.date <- date()
if (!file.exists(output.pathname)) {
  dump(c("X", "group.id", "creation.date"), output.pathname)
  cat("\nMade simplex-10")
}

##-----------------------------------------------------------------
## Smaller version of the simplex data

## Three dimensions
p <- 3
clust.sizes <- seq(from = 30, by = 5, length = p)
clust.sds <- rep(0.25, p)
msd.out <- make.simplex.data(clust.sizes, clust.sds)
X <- msd.out$X
group.id <- msd.out$group.id
output.filename <- "simplex-3-small.R"
output.pathname <- paste(output.directory, output.filename, sep = "\\")
creation.date <- date()
if (!file.exists(output.pathname)) {
  dump(c("X", "group.id", "creation.date"), output.pathname)
  cat("\nMade simplex-3-small")
}

## Five dimensions
p <- 5
clust.sizes <- seq(from = 30, by = 5, length = p)
clust.sds <- rep(0.25, p)
msd.out <- make.simplex.data(clust.sizes, clust.sds)
X <- msd.out$X
group.id <- msd.out$group.id
output.filename <- "simplex-5-small.R"
output.pathname <- paste(output.directory, output.filename, sep = "\\")
creation.date <- date()
if (!file.exists(output.pathname)) {
  dump(c("X", "group.id", "creation.date"), output.pathname)
  cat("\nMade simplex-5-small")
}

## Ten dimensions
p <- 10
clust.sizes <- seq(from = 30, by = 5, length = p)
clust.sds <- rep(0.25, p)
msd.out <- make.simplex.data(clust.sizes, clust.sds)
X <- msd.out$X
group.id <- msd.out$group.id
output.filename <- "simplex-10-small.R"
output.pathname <- paste(output.directory, output.filename, sep = "\\")
creation.date <- date()
if (!file.exists(output.pathname)) {
  dump(c("X", "group.id", "creation.date"), output.pathname)
  cat("\nMade simplex-10-small")
}

##-----------------------------------------------------------------
## There are also 1999 versions of the simplex data in 3d and 7d
## Convert those to R format so we can see if we can reproduce
## old simulation results

data.restore("z:\\Werner\\Clustering\\Data\\simplex-data-12-21-99.dmp", print = T)
X <- X.simplex.3
group.id <- labels.simplex.3
output.filename <- "simplex-3-1999.R"
output.pathname <- paste(output.directory, output.filename, sep = "\\")
creation.date <- date()
if (!file.exists(output.pathname)) {
  dump(c("X", "group.id", "creation.date"), output.pathname)
  cat("\nMade simplex-3-1999")
}

X <- X.simplex.7
group.id <- labels.simplex.7
output.filename <- "simplex-7-1999.R"
output.pathname <- paste(output.directory, output.filename, sep = "\\")
creation.date <- date()
if (!file.exists(output.pathname)) {
  dump(c("X", "group.id", "creation.date"), output.pathname)
  cat("\nMade simplex-7-1999")
}


##=================================================================
## Bullseye data

rbull <- function(n, eye.radius, inner.ring.radius, outer.ring.radius) {
  square.area <- (2 * outer.ring.radius)^2
  eye.area <- pi * eye.radius^2
  ring.area <- pi * (outer.ring.radius^2 - inner.ring.radius^2)
  efficiency <- (eye.area + ring.area) / square.area
  n.square <- trunc(2 * n / efficiency)
  X.square <- matrix(outer.ring.radius * (2 * runif(2*n.square) - 1), ncol = 2)
  r <- sqrt(X.square[,1]^2 + X.square[,2]^2)
  in.bullseye <- (r < eye.radius) | 
  ((r > inner.ring.radius) & (r < outer.ring.radius))
  X.bullseye <- X.square[in.bullseye,]
  r.bullseye <- r[in.bullseye]
  labels <- rep(2, length(r.bullseye))
  labels[r.bullseye < eye.radius] <- 1
  return(list(X = X.bullseye[1:n,], group.id = labels[1:n]))	               
}

eye.radius <- 1
inner.ring.radius <- 1.5
outer.ring.radius <- 2
n.bull <- 500
rbull.out <- rbull(n.bull, eye.radius, inner.ring.radius, outer.ring.radius)
X <- rbull.out$X
group.id <- rbull.out$group.id
plot(X, pch = 20)
output.filename <- "bullseye-2.R"
output.pathname <- paste(output.directory, output.filename, sep = "\\")
creation.date <- date()
if (!file.exists(output.pathname)) {
  dump(c("X", "group.id", "creation.date"), output.pathname)
  cat("\nMade bullseye-2")
}

## There also is a 2000 version of the bullseye data
data.restore("z:\\Werner\\Clustering\\Data\\bullseye-data-6-16-00.dmp", print = T)
X <- X.bull.2
group.id <- labels.bull
output.filename <- "bullseye-2-2000.R"
output.pathname <- paste(output.directory, output.filename, sep = "\\")
creation.date <- date()
if (!file.exists(output.pathname)) {
  dump(c("X", "group.id", "creation.date"), output.pathname)
  cat("\nMade bullseye-2-2000")
}

## 10d data has 8 noise dimensions with sd = 0.1
X <- X.bull.10
group.id <- labels.bull
output.filename <- "bullseye-10-2000.R"
output.pathname <- paste(output.directory, output.filename, sep = "\\")
creation.date <- date()
if (!file.exists(output.pathname)) {
  dump(c("X", "group.id", "creation.date"), output.pathname)
  cat("\nMade bullseye-10-2000")
}

##=================================================================
## Smooth bullseye (added 11-23-06)
## Convolute with Gaussian kernel

eye.radius <- 1
inner.ring.radius <- 1.5
outer.ring.radius <- 2
kernel.sd <- 0.15
n.bull <- 500
rbull.out <- rbull(n.bull, eye.radius, inner.ring.radius, outer.ring.radius)
Y <- kernel.sd * matrix(rnorm(2*n.bull), ncol = 2)
X <- rbull.out$X + Y
group.id <- rbull.out$group.id
plot(X, pch = 20)
output.filename <- "bullseye-2-smooth-015.R"
output.pathname <- paste(output.directory, output.filename, sep = "\\")
creation.date <- date()
if (!file.exists(output.pathname)) {
  dump(c("X", "group.id", "creation.date"), output.pathname)
  cat("\nMade bullseye-2-smooth")
}


##=================================================================
## 1. Areas 5-9 of olive oil data, projected on first two
## discriminant coordinates
## 1-31-08: It seems the noise for the noisy data is much too large

olive.pathname <- paste(data.directory, "olive.dat", sep = "\\")
olive.dat <- matrix(scan(olive.pathname), ncol = 10, byrow=T)
region <- olive.dat[,1]
area <- olive.dat[,2]
reg23 <- (area >= 5)
X.olive.reg23 <- olive.dat[reg23, 3:10]
area.reg23 <- area[reg23]

lda.out <- lda(X.olive.reg23, area.reg23)
pairs(lda.out, dimen = 2)
X.proj <- X.olive.reg23 %*% lda.out$scaling
plot(X.proj[,1:2], type = "n")
text(X.proj[,1:2], labels = area.reg23, cex = 0.7)

X <- X.proj[,1:2]
group.id <- area.reg23

output.filename <- "olive-5-lda-2d.R"
output.pathname <- paste(output.directory, output.filename, sep = "\\")
creation.date <- date()
if (!file.exists(output.pathname)) {
  dump(c("X", "group.id", "creation.date"), output.pathname)
  cat("\nMade olive-5-lda-2d")
}

## Add 3 Gaussian noise variables
input.filename <-  "olive-5-lda-2d.R"
input.pathname <- paste(data.directory, input.filename, sep = "\\")
source(input.pathname)
noise.sd <- sqrt(sum(diag(var(X)))/2)
noise.dim <- 3
noise.mat <- matrix(rnorm(nrow(X) * noise.dim, sd = noise.sd),
                    ncol = noise.dim)
X <- cbind(X, noise.mat)
output.filename <- "olive-5-lda-5d.R"
output.pathname <- paste(output.directory, output.filename, sep = "\\")
creation.date <- date()
if (!file.exists(output.pathname)) {
  dump(c("X", "group.id", "creation.date"), output.pathname)
  cat("\nMade olive-5-lda-5d")
}

## Add 8 Gaussian noise variables
input.filename <-  "olive-5-lda-2d.R"
input.pathname <- paste(data.directory, input.filename, sep = "\\")
source(input.pathname)
noise.sd <- sqrt(sum(diag(var(X)))/2)
noise.dim <- 8
noise.mat <- matrix(rnorm(nrow(X) * noise.dim, sd = noise.sd),
                    ncol = noise.dim)
X <- cbind(X, noise.mat)
output.filename <- "olive-5-lda-10d.R"
output.pathname <- paste(output.directory, output.filename, sep = "\\")
creation.date <- date()
if (!file.exists(output.pathname)) {
  dump(c("X", "group.id", "creation.date"), output.pathname)
  cat("\nMade olive-5-lda-10d")
}


##=================================================================
## The olive oil data

olive.pathname <- paste(data.directory, "olive.dat", sep = "\\")
olive.dat <- matrix(scan(olive.pathname), ncol = 10, byrow=T)
group.id <- olive.dat[,2]
X <- olive.dat[,3:10]
output.filename <- "olive.R"
output.pathname <- paste(output.directory, output.filename, sep = "\\")
creation.date <- date()
if (!file.exists(output.pathname)) {
  dump(c("X", "group.id", "creation.date"), output.pathname)
  cat("\nMade olive")
}


##=================================================================
## The crab data

crab.pathname <- paste(data.directory, "crab-matrix.dat", sep = "\\")
crab.all <- matrix(scan(crab.pathname), ncol = 7, byrow = T)
species <- crab.all[,1]
gender <- crab.all[,2]
group.id <- 2 * (species-1) + gender -1
X <- crab.all[, 3:7]
output.filename <- "crab.R"
output.pathname <- paste(output.directory, output.filename, sep = "\\")
creation.date <- date()
if (!file.exists(output.pathname)) {
  dump(c("X", "group.id", "creation.date"), output.pathname)
  cat("\nMade crab")
}


##=================================================================
## The ZIP handwritten digit training data

input.pathname <- paste(data.directory, "zip-train.dat", sep = dir.sep)
zip.train <- matrix(scan(input.pathname), ncol = 257, byrow = T)
## group.id <- zip.train[,1]
group.id <- zip.train[,1] + 1  ## Some code my assume that group.id = 1...k
X <- zip.train[,2:257]
output.pathname <- paste(output.directory, "zip-train-256.R", sep = dir.sep)
creation.date <- date()
if (!file.exists(output.pathname)) {
  dump(c("X", "group.id", "creation.date"), output.pathname)
  cat("\nMade zip-train-256")
}
## source(output.pathname, echo = F)   ## Takes about 5min to source
XX <- X

## Projection on 16 largest PCs
zip.train.eigen <- eigen(var(XX))
X <- XX %*% zip.train.eigen$vectors[,1:16]
output.filename <- "zip-train-16.R"
output.pathname <- paste(output.directory, output.filename, sep = "\\")
creation.date <- date()
if (!file.exists(output.pathname)) {
  dump(c("X", "group.id", "creation.date"), output.pathname)
  cat("\nMade zip-train-16")
}

## Projection on 64 largest PCs
X <- XX %*% zip.train.eigen$vectors[,1:64]
output.filename <- "zip-train-64.R"
output.pathname <- paste(output.directory, output.filename, sep = "\\")
creation.date <- date()
if (!file.exists(output.pathname)) {
  dump(c("X", "group.id", "creation.date"), output.pathname)
  cat("\nMade zip-train-64")
}

##=================================================================
## Make smaller version of the ZIP data with 500 obs

input.pathname <- paste(data.directory, "zip-train.dat", sep = "\\")
zip.train <- matrix(scan(input.pathname), ncol = 257, byrow = T)
group.id.all <- zip.train[,1]
X.all <- zip.train[,2:257]
n <- length(group.id)
m <- 500
sub <- sample(1:n, m)
X <- X.all[sub,]
group.id <- group.id.all[sub]
output.filename <- "zip-train-500-256.R"
output.pathname <- paste(output.directory, output.filename, sep = "\\")
creation.date <- date()
if (!file.exists(output.pathname)) {
  dump(c("X", "group.id", "creation.date"), output.pathname)
  cat("\nMade zip-train-500-256")
}
XX <- X

## Projection on 16 largest PCs
zip.train.eigen <- eigen(var(XX))
X <- XX %*% zip.train.eigen$vectors[,1:16]
output.filename <- "zip-train-500-16.R"
output.pathname <- paste(output.directory, output.filename, sep = "\\")
creation.date <- date()
if (!file.exists(output.pathname)) {
  dump(c("X", "group.id", "creation.date"), output.pathname)
  cat("\nMade zip-train-500-16")
}

## Projection on 64 largest PCs
X <- XX %*% zip.train.eigen$vectors[,1:64]
output.filename <- "zip-train-500-64.R"
output.pathname <- paste(output.directory, output.filename, sep = "\\")
creation.date <- date()
if (!file.exists(output.pathname)) {
  dump(c("X", "group.id", "creation.date"), output.pathname)
  cat("\nMade zip-train-500-64")
}

##=================================================================
## LDA for ZIP data

input.filename <- "zip-train-500-256.R"
input.pathname <- paste(output.directory, input.filename, sep = "\\")
output.filename <- "zip-train-500-256-lda-9d.R"
output.pathname <- paste(output.directory, output.filename, sep = "\\")

source(input.pathname)
lda.out <- lda(X, group.id)
pairs(lda.out, dimen = 2)
X.proj <- X %*% lda.out$scaling

plot(X.proj[,1:2])
plot(X.proj[,1:2], type = "n")
text(X.proj[,1:2], labels = group.id, cex = 0.7)

X <- X.proj[,1:9]
plot(X[,1:2])
creation.date <- date()
if (!file.exists(output.pathname)) {
  dump(c("X", "group.id", "creation.date"), output.pathname)
  cat("\nMade zip-train-500-256-lda-9d")
}


##=================================================================
## TDT data

label.data <- read.table("Z:\\Werner\\Nsa\\Docclust\\tdt-labels.dat",
			 header = F, col.names = c("doc.id", "topic"))

# Document TDT008481 appears with two different topics, 22 and 3. 

doc.id <- label.data$doc.id
topic <- label.data$topic
length(doc.id)
length(topic)
topic[doc.id=="TDT008481"]
n <- length(doc.id)

# Arbitrarily delete one entry.

index.of.duplicate <- (1:n)[doc.id=="TDT008481"][2]
topic <- topic[- index.of.duplicate]
doc.id <- doc.id[- index.of.duplicate]
length(doc.id)
length(topic)
length(unique(topic))
n <- length(doc.id)
group.id <- topic



## Now it checks out.
## Read in projections of documents onto first 500 PCs
## I asume that XXt.log.ldf.pca contains the projections of the
## documents on the 500 largest PCs.
## data.restore("d:\\Werner\\Nsa\\Docclust\\pca.log.ldf", print=T)
## Does not work in R
## 4-10-06: Make new file pca-log-ldf.dat which has only the matrix entries.
## I assume the matrix is written in column-major order

X <- matrix(scan("Z:\\Werner\\NSA\\Docclust\\pca-log-ldf.dat"), nrow = n, ncol = 500)
output.filename <- "tdt-pca-log-ldf-500.R"
output.pathname <- paste(output.directory, output.filename, sep = "\\")
creation.date <- date()
if (!file.exists(output.pathname)) {
  dump(c("X", "group.id", "creation.date"), output.pathname)
  cat("\nMade tdt-pca-log-ldf-500.R")
}
XX <- X

## 16 Dimensions
X <- XX[,1:16]
output.filename <- "tdt-pca-log-ldf-16.R"
output.pathname <- paste(output.directory, output.filename, sep = "\\")
creation.date <- date()
if (!file.exists(output.pathname)) {
  dump(c("X", "group.id", "creation.date"), output.pathname)
  cat("\nMade tdt-pca-log-ldf-16.R")
}

## 64 Dimensions
X <- XX[,1:64]
output.filename <- "tdt-pca-log-ldf-64.R"
output.pathname <- paste(output.directory, output.filename, sep = "\\")
creation.date <- date()
if (!file.exists(output.pathname)) {
  dump(c("X", "group.id", "creation.date"), output.pathname)
  cat("\nMade tdt-pca-log-ldf-64.R")
}

## 256 Dimensions
X <- XX[,1:256]
output.filename <- "tdt-pca-log-ldf-256.R"
output.pathname <- paste(output.directory, output.filename, sep = "\\")
creation.date <- date()
if (!file.exists(output.pathname)) {
  dump(c("X", "group.id", "creation.date"), output.pathname)
  cat("\nMade tdt-pca-log-ldf-256.R")
}


##=================================================================
## The ALL Leukemia data

all.folder <- "ALL-cell-cycle-2-27-07"
all.X.filename <- "all7_cancer.txt"
all.classes.filename <- "all7_cancer_classes.txt"
all.X.pathname <- paste(data.directory, all.folder, all.X.filename,
                        sep = "\\")
all.classes.pathname <- paste(data.directory, all.folder, all.classes.filename,
                        sep = "\\")
bla <- scan(all.X.pathname)
n <- bla[1]
n.genes <- bla[2]
XX <- matrix(bla[3:length(bla)], nrow = n, byrow = T)
classes <- factor(scan(all.classes.pathname, what = ""))
attr(classes, "levels")
group.id <- as.integer(classes)
group.id
table(group.id)
attr(classes, "levels")

## The classes are
## [1] "BCR-ABL"     "E2A-PBX1"    "Hyperdip>50" "MLL"         "OTHERS"     
## [6] "T-ALL"       "TEL-AML1"   


##-----------------------------------------------------------------
## Dump all the data
X <- XX
data.filename <- "all7-215-12558.R"
output.pathname <- paste(output.directory, data.filename, sep = "\\")
output.pathname
creation.date <- date()
if (!file.exists(output.pathname)) {
  dump(c("X", "group.id", "creation.date"), output.pathname)
  cat("\nMade all7-215-12558.R")
}

##-----------------------------------------------------------------
## Make data sets containing only the genes with highest variance
sd.genes <- apply(XX, 2, sd)
XX.sd.ordered <- XX[,rev(order(sd.genes))]

X <- XX.sd.ordered[,1:500]
data.filename <- "all7-215-500.R"
output.pathname <- paste(output.directory, data.filename, sep = "\\")
output.pathname
creation.date <- date()
if (!file.exists(output.pathname)) {
  dump(c("X", "group.id", "creation.date"), output.pathname)
  cat("\nMade all7-215-500.R")
}

X <- XX.sd.ordered[,1:1000]
data.filename <- "all7-215-1000.R"
output.pathname <- paste(output.directory, data.filename, sep = "\\")
output.pathname
creation.date <- date()
if (!file.exists(output.pathname)) {
  dump(c("X", "group.id", "creation.date"), output.pathname)
  cat("\nMade all7-215-1000.R")
}

X <- XX.sd.ordered[,1:2000]
data.filename <- "all7-215-2000.R"
output.pathname <- paste(output.directory, data.filename, sep = "\\")
output.pathname
creation.date <- date()
if (!file.exists(output.pathname)) {
  dump(c("X", "group.id", "creation.date"), output.pathname)
  cat("\nMade all7-215-2000.R")
}

## Pick 1000 genes with highest variance, normalize rows
XXX <- XX.sd.ordered[,1:1000]
row.means <- apply(XXX, 1, mean)
XXX.row.centered <- sweep(XXX, 1, row.means)
row.sds <- apply(XXX, 1, sd)
XXX.row.normalized <- sweep(XXX.row.centered, 1, row.sds, FUN = "/")
X <- XXX.row.normalized
data.filename <- "all7-215-1000-rn.R"
output.pathname <- paste(output.directory, data.filename, sep = "\\")
output.pathname
creation.date <- date()
if (!file.exists(output.pathname)) {
  dump(c("X", "group.id", "creation.date"), output.pathname)
  cat("\nMade all7-215-1000-rn.R")
}

xrn.eigen <- eigen(var(XXX.row.normalized))
evals <- xrn.eigen$values[1:(n-1)]
dim(xrn.eigen$vectors)
evects <- xrn.eigen$vectors[,1:(n-1)]
plot(cumsum(evals)/sum(evals), type = "l")
m <- 10
frac.explained[m]
X <- XXX.row.normalized %*% evects[,1:m]
data.filename <- "all7-215-1000-rn-10PC.R"
output.pathname <- paste(output.directory, data.filename, sep = "\\")
output.pathname
creation.date <- date()
if (!file.exists(output.pathname)) {
  dump(c("X", "group.id", "creation.date"), output.pathname)
  cat("\nMade all7-215-1000-rn-10PC.R")
}

##=================================================================
## The banana-clump data (11-23-07)

## Here is the code to make the banana-clump data
## n.banana <- 100
## sd.banana <- 0.04
## banana.curve <- define.curve(loop = F)
## banana.obs <- create.points.around.curve(banana.curve, n.banana, sd.banana)

## n.clump <- 100
## sd.clump <- 0.1
## mu.clump <- c(0.6, 0.4)
## shift.mat <- matrix(mu.clump, nrow = n.clump, ncol = 2, byrow = T)
## clump.obs <- sd.clump * matrix(rnorm(2*n.clump), ncol = 2) + shift.mat
## X <- rbind(banana.obs, clump.obs)
## plot(X, pch = 20)
## banana.clump.obs.pathname <- paste(figure.directory, "banana-clump-obs.R", sep = "\\")
##dump("X", banana.clump.obs.pathname)


data.filename <- "banana-clump.R"
output.pathname <- paste(output.directory, data.filename, sep = "\\")
source(paste(data.directory, "\\banana-clump-obs.R", sep = ""))
group.id <- c(rep(1, 100), rep(2, 100))
creation.date <- date()
if (!file.exists(output.pathname)) {
  dump(c("X", "group.id", "creation.date"), output.pathname)
  cat("\nMade banana-clump.R")
}

##=================================================================
## A more extreme version of the banana-clump data
## A spherical Gaussian cluster at the origin and a spherical arc

generate.points.on.arc <- function(n, r, phi.min, phi.max) {
  ## phi <- runif(n, phi.min, phi.max)
  phi <- seq(from = phi.min, to = phi.max, length.out = n)
  X <- matrix(0, nrow = n, ncol = 2)
  X[,1] <- r * sin(phi)
  X[,2] <- r * cos(phi)
  return(X)
}

n.arc <- 300
r <- 5
phi.min <- 0
phi.max <- 2 * pi
noise.sd <- 0.5
n.clump <- 100

arc <- generate.points.on.arc(n.arc, r, phi.min, phi.max)
noise.sd <- matrix(rnorm(2 * n.arc, sd = noise.sd), nrow = n.arc, ncol = 2)
clump <- matrix(rnorm(2 * n.clump), nrow = n.clump, ncol = 2)
X <- rbind(clump, arc + noise.sd)
plot(X, pch = 19)
group.id <- c(rep(1, n.clump), rep(2, n.arc))

data.filename <- "circle-clump.R"
output.pathname <- paste(output.directory, data.filename, sep = "\\")
creation.date <- date()
if (!file.exists(output.pathname)) {
  dump(c("X", "group.id", "creation.date"), output.pathname)
  cat("\nMade circle-clump.R")
}

##=================================================================
## Nice univariate trimodal data to illustrate CWC

data.filename <- "nice-univariate-trimodal.R"
output.pathname <- paste(output.directory, data.filename, sep = "\\")

if (!file.exists(output.pathname)){
  n.comp <- 3
  n <- 100
  mu <- c(0, 4, 10)
  sd <- c(1, 1, 1)
  mix <- c(0.4, 0.4, 0.2)
  x <- NULL
  for (i in 1:n.comp) {
    x <- c(x, rnorm(n * mix[i], mean = mu[i], sd = sd[i]))
  }
  X <- as.matrix(sort(x), ncol = 1)
  X <- sphere(X)
  group.id <- c(rep(1, n * mix[1]), rep(2, n * mix[2]), rep(3, n * mix[3]))
  creation.date <- date()
  dump(c("X", "group.id", "creation.date"), output.pathname)
  cat("\nMade nice-univariate-trimodal,R")
}

##=================================================================
## Swiss banknote data

data.filename <- "swiss-banknote.R"
output.pathname <- paste(output.directory, data.filename, sep = dir.sep)

library("alr3")
data(banknote)
X <- as.matrix(banknote)[,1:6]
group.id <- banknote$Y + 1
creation.date <- date()

dump(c("X", "group.id", "creation.date"), output.pathname)
cat("\nMade swiss-banknote.R")

##=================================================================
## Iris data

data(iris)
n <- nrow(iris)
X <- as.matrix(iris[,1:4])
species <- iris[,5]
group.id <- rep(0, n)
group.id[species == "setosa"] <- 1
group.id[species == "versicolor"] <- 2
group.id[species == "virginica"] <- 3
creation.date <- date()
output.pathname <- paste(output.directory, "iris.R", sep = dir.sep)
dump(c("X", "group.id", "creation.date"), file = output.pathname)
cat("\nMade iris.R")


##=================================================================
##=================================================================
## Functions to fit Gaussian mixture model
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
  for (i in 1:ncomp) {
    mpi <- mixture.par[[i]]
    ni <- comp.n[i]
    Xi <- mvrnorm(ni, mpi$mu, mpi$Sigma, empirical = T)
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

