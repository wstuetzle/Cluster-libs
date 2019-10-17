## Try Clustering with Confidence (CWC)
## 5-19-09, 7-13-2010
##=================================================================

## There will be two versions, differing in how they decide when to
## split. The "vertex" version splits if (i) the threshold graph
## Gplus(lambda) has two connected components; (ii) in each connected
## component there is a vertex for which the lower confidence bound
## for the density is > lambda (a vertex which is also in
## Gminus(lambda)).
## The "edge" version (to be used with the 1nn
## density estimate) splits if (i) (as above) and if (ii) in each
## connected component there is an edge which is also Gminus(lambda).

## The vertex sets of the connected components of Gplus(lambda) are
## the same as the vertex sets of the connected components of
## Mastplus(lambda), where Mastplus is the maximal spanning tree of
## Gplus. However, the edge set of Gplus(lambda) is obviously
## different from the edge set of Mastplus(lambda). So the result of
## the "edge" version depends on whether we apply test (ii) to
## Gplus(lambda) or to Mastplus(lambda). It's easier to do the latter,
## so that's what we will do. We can compute Mastplus using the GSL
## code.

##=================================================================

## OS <- "Windows"
## ## OS <- "Linux"

## if (OS == "Windows") {
##   cwc.dir <- "Z:\\Werner\\Clustering\\CWC"
##   dir.sep <- "\\"
## }
## if (OS == "Linux") {
##   cwc.dir <- "/homes/wxs/Werner/Clustering/CWC"
##   dir.sep <- "/"
## }

OS <- "OSX"
echo <- T

if (OS == "OSX") {
  cluster.dir <- "/Users/wxs/Dropbox/Clustering"
  dir.sep = "/"
}

libs.dir <- paste(cluster.dir, "A-Git-repos", "Cluster-libs", sep = "/")
source(paste(libs.dir, "gsl-functions-6-8-2018.R", sep = dir.sep), echo = echo)
source(paste(libs.dir, "cwc-functions-8-28-2010-v2017.R", sep = dir.sep), echo = echo)
source(paste(libs.dir, "gaussian-mixture-model-functions.R", sep = dir.sep), echo = echo)

cwc.dir <- paste(cluster.dir, "CWC-via-confidence-bands-before-Germany-2011",
                 sep = dir.sep)
code.dir <- paste(cwc.dir, "Code", sep = dir.sep)
data.dir <- paste(cwc.dir, "Data", sep = dir.sep)

library(MASS)
library(mvtnorm)
quartz()

##=================================================================
##=================================================================
##=================================================================
## Chapter 1: Make some data
##
## Make univariate Gaussian mixture
##

data.filename <- "nice-univariate-trimodal.R"
data.pathname <- paste(data.dir, data.filename, sep = dir.sep)

if(file.exists(data.pathname)) source(data.pathname, echo = echo)

if (!file.exists(data.pathname)){
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
  ## dump(c("X", "group.id"), data.pathname)
}

##-----------------------------------------------------------------
## Make uniforms in various dimensions
##
n1 <- 50
ms <- c(1, 2, 5, 10)

for (m in ms) {
  data.filename <- paste("uniform-n", n1 * m, "-m", m, ".R", sep = "")
  data.pathname <- paste(data.dir, data.filename, sep = dir.sep)
  if (!file.exists(data.pathname)) {
    X <- matrix(runif(n1 * m^2), ncol = m)
    if (m == 1) X <- matrix(sort(runif(n1 * m^2)), ncol = m)
    X <- sphere(X)
    group.id <- rep(1, n1 * m)
    ## dump(c("X", "group.id"), data.pathname)
  }
}

## n <- 560
## m <- 7

n <- 400
m <- 2

n <- 572
m <- 8

data.filename <- paste("uniform-n", n, "-m", m, ".R", sep = "")
data.pathname <- paste(data.dir, data.filename, sep = dir.sep)
file.exists(data.pathname)
if (!file.exists(data.pathname)) {
  X <- matrix(runif(n * m), ncol = m)
  if (m == 1) X <- matrix(sort(runif(n * m)), ncol = m)
  X <- sphere(X)
  group.id <- rep(1, n)
  dump(c("X", "group.id"), data.pathname)
}

##-----------------------------------------------------------------
## Make standard Gaussians in various dimensions
##
n1 <- 50
ms <- c(1, 2, 5, 10)

for (m in ms) {
  X <- matrix(rnorm(n1 * m^2), ncol = m)
  X <- sphere(X)
  group.id <- rep(1, n1 * m)
  data.filename <- paste("gaussian-n", n1 * m, "-m", m, ".R", sep = "")
  data.pathname <- paste(data.dir, data.filename, sep = dir.sep)
  data.pathname <- paste(data.dir, data.filename, sep = dir.sep)
  cat(file.exists(data.pathname), " ")
  ## dump(c("X", "group.id"), data.pathname)
}

##-----------------------------------------------------------------
## Make fake olive oil data from gaussian mixture fitted to the olive
## data.

## Version (1): Use area labels.

library(MASS)
data.filename <- "olive.R"
data.pathname <- paste(data.dir, data.filename, sep = dir.sep)
source(data.pathname, echo = echo)
Y <- NULL
for (i in 1:max(group.id)) {
  ind <- (group.id == i)
  mu <- apply(X[ind, ], 2, mean)
  Sigma <- var(X[ind,])
  Z <- mvrnorm(sum(ind), mu, Sigma)
  Y <- rbind(Y, Z)
}
X <- Y
creation.date <- date()
data.filename <- "olive-fake-v1.R"
data.pathname <- paste(data.dir, data.filename, sep = dir.sep)
file.exists(data.pathname)
## dump(c("X", "group.id", "creation.date"), file = data.pathname)
dim(X)


##=================================================================
##=================================================================
##=================================================================
## Chapter 2: Play with univariate Gaussian mixture
##
data.name <- "nice-univariate-trimodal"
data.filename <- "nice-univariate-trimodal.R"
data.pathname <- paste(data.dir, data.filename, sep = dir.sep)
source(data.pathname, echo = echo)
n <- nrow(X)
n
kmax <- 1

hist(X, breaks = 20, col = "green")

 mu <- c(0, 4, 10)
  sd <- c(1, 1, 1)
  mix <- c(0.4, 0.4, 0.2)

## nut.density <- function(x) {
##   n <- length(x)
##   p <- rep(0, n)
##   for (i in 1:3) {
##     p <- p + mix[i] * dnorm(x, mu[i], sd[i])
##   }
##   return(p)
## }


##-----------------------------------------------------------------
## Try kernel density estimate with LS CV

trial.h <- c(0.01, 0.05, 0.10, 0.15, 0.20, 0.30, 0.40, 0.60, 0.80, 1.00, 1.50, 2.00)
cvs.out <- cv.search(X, trial.par = trial.h)
cvs.out
h <- cvs.out$opt.smopar
h

## h <- 0.05  ## For illustration, make estimate more noisy
density <- make.gaussian.kernel.density.estimate(X, h)

## Compute vertex and edge weights for test graph. Vertices equally spaced
## in range of data

n.vertices <- 100
nv <- n.vertices
ne <- nv - 1
tg.vertices <- as.matrix(seq(min(X), max(X), length = n.vertices), ncol = 1)
tg.edges <- cbind((1:(n.vertices - 1)), (2:n.vertices))

## Alternatively, use observations as vertices
tg.vertices <- X   ## already sorted
tg.edges <- cbind((1:(n - 1)), (2:n))

ngrid <- 10
cwcvw.out <- cwc.ve.weights(tg.vertices, tg.edges, density, ngrid)
vw <- cwcvw.out$vertex.weights
ew <- cwcvw.out$edge.weights

## Plot vertex and edge weights
plot(tg.vertices, vw, pch = 20, col = "red")
for (i in 1:nrow(tg.edges)) {
  e <- tg.edges[i,]
  v1 <- tg.vertices[e[1]]
  v2 <- tg.vertices[e[2]]
  lines(c(v1, v2), rep(ew[i], 2), lwd = 3)
}
lines(tg.vertices, density(tg.vertices))
title("Vertex and edge weights")


##-----------------------------------------------------------------
## Try computing resample vertex and edge weights with CV for both
## regular bootstrap and half-sampling

##-----------------------------------------------------------------
## Generate or read resample weights

resample.type <- "half.sampling"
## resample.type <- "bootstrap"
nres <- 100  ## number of resamples
ngrid <- 10
bandwidth <- "cv"
data.name <- "nice-univariate-trimodal"

rw.filename <- resample.weights.filename(data.name, resample.type, nres,
                                         bandwidth, kmax, ngrid)
rw.filename <- paste(rw.filename, ".R", sep = "")
rw.pathname <- paste(data.dir, rw.filename, sep = dir.sep)
rw.pathname

file.exists(rw.pathname)
if (file.exists(rw.pathname)) source(rw.pathname, echo = echo)
if (!file.exists(rw.pathname)) {
  resample.kernel.ve.weights <-
    cwc.resample.kernel.ve.weights(X, kmax = kmax, nres = nres,
                                   ngrid = ngrid, bandwidth = bandwidth,
                                   resample.type = resample.type,
                                   debug.filename = stdout())
  dump("resample.kernel.ve.weights", file = rw.pathname)
}
names(resample.kernel.ve.weights)
 ## [1] "tg.vertices"             "tg.edges"               
 ## [3] "vertex.weights"          "edge.weights"           
 ## [5] "orig.cvs"                "resample.vertex.weights"
 ## [7] "resample.edge.weights"   "resamples"              
 ## [9] "resample.h"              "resample.cvs"      
cwcbvw.hs <- resample.kernel.ve.weights
cwcbvw.hs$tg.vertices
cwcbvw.hs$vertex.weights
dim(cwcbvw.hs$resample.vertex.weights)
dim(cwcbvw.hs$resample.edge.weights)
length(cwcbvw.hs$resample.h)

hist(cwcbvw.hs$resample.h, nclass = 20, col = "green", main = "")
title("Bandwidths for resamples")

vw.hs <- cwcbvw.hs$resample.vertex.weights
dim(vw.hs)
ew.hs <- cwcbvw.hs$resample.edge.weights


## Bootstrap
## Note: LS CV does not work on bootstrap samples becuase of ties.
## We use the bandwidth estimated from the original sample for all
## the bootstrap samples

## resample.type <- "half.sampling"
resample.type <- "bootstrap"
nres <- 100  ## number of resamples
ngrid <- 10
bandwidth <- "cv"

rw.filename <- resample.weights.filename(data.name, resample.type, nres,
                                         bandwidth, kmax, ngrid)
rw.filename <- paste(rw.filename, ".R", sep = "")
rw.pathname <- paste(data.dir, rw.filename, sep = dir.sep)
rw.pathname

file.exists(rw.pathname)
if (file.exists(rw.pathname)) source(rw.pathname, echo = echo)
if (!file.exists(rw.pathname)) {
  resample.kernel.ve.weights <-
    cwc.resample.kernel.ve.weights(X, kmax = kmax, nres = nres,
                                   ngrid = ngrid, bandwidth = bandwidth,
                                   resample.type = resample.type,
                                   debug.filename = stdout())
  dump("resample.kernel.ve.weights", file = rw.pathname)
}
names(resample.kernel.ve.weights)
cwcbvw.bt <- resample.kernel.ve.weights

vw.bt <- cwcbvw.bt$resample.vertex.weights
ew.bt <- cwcbvw.bt$resample.edge.weights


## Stick bt and hs weights into arrays

nv <- nrow(vw.bt)
rvw <- array(0, c(nv, nres, 2),
            dimnames = list(NULL, NULL, c("hs", "bt")))
rvw[,,"hs"] <- vw.hs
rvw[,,"bt"] <- vw.bt

ne <- nrow(ew.bt)
rew <- array(0, c(ne, nres, 2),
            dimnames = list(NULL, NULL, c("hs", "bt")))
rew[,,"hs"] <- ew.hs
rew[,,"bt"] <- ew.bt


##-----------------------------------------------------------------
## Plot bootstrap and halfsampling vertex weights

rt <- "hs"
## rt = "bt"
ylim <- c(0, max(rvw[,,rt]))
for (i in 1:nres) {
  plot(tg.vertices, rvw[,i,rt], ylim = ylim, type = "l", lwd = 3, col = "red")
  for (j in 1:i) {
    lines(tg.vertices, rvw[,j,rt])
  }
  lines(tg.vertices, rvw[,i,rt], lwd = 3, col = "red")
  bla <- readline("\nHit return to see next curve or any character to exit")
  if (bla != "") break
}

ylim <- range(c(vw.hs, vw.bt))
xlim = range(tg.vertices)

plot(0, 0, xlim = xlim, ylim = ylim, type = "n")
for (i in 1:nres) {
  points(tg.vertices, rvw[, i, "bt"], pch = ".")
  ## points(tg.vertices, rvw[, i, "hs"], pch = ".", col = "red")
}

##-----------------------------------------------------------------
## Compute Gaussian non-simultaneous confidence intervals for edge weights
## Check for normality of bvw

## resample.type <- "bt"
resample.type <- "hs"
while (TRUE) {
  i <- sample(1:nv, 1)
  qqnorm(rvw[i,,resample.type], main = "", pch = 20)
  title(paste("Normal Q-Q plot for resamples, vertex ", i))
  bla <- readline("\nHit return to see next plot or any character to exit")
  if (bla != "") break
}

while (TRUE) {
  i <- sample(1:ne, 1)
  qqnorm(rew[i,,resample.type], main = "", pch = 20)
  title(paste("Normal Q-Q plot for resamples, edge ", i))
  bla <- readline("\nHit return to see next plot or any character to exit")
  if (bla != "") break
}

rvw.means <- apply(rvw, c(1, 3), mean)
dim(rvw.means)
rvw.sds <- apply(rvw, c(1, 3), sd)

plot(rvw.means[,"hs"], rvw.means[,"bt"], pch = 20)
abline(0, 1, col = "red", lwd = 3)

plot(rvw.sds[,"hs"], rvw.sds[,"bt"], pch = 20)
abline(0, 1, col = "red", lwd = 3)
##
## Not a fair comparison, because the bootstrap uses the bandwidth of the
## original sample for all the resamples. LS CV fails for bootstrap resamples
## because of ties

plot(tg.vertices, rvw.sds[,"hs"], ylim = c(0, max(rvw.sds)), type = "l")
lines(tg.vertices, rvw.sds[,"bt"], col = "red")

ylim <- range(rvw.sds / sqrt(rvw.means))
plot(tg.vertices, rvw.sds[,"hs"] / sqrt(rvw.means[,"hs"]), type = "l",
     ylim = ylim)
lines(tg.vertices, rvw.sds[,"bt"] / sqrt(rvw.means[,"bt"]), col = "red")

## The difference between resampling with and without replacement needs
## some more looking into

## Gaussian confidence intervals could be centered around mean of bootstrap
## distributions or about the weights for the original sample
## Look at the difference

plot(tg.vertices, vw, type = "l")
lines(tg.vertices, rvw.means[,"hs"], col = "red")
lines(tg.vertices, rvw.means[,"bt"], col = "blue")

## The mean of the bootstrap distributions is smoother (at least in this
## example)

##-----------------------------------------------------------------
## Compute and plot empirical and Gaussian intervals
##
## 4-12-2018
## For this example with nrep = 100 it's impossible to get empirical
## confidence bounds with simultaneous coverage probability = 0.8 or bigger.
## Non-simultaeous going from the 2nd smallest to the 2nd largest weights
## only have simultaneous coverage prob 0.74

target.cov.prob <- 0.1
target.cov.prob <- 0.5
target.cov.prob <- 0.8
target.cov.prob <- 0.85
target.cov.prob <- 0.9
target.cov.prob <- 0.95
simultaneous <- T
## simultaneous <- F
rt <- "hs"

if (rt == "hs") resample.kernel.ve.weights <- cwcbvw.hs
if (rt == "bt") resample.kernel.ve.weights <- cwcbvw.bt

band.type <- "empirical"
empirical.confidence.bounds <-
  cwc.confidence.bounds(resample.kernel.ve.weights, target.cov.prob,
                        band.type = band.type, simultaneous = simultaneous)
cov.parameter <- empirical.confidence.bounds$cov.parameter
cov.parameter
# cov.parameter  <- 1
cwc.simul.cov.prob(resample.kernel.ve.weights, cov.parameter,
                   band.type = band.type)

band.type <- "gaussian"
gaussian.confidence.bounds <-
  cwc.confidence.bounds(resample.kernel.ve.weights, target.cov.prob,
                        band.type = band.type, simultaneous = simultaneous)
cov.parameter <- gaussian.confidence.bounds$cov.parameter
cov.parameter
cwc.simul.cov.prob(resample.kernel.ve.weights, cov.parameter,
                   band.type = band.type)
2 * (1-pnorm(cov.parameter)) 

band.type <- "poisson"
poisson.confidence.bounds <-
  cwc.confidence.bounds(resample.kernel.ve.weights, target.cov.prob,
                        band.type = band.type, simultaneous = simultaneous)
cov.parameter <- poisson.confidence.bounds$cov.parameter
cov.parameter
cwc.simul.cov.prob(resample.kernel.ve.weights, cov.parameter,
                   band.type = band.type)
## Something is not right here for poisson bounds
## 10-25-2018  Looks ok


## Plot empirical and Gaussian intervals

upper.empirical <- empirical.confidence.bounds$upper.vw
lower.empirical <- empirical.confidence.bounds$lower.vw
upper.gaussian <- gaussian.confidence.bounds$upper.vw
lower.gaussian <- gaussian.confidence.bounds$lower.vw
upper.poisson <- poisson.confidence.bounds$upper.vw
lower.poisson <- poisson.confidence.bounds$lower.vw


ylim <- c(min(rvw[,,rt]), max(rvw[,,rt]))
plot(tg.vertices, upper.empirical, ylim = ylim,
     type = "l", col = "red")
lines(tg.vertices, lower.empirical, col = "red")
lines(tg.vertices, vw, type = "l")
lines(tg.vertices, lower.gaussian, col = "blue")
lines(tg.vertices, upper.gaussian, col = "blue")
lines(tg.vertices, lower.poisson, col = "green")
lines(tg.vertices, upper.poisson, col = "green")
## Add resample curves
for (i in 1:nres) points(tg.vertices, rvw[,i,rt], pch = ".")

## What we really want to plot is the upper edges weights and the lower
## vertex.weights

## plot(tg.vertices, vw, pch = 19, type = "n")
## for (i in 1:nrow(tg.edges)) {
##   e <- tg.edges[i,]
##   v1 <- tg.vertices[e[1]]
##   v2 <- tg.vertices[e[2]]
##   uew <- ew[i] + z.cov * bew.sds[i]
##   lines(c(v1, v2), rep(uew, 2), lwd = 3)
## }
## ## points(tg.vertices, vw - z.cov * bvw.sds, pch = 20, col = "red")
## lines(tg.vertices, vw - z.cov * bvw.sds, col = "red")

## -----------------------------------------------------------------
## Compute simultaneous coverage prob of vertex intervals around
## vertex weights of original sample


## Non-simultaneous coverage pro:   0.99     0.95
## Simultaneous Gaussian:           0.89     0.53
## Simultaneous empirical:          0.87     0.51

## If we wanted simultaneous coverage prob = 0.8, Bonferroni would
## suggest non-simultaneous coverage prob 1 - 0.2 / 100 = 0.998
## So Bonferroni is pessimistic

## We are wasting a lot on areas in whch we are not really interested.
## We are much too conservative. We should focus on minima and peaks.

##-----------------------------------------------------------------
## Try intervals of the form
## [phat - gamma.minus * sqrt(phat), phat + gamma.plus * sqrt(phat)]
##
## Let's find gamma.minus, gamma.plus by fitting lower empirical vertex
## weight confidence bounds and upper empirical edge weight confidence
## bounds, respectively

## If y_i are the confidence bounds and phat_i are the weights, then
## the LS estimate for gamma is
##
## gamma = [sum(y_i sqrt(phat_i)) - sum(phat_i^{3/2)] / sum(phat_i)

confidence.level = 0.99
rt <- "hs"
cwc.compare.confidence.bounds(vw, ew, rvw[,,rt], rew[,,rt],
                              confidence.level = confidence.level)




##-----------------------------------------------------------------
## Try simple GSL clustering

trial.h <- c(0.01, 0.05, 0.10, 0.15, 0.20, 0.30, 0.40, 0.60, 0.80, 1.00,
             1.50, 2.00)
cvs.out <- cv.search(X, trial.par = trial.h)
cvs.out
h <- cvs.out$opt.smopar
h
## h <- 0.05  ## For illustration, make estimate more noisy
density <- make.gaussian.kernel.density.estimate(X, h)
plot(X, density(X), type = "l")

## mdsfun <- gsl.make.kernel.mdsfun(X, h)
mdsfun <- gsl.make.kernel.mdsfun(X, h, kmax = 1)
gsl.cluster.out <- gsl.cluster(X, mdsfun, pruning.criterion = "size",
                        pruning.threshold = 0, assign.fluff = F)
gsl.runt.size(gsl.cluster.out)[1:10]
gsl.runt.rise(gsl.cluster.out)[1:10]
gsl.cluster.out <- gsl.cluster(X, mdsfun, pruning.criterion = "size",
                        pruning.threshold = 18, assign.fluff = T)

cluster.id <- gsl.observation.labels(gsl.cluster.out)
unique(cluster.id)
in.core <- gsl.observation.in.core(gsl.cluster.out)

##-----------------------------------------------------------------
## Plot density

nboot <- 100
pboot <- rep(0, n)
for (i in 1:nboot) {
  for (j in 1:n) pboot[j] <- bootstrap.kernel.mdsfun(j, j)[i]
  ## plot(X[,1], density(X), type = "l")
  plot(X[,1], density(X), pch = 20, ylim = c(0, max(pboot, density(X))))
  ## lines(X[,1], pboot, col = "red")
  points(X[,1], pboot, col = "red", pch = 20)
  title(paste("Bootstrap sample ", i))
  if (!is.null(opt.level)) abline(h = opt.level)
  ##browser()
  lines(c(X[split.edge[1], 1], X[split.edge[2], 1]),
        rep(split.edge.weights[i], 2), lwd = 7, col = "green")
  readline("\nHit return to see next curve ")
}

plot(X[,1], density(X), type = "l")
abline(h = opt.level)
for (i in supporting.bs) {
  for (j in 1:n) pboot[j] <- bootstrap.kernel.mdsfun(j, j)[i]
  lines(X[,1], pboot, col = "red")
  readline("\nHit return to see next curve ")
}

##-----------------------------------------------------------------
## Try the dip test

library(diptest)
x <- X[,1]
dip.out <- dip(x, full.result = T, debug = T)
dip.out

n <- length(x)
ux <- unique(x)
nux <- length(ux)
ecdf <- rep(0, n)
for (i in 1:nux) ecdf[x == ux[i]] <- sum(x <= ux[i])/n
ecdf
plot(x, ecdf, type = "n")
for (i in 1:n) {
  lines(c(x[i], x[i+1]), rep(ecdf[i], 2))
}
lo <- dip.out$lo.hi[1] ## 44
hi <- dip.out$lo.hi[2] ## 74

plot(x, ecdf, pch = ".")
lines(c(x[lo], x[hi]), c(ecdf[lo], ecdf[hi]), col = "red")

gcm <- dip.out$gcm   ## greatest convex minorant
lcm <- dip.out$lcm   ## least concave majorant

m <- sum(gcm > 0)
for (i in 1:(m-1)) {
  index.int <- c(gcm[i], gcm[i+1])
  lines(x[index.int], ecdf[index.int], col = "red")
}
lines(c(x[74], x[100]), c(ecdf[74], ecdf[100]), col = "blue")
        

bla



##=================================================================
##=================================================================
##=================================================================
##=================================================================
## Chapter 3: Try CWC

trial.h <- c(0.01, 0.05, 0.10, 0.15, 0.20, 0.30, 0.40, 0.60, 0.80, 1.00,
             1.50, 2.00)

data.name <- "nice-univariate-trimodal"
kmax <- 1

data.name <- "banana-clump-cv"
kmax <- 20

data.name <- "uniform-n50-m1"
kmax <- 1

data.name <- "uniform-n100-m2"
kmax <- 20

data.name <- "uniform-n250-m5"
kmax <- 20

data.name <- "uniform-n500-m10"
kmax <- 20

data.name <- "gaussian-n50-m1"
kmax <- 1

data.name <- "gaussian-n100-m2"
kmax <- 20

data.name <- "gaussian-n250-m5"
kmax <- 20

data.name <- "gaussian-n500-m10"
kmax <- 20

data.name <- "bullseye"
kmax <- 20

data.name <- "olive"
kmax <- 20

data.name <- "olive-fake-v1"
kmax <- 20

data.name <- "simplex-7-1999"
kmax <- 20

data.filename <- paste(data.name, ".R", sep = "")
data.pathname <- paste(data.dir, data.filename, sep = dir.sep)
source(data.pathname,echo = T)
n <- nrow(X)
X <- sphere(X)
diag(var(X))


##-----------------------------------------------------------------
## Generate or read resample weights

resample.type <- "half.sampling"
## resample.type <- "bootstrap"
nres <- 100  ## number of resamples
kmax <- 20
## kmax <- 1
ngrid <- 10
bandwidth <- "cv"
## bandwidth <- "cv.mean"

rw.filename <- paste(resample.weights.filename(data.name, resample.type, nres,
                                         bandwidth, kmax, ngrid),
                     ".R", sep = "")
rw.pathname <- paste(data.dir, rw.filename, sep = dir.sep)
rw.pathname

## if (file.exists(rw.pathname)) source(rw.pathname, echo = echo)
## Don't source - It's way too slow. Look below!

if (!file.exists(rw.pathname)) {
  resample.kernel.ve.weights <-
    cwc.resample.kernel.ve.weights(X, kmax = kmax, nres = nres,
                                   ngrid = ngrid, bandwidth = bandwidth,
                                   resample.type = resample.type,
                                   debug.filename = stdout())
  dump("resample.kernel.ve.weights", file = rw.pathname)
}
names(resample.kernel.ve.weights)
 ## [1] "tg.vertices"             "tg.edges"               
 ## [3] "vertex.weights"          "edge.weights"           
 ## [5] "orig.cvs"                "resample.vertex.weights"
 ## [7] "resample.edge.weights"   "resamples"              
 ## [9] "resample.h"              "resample.cvs"


#################################################################
## 4-17-2018
##
## Sourcing large files is verrrrry slow. So I sourced the old .R file
## and then used "save" to write a new .sav file

## I did this for
"banana-clump-cv-rw-hs-100-cv-20-10"
"olive-rw-hs-100-cv-20-10"

rw.save.filename <- paste(resample.weights.filename(data.name, resample.type, nres,
                                         bandwidth, kmax, ngrid),
                     ".sav", sep = "")
rw.save.pathname <- paste(data.dir, rw.save.filename, sep = dir.sep)

## save(resample.kernel.ve.weights, file = rw.save.pathname)

load (rw.save.pathname)
names(resample.kernel.ve.weights)

## For olive oil data read old version of bootstrap.weights
## rw.filename <- "olive-bootstrap-weights-k20.R"
## rw.pathname <- paste(data.dir, rw.filename, sep = dir.sep)
## source(rw.pathname, echo = T)
## names(bootstrap.weights)
## resample.kernel.ve.weights <-
##   list(tg.vertices = bootstrap.weights$tg.vertices,
##        tg.edges = bootstrap.weights$tg.edges,
##        resample.vertex.weights = bootstrap.weights$bootstrap.vertex.weights,
##        resample.edge.weights = bootstrap.weights$bootstrap.edge.weights,
##        resamples = bootstrap.weights$bootstrap.samples,
##        resample.h <- bootstrap.weights$bootstrap.h,
##        resample.cvs <- bootstrap.weights$cvs.out
##        )

##-----------------------------------------------------------------
## Look at vertex and edge weights
vw <- resample.kernel.ve.weights$vertex.weights
ew <- resample.kernel.ve.weights$edge.weights
rvw <- resample.kernel.ve.weights$resample.vertex.weights
rew <- resample.kernel.ve.weights$resample.edge.weights

dim(rvw)
for (i in 1:nrow(rvw)) {
  ## qqnorm(log(rvw[i,]), pch = 20)
  ## hist(log(rvw[i,]), col = "green")
  qqnorm(rvw[i,], pch = 20)
  ## hist(rvw[i,], col = "green")
  ## abline(v = vw[i], lwd = 3, col = "red")
  print(sort(rvw[i,]))
  bla <- readline("\nHit return to see next curve or any character to exit")
  if (bla != "") break
}


##-----------------------------------------------------------------
## Calculate confidence bounds

target.cov.prob <- 0.1
target.cov.prob <- 0.5
target.cov.prob <- 0.8
target.cov.prob <- 0.85
target.cov.prob <- 0.9
target.cov.prob <- 0.95
simultaneous = T

band.type <- "empirical"
empirical.confidence.bounds <-
  cwc.confidence.bounds(resample.kernel.ve.weights, target.cov.prob,
                        band.type = band.type, simultaneous = simultaneous)
cov.parameter <- empirical.confidence.bounds$cov.parameter
cov.parameter 
cov.parameter <- 1
cwc.simul.cov.prob(resample.kernel.ve.weights, cov.parameter,
                   band.type = band.type)

band.type <- "gaussian"
gaussian.confidence.bounds <-
  cwc.confidence.bounds(resample.kernel.ve.weights, target.cov.prob,
                        band.type = band.type, simultaneous = simultaneous)
cov.parameter <- gaussian.confidence.bounds$cov.parameter
cov.parameter
cwc.simul.cov.prob(resample.kernel.ve.weights, cov.parameter,
                   band.type = band.type)

band.type <- "poisson"
poisson.confidence.bounds <-
  cwc.confidence.bounds(resample.kernel.ve.weights, target.cov.prob,
                        band.type = band.type, simultaneous = simultaneous)
cov.parameter <- poisson.confidence.bounds$cov.parameter
cov.parameter
cwc.simul.cov.prob(resample.kernel.ve.weights, cov.parameter,
                   band.type = band.type)

##-----------------------------------------------------------------
## Grapically compare bounds

plot(poisson.confidence.bounds$upper.vw, gaussian.confidence.bounds$upper.vw,
     pch = 20)
abline(0, 1, lwd = 3, col = "red")

plot(poisson.confidence.bounds$upper.vw, empirical.confidence.bounds$upper.vw,
     pch = 20)
abline(0, 1, lwd = 3, col = "red")

plot(poisson.confidence.bounds$upper.ew, gaussian.confidence.bounds$upper.ew,
     pch = 20)
abline(0, 1, lwd = 3, col = "red")

plot(poisson.confidence.bounds$lower.vw, gaussian.confidence.bounds$lower.vw,
     pch = 20)
abline(0, 1, lwd = 3, col = "red")

plot(poisson.confidence.bounds$lower.vw, empirical.confidence.bounds$lower.vw,
     pch = 20)
abline(0, 1, lwd = 3, col = "red")

plot((poisson.confidence.bounds$upper.vw - poisson.confidence.bounds$lower.vw),
     (gaussian.confidence.bounds$upper.vw - gaussian.confidence.bounds$lower.vw),
     pch = 20)
abline(0, 1, lwd = 3, col = "red")

plot((poisson.confidence.bounds$upper.vw - poisson.confidence.bounds$lower.vw),
     (empirical.confidence.bounds$upper.vw - empirical.confidence.bounds$lower.vw),
     pch = 20)
abline(0, 1, lwd = 3, col = "red")

## Gaussian bounds seem to be much wider than poisson bounds ????
##
## There seems to be a fundamental problem for olive oil data. Even with
## nout = 1, empirical intervals have simultaneous coverage probability = 0.
## For Gaussian and poisson intervals, lower bounds are all negative.
## So edge and vertex weights must be very non-normal.

##-----------------------------------------------------------------
## Do some exploratory analysis of edge weights for olive oil data.

## > names(resample.kernel.ve.weights)
##  [1] "tg.vertices"             "tg.edges"               
##  [3] "vertex.weights"          "edge.weights"           
##  [5] "orig.cvs"                "resample.vertex.weights"
##  [7] "resample.edge.weights"   "resamples"              
##  [9] "resample.h"              "resample.cvs"

rvw <- resample.kernel.ve.weights$resample.vertex.weights
dim(rvw)
vw <- resample.kernel.ve.weights$vertex.weights
n <- nrow(rvw)

mean.rvw <- apply(rvw, 1, mean)
plot(vw, mean.rvw, pch = 20)
hist(log(vw), col = "green")
sum(vw > 0.15)
## Highly correlated for large values of vw
## The lower bound for the estimated density at an observation is
## 1/n (value of kernel at 0) which in our case is approximately
## min(vw) ~0.15. There are only 219 obs with phat > 0.15.

while (T) {
  i <- sample(1:n, 1)
  qqnorm(rvw[i,], pch = 20)
  title(i)
  bla <- readline("\nHit return to see next curve or any character to exit")
  if (bla != "") break
}
## Most plots look very discrete. In low density regions, estimated density
## takes essentially only two values, depending on whether or not the obs is
## in the bootstrap sample.

## Let's look at the vertices with high vertex weights

ovw <- rev(order(mean.rvw))
ovw <- order(mean.rvw)
ylim = c(0, max(rvw))
for (i in 1:n) {
  qqnorm(rvw[ovw[i],], ylim = ylim, pch = 20, main = i)
  bla <- readline("\nHit return to see next curve or any character to exit")
  if (bla != "") break
}

##=================================================================
## Clustering with confidence
## First try the original idea:
## - Calculate confidence band (simultaneous or non-simultaneous)
## - Construct cluster tree from upper confidence band
## - Use bootstrap rise as pruning criterion with threshold 0

target.cov.prob <- 0.8
band.type <- "gaussian"
band.type <- "empirical"
band.type = "poisson"
simultaneous = FALSE

confidence.bounds <-
  cwc.confidence.bounds(resample.kernel.ve.weights, target.cov.prob,
                        band.type = band.type, simultaneous = simultaneous)
cov.parameter <- confidence.bounds$cov.parameter
cov.parameter
cwc.simul.cov.prob(resample.kernel.ve.weights, cov.parameter,
                   band.type = band.type)
names(confidence.bounds)

uc.mdsfun <- cwc.make.upper.confidence.mdsfun(resample.kernel.ve.weights, confidence.bounds)
obs.lower <- confidence.bounds$lower.vw
## obs.lower
pruning.criterion <- "bootstrap.rise"
##pruning.criterion <- "size"
gsl.cluster.out <- gsl.cluster(X, uc.mdsfun, pruning.criterion = pruning.criterion,
                               pruning.threshold = 0, gsl.cluster.out = NULL,
                               assign.fluff = T,
                               obs.density.lower.confidence.bounds = obs.lower)
bla <- gsl.runt.stats(gsl.cluster.out)
bla[1:min(20, nrow(bla)),]

table(group.id, gsl.cluster.out$leaf.code)
rand.index(group.id, gsl.cluster.out$leaf.code)
names(gsl.cluster.out)


##-----------------------------------------------------------------


## For poisson bounds we have two options:
##
## 1. Pick a simultaneous coverage prob. Compute upper confidence bounds
## for edge weights and lower confidence bounds for vertex weights.
## Compute mast for upper.confidence.mdsfun. Use "bootstrap.rise" as the
## pruning criterion.
## 
## 2. Compute mast for resample.mean.mdsfun. Use runt.poisson.gamma
## as the pruning.criterion.

## Start with version 2 and compute values of runt.poisson.gamma and the
## corresponding coverage probs

resample.mean.mdsfun <- cwc.make.resample.mean.mdsfun(resample.kernel.ve.weights)
obs.density <- rep(0, n)
for (i in 1:n) obs.density[i] <- resample.mean.mdsfun(i, i)
mast <- gsl.mast(resample.mean.mdsfun)
mast.dendogram <- gsl.mast.dendogram(mast)
pc.out <- gsl.pruning.criteria(mast.dendogram, obs.density)
cct.out <- gsl.compute.cluster.tree(mast.dendogram, pc.out,
                                    pruning.criterion = "poisson.gamma",
                                    pruning.threshold = 0)
runt.poisson.gamma <- cct.out[,"runt.poisson.gamma"]
bla <- rev(sort(runt.poisson.gamma))[1:20]
bla

cct.out <- gsl.compute.cluster.tree(mast.dendogram, pc.out,
                                    pruning.criterion = "size",
                                    pruning.threshold = 0)
runt.size <- cct.out[,"runt.size"]
bla <- rev(sort(runt.size))[1:20]
bla

## For bandwidth = cv.mean, the two largest values of runt.poisson.gamma are
## [1] 0.375976562 0.151367188. Let's compute the corresponding coverage probs:

cwc.cov.prob(resample.kernel.ve.weights, bla[1], band.type = "poisson")
## 1

cwc.cov.prob(resample.kernel.ve.weights, bla[2], band.type = "poisson")
## 0


## Now version (2)
## For any target.cov.prob <= 1 we should get a tree with two leaves:

target.cov.prob <- 0.8
band.type <- "poisson"
poisson.confidence.bounds <-
  cwc.confidence.bounds(resample.kernel.ve.weights, target.cov.prob,
                        band.type = band.type)
cov.parameter <- poisson.confidence.bounds$cov.parameter
cov.parameter
cwc.cov.prob(resample.kernel.ve.weights, cov.parameter,
                   band.type = band.type)
upper.confidence.mdsfun <-
  cwc.make.upper.confidence.mdsfun(resample.kernel.ve.weights,
                                   poisson.confidence.bounds) 
obs.density.lower.confidence.bounds <- poisson.confidence.bounds$lower.vw
obs.density <- poisson.confidence.bounds$upper.vw

mast <- gsl.mast(upper.confidence.mdsfun)
mast.dendogram <- gsl.mast.dendogram(mast)
pc.out <- gsl.pruning.criteria(mast.dendogram, obs.density,
                               obs.density.lower.confidence.bounds)
cct.out <- gsl.compute.cluster.tree(mast.dendogram, pc.out,
                                    pruning.criterion = "bootstrap.rise",
                                    pruning.threshold = 0)
cct.out
## Indeed



##-----------------------------------------------------------------
## Compare pruning criteria

runt.size <- pmin(pc.out[,"left.size"], pc.out[,"right.size"])
runt.excess.mass <- pmin(pc.out[,"left.excess.mass"], pc.out[,"right.excess.mass"])
runt.rise <- pmin(pc.out[,"left.rise"], pc.out[,"right.rise"])
runt.bootstrap.rise <- pmin(pc.out[,"left.bootstrap.rise"],
                           pc.out[,"right.bootstrap.rise"])
runt.poisson.gamma <- pmin(pc.out[,"left.poisson.gamma"],
                           pc.out[,"right.poisson.gamma"])

plot(runt.size, runt.excess.mass, pch = 20)
plot(runt.size, runt.rise, pch = 20)
plot(runt.size, runt.bootstrap.rise, pch = 20)
plot(runt.size, runt.poisson.gamma, pch = 20)
plot(runt.rise, runt.poisson.gamma, pch = 20)
plot(runt.bootstrap.rise, runt.poisson.gamma, pch = 20)

hist(runt.poisson.gamma, col = "green")
runt.poisson.gamma
confidence.levels <- rep(0, n-1)
for (i in 1:(n-1))
  confidence.levels[i] <- cwc.cov.prob(resample.kernel.ve.weights,
                                             runt.poisson.gamma[i],
                                             band.type = "poisson")
confidence.levels
## All "0" except one "1"

plot(rev(sort(runt.rise))[1:20], pch = 20)
plot(rev(sort(runt.poisson.gamma))[1:20], pch = 20)
plot(rev(sort(runt.size))[1:20], pch = 20)


##-----------------------------------------------------------------
## For 1d data
## Plot edge and vertex weights from upper bootstrap confidence
## bounds

ylim <- c(0, max(diag(uc.mdsmat)))
plot(X[,1], diag(uc.mdsmat), ylim = ylim, pch = 20, type = "n")
for (i in 1:n.edges) {
  v1 <- tg.edges[i, 1]
  v2 <- tg.edges[i, 2]
  lines(c(X[v1, 1], X[v2, 1]), rep(uc.mdsmat[v1, v2], 2))
}
points(X[,1], obs.density.lower.confidence.bounds, pch = 20, col = "red")


##-----------------------------------------------------------------
## For 2d data
## For given level, plot all edges with (upper confidence limit) of
## edge weight >= level and all vertices with lower confidence bound
## of vertex weight >= level

level <- 0.0144   ## banana-clump
level <- 0.053    ## bullseye

plot(X, pch = 20)
points(X[(obs.density.lower.confidence.bounds >= level),], pch = 20,
       col = "red")
for (i in 1:(n-1)) {
  v1 <- mast$edges[i,1]
  v2 <- mast$edges[i,2]
  edge.weight <- upper.confidence.mdsfun(v1, v2)
  if (edge.weight >= level) lines(rbind(X[v1,], X[v2,]))
}

##-----------------------------------------------------------------
## Look at bootstrap edge and vertex weights
##

par(mfrow = c(2,1))
i <- sample(1:nrow(bew), 1)
hist(bew[i,], col = "green", nclass = 15)
hist(log(bew[i,]), col = "green", nclass = 15)

par(mfrow = c(2,1))
i <- sample(1:nrow(bvw), 1)
hist(bvw[i,], col = "green", nclass = 15)
hist(log(bvw[i,]), col = "green", nclass = 15)


##=================================================================
##=================================================================
##=================================================================
##=================================================================
## Chapter 4: Simulations to assess level and power of CWC
##
## For uniform and Gaussian, plot true vs nominal significance level:
##
## - For each replication i, find smallest alpha^n_i for which CWC with
##   non-simultaneous confidence level 1-alpha_i produces non-trivial tree.
##   Determine corresponding simultaneous alpha^s_i.
## - If calibration was perfect, the fraction of replications with
##   alpha_i < alpha should be alpha


## Two approaches for pruning:
## ------------------------------
## 1. Hypotheses testing:
## - Pick a criterion (size, excess mass, or persistence)
## - Simulate distribution of criterion for uniform distribution
##   (or some other unimodal distribution that is closest to the
##   data distribution)
## - For given significance level alpha, find (1-alpha) quantile
##   gamma of criterion. Use gamma as the pruning threshold.
##
## 2. Clustering with confidence

## What takes the time is computing the bootstrap vertex and edge
## weights. We want to save them so that we can, for example, compute
## trees for different confidence levels. We could easily end up with
## thousands of files

## (5-26-2010): Use AWS to compute bootstrap vertex and edge
## weights. Try to do remaining analysis locally. (This can always be
## changed later).

## The file name of each data file and each weight file should
## uniquely identify all the experimental variables. Each file should
## be a dump file that is self descriptive

## For each data / bootstrap samples file we need to know
## - distribution: sampling distribution from which data set was generated
## - n: sample size
## - m: dimension
## - resampling.type (half or boot)
## - sample.id

## For each bootstrap weight file we also need to know
## - kmax
## - include.mst
## - ngrid
## - bandwidth

## Each simulation experiment consists of several steps:
## 1. For each replication, generate an R dump file with defining data
##    matrix X and bootstrap samples matrix bootstrap.samples;
## 2. compute bootstrap vertex and edge weights;
## 3. Find




## What takes the time is computing the bootstrap vertex and edge
## weights. We want to save them so that we can, for example, compute
## trees for different confidence levels.  We could easily end up with
## thousands of files

## cwc.boot.weights <- function(X, kmax, bootstrap.samples = NULL, nboot =
##                              NULL, trial.h = NULL, ngrid = 10,
##                              bandwidth = "cv.mean", resampling.type =
##                              "half.sampling", debug = F) {
##   if (is.null(trial.h))
##     trial.h <- c(0.01, 0.05, 0.10, 0.15, 0.20, 0.30, 0.40, 0.60, 0.80,
##                  1.00, 1.50, 2.00)
##     tg.edges <- cwc.determine.edges(X, kmax = kmax, include.mst = T)

##   n.edges <- nrow(tg.edges)

## if (file.exists(bootstrap.weights.pathname))
##   source(bootstrap.weights.pathname, echo = T)

## if (!file.exists(bootstrap.weights.pathname)) {
##   bootstrap.weights <- cwc.boot.kernel.ve.weights(X, tg.edges, X, nboot, 
##                                          trial.h = trial.h, interactive = T,
##                                          bandwidth = "cv.mean")
##   dump("bootstrap.weights", bootstrap.weights.pathname)
## }


## ## Determining level: For each sample, find largest non- simutaneous
## ## confidence level that results in a split, and the corresponding
## ## simultaneous confidence level. Do that for empirical and Gaussian
## ## confidence intervals.

## ## Parameters:
## ## - k (for nearest neighbor graph)
## ## - bandwidth (for bandwidth selection)
## ## - bootstrap.samples (a matrix with nboot columns and either
## ##   n or n/2 rows, depending on whether we sample with or without
## ##   replacement

## ##-----------------------------------------------------------------
## ## Here is the rscript:
## ## We assume that the input file (copied from S3) is an R dump file
## ## defining
## ## X          data matrix

## args <- commandArgs(TRUE)
## input.filename <- args[1]



## nboot <- args[2]
## data <- readLines(con = data.filename)
## input.filename <- 



###################################################################  
## Chapter 5: Try idea for assessing significance of splits

trial.h <- c(0.01, 0.05, 0.10, 0.15, 0.20, 0.30, 0.40, 0.60, 0.80, 1.00,
             1.50, 2.00)
kmax <- 20
n.noot <- 100
boot.bandwidth = "cv.mean"

data.name <- "nice-univariate-trimodal"
kmax <- 1

data.name <- "banana-clump-cv"
kmax <- 20

data.name <- "uniform-n50-m1"
kmax <- 1

data.name <- "uniform-n100-m2"
kmax <- 20

data.name <- "uniform-n250-m5"
kmax <- 20

data.name <- "uniform-n500-m10"
kmax <- 20
## The bootstrap distribution of split edge weights is clearly
## bimodal for most of the nodes. Conjecture: Split edge weight
## is high if one or both end points are in bootstrap sample.
##

data.name <- "gaussian-n50-m1"
kmax <- 1

data.name <- "gaussian-n100-m2"
kmax <- 20

data.name <- "gaussian-n250-m5"
kmax <- 20

data.name <- "gaussian-n500-m10"
kmax <- 20

## Seems like there are problems in dimension 10 even with
## gaussian data. For n = 500, m = 10, runt size threshold = 3
## we get
## sig.ass = [1] 1.00 0.63 0.77 0.84 0.96

data.name <- "bullseye"
kmax <- 20
## Three large runt sizes -> three highly significant splits

data.name <- "olive"
kmax <- 20
## Doesn't work. Results in many highly significant splits. Split
## edge weights are discrete for some splits
## Maybe partly a multiplicity problem?

data.filename <- paste(data.name, ".R", sep = "")
data.pathname <- paste(data.dir, data.filename, sep = dir.sep)
source(data.pathname,echo = T)
X <- sphere(X)
diag(var(X))

mdsmat.filename <- paste(data.name, "-mdsmat.R", sep = "")
mdsmat.pathname <- paste(data.dir, mdsmat.filename, sep = dir.sep)
if (file.exists(mdsmat.pathname))
  source(mdsmat.pathname, echo = T)

if (!file.exists(mdsmat.pathname)) {
  ## trial.h <- c(0.01, 0.05, 0.10, 0.15, 0.20, 0.30, 0.40, 0.60, 0.80, 1.00, 1.50, 2.00)
  cvs.out <- cv.search(X, trial.par = trial.h)
  cvs.out
  h <- cvs.out$opt.smopar
  density <- make.gaussian.kernel.density.estimate(X, h)
  mdsmat <- gsl.mdsmat(X, density, kmax = kmax, interactive = T)
  dump("mdsmat", mdsmat.pathname)
}

mdsfun <- gsl.make.mdsfun.from.mdsmat(mdsmat)
gsl.cluster.out <- gsl.cluster(X, mdsfun, pruning.criterion = "size",
                        pruning.threshold = 0, assign.fluff = F)
gsl.runt.size(gsl.cluster.out)[1:20]

bootstrap.weights.filename <- paste(data.name, "-bootstrap-weights.R", sep = "")
bootstrap.weights.pathname <- paste(data.dir, bootstrap.weights.filename,
                                    sep = dir.sep)
if (file.exists(bootstrap.weights.pathname))
  source(bootstrap.weights.pathname, echo = T)

if (!file.exists(bootstrap.weights.pathname)) {
  mast.dendogram <- gsl.cluster.out$mast.dendogram
  mast.edges <- mast.dendogram$edges
  bootstrap.weights <-
    cwc.boot.kernel.ve.weights(X, mast.edges, X, nboot, 
                               trial.h = trial.h, bandwidth = boot.bandwidth)
  dump("bootstrap.weights", bootstrap.weights.pathname)
}

bootstrap.vertex.weights <- bootstrap.weights$bootstrap.vertex.weights
bootstrap.edge.weights <- bootstrap.weights$bootstrap.edge.weights

bootstrap.kernel.mdsfun <-
  cwc.make.bootstrap.kernel.mdsfun(mast.edges, bootstrap.vertex.weights,
                                   bootstrap.edge.weights)

gsl.runt.size(gsl.cluster.out)[1:20]

pruning.threshold <- 20
gsl.cluster.out <- gsl.cluster(X, mdsfun, pruning.criterion = "size",
                               pruning.threshold = pruning.threshold,
                               assign.fluff = T)

##-----------------------------------------------------------------
## Visually select cluster tree nodes for significance assessment
inode <- gsl.draw.tree.pick.node(gsl.cluster.out, tree.type = "cluster.tree") 

uniform.level <- T
ass.out <- cwc.assess.split.significance(gsl.cluster.out, inode,
                                         bootstrap.kernel.mdsfun,
                                         uniform.level = uniform.level)
split.edge.weights <- ass.out$split.edge.weights
left.maxes <- ass.out$left.maxes
right.maxes <- ass.out$right.maxes
min.core.maxes <- apply(rbind(left.maxes, right.maxes), 2, min)
plot(split.edge.weights, min.core.maxes, pch = 20)
opt.level <- ass.out$opt.level
abline(h = opt.level)
abline(v = opt.level)
opt.level
supporting.bs <- ass.out$supporting.bs
supporting.bs
length(supporting.bs) / nboot
inode

##-----------------------------------------------------------------
## For 10d uniform, where split edge weights are bimodal, calculate
## how many edges are in each mode.
## inode = 1
sum(split.edge.weights < 2e-6)  ## 74
## inode = 2
sum(split.edge.weights < 4e-6)  ## 81
## inode = 3
sum(split.edge.weights < 4e-6)  ## 77
## inode = 4
sum(split.edge.weights < 4e-6)  ## 73

## What is the probability that at least one of two observations
## is in the bootstrap sample?
## P = 1- (1 - 2/nboot)^nboot ~ 1 - e^(-2)
1 - exp(-2)  ## 0.86
## Plausible but not totally convincing


##-----------------------------------------------------------------
## Assess significance for all cluster tree nodes with runt size
## >= rs.threshold and runt.excess.mass >= rem.threshold

rs.threshold <- 3
gsl.cluster.out <- gsl.cluster(X, mdsfun, pruning.criterion = "size",
                               pruning.threshold = rs.threshold,
                               assign.fluff = T)
rem.threshold <- 0
cluster.tree <- gsl.cluster.out$cluster.tree
n.node <- nrow(cluster.tree)
runt.size <- cluster.tree[, "runt.size"]
runt.excess.mass <- cluster.tree[, "runt.excess.mass"]
ors <- rev(order(runt.size))
sorted.inode <- (1:n.node)[ors]
sorted.rs <- runt.size[ors]
sorted.rem <- runt.excess.mass[ors]
ass <- (sorted.rs >= rs.threshold) & (sorted.rem >= rem.threshold)
n.ass <- sum(ass)
sorted.inode.ass <- sorted.inode[ass]
sorted.rem.ass <- sorted.rem[ass]
sorted.rs.ass <- sorted.rs[ass]
sig.ass <- rep(0, n.ass)

for (i in 1:n.ass) {
  ass.out <- cwc.assess.split.significance(gsl.cluster.out, sorted.inode[i],
                                           bootstrap.kernel.mdsfun,
                                           uniform.level = uniform.level)
  sig.ass[i] <- length(ass.out$supporting.bs) / nboot
}
sig.ass
plot(sorted.rem.ass, sig.ass, pch = 20)
plot(sorted.rs.ass, sig.ass, pch = 20)
plot(sorted.rem.ass, sorted.rs.ass, pch = 20)

## This idea for testing significance of splits may have some potential
## for low-dimensional data, but does not seem to work for high dimensions.


























##=================================================================
##=================================================================
##=================================================================
##=================================================================
## Chapter 6: Now the same thing for 1-nn estimate

density <- cwc.make.nn.density.estimate(X)

## Compute vertex and edge weights for test graph. Vertices equally spaced
## in range of data

n.vertices <- 100
tg.vertices <- as.matrix(seq(min(X), max(X), length = n.vertices), ncol = 1)
tg.edges <- cbind((1:(n.vertices - 1)), (2:n.vertices))

## For 1-nn estimate, it's probably better to choose the observations as
## vertices

n.vertices <- nrow(X)
n.edges <- n.vertices - 1
tg.vertices <- as.matrix(sort(X), ncol = 1)
tg.edges <- cbind((1:(n.vertices - 1)), (2:n.vertices))

ngrid <- 10
cwcvw.out <- cwc.ve.weights(tg.vertices, tg.edges, density, ngrid)
vw <- cwcvw.out$vertex.weights
ew <- cwcvw.out$edge.weights

## Plot edge weights
ylim <- c(0, max(sqrt(ew)))
plot(0, 0, xlim = range(tg.vertices), ylim = ylim, type = "n")
for (i in 1:nrow(tg.edges)) {
  e <- tg.edges[i,]
  v1 <- tg.vertices[e[1]]
  v2 <- tg.vertices[e[2]]
  lines(c(v1, v2), sqrt(rep(ew[i], 2)), lwd = 3)
}

plot(tg.vertices[1:(n.vertices - 1)], log(ew), type = "l")
plot(tg.vertices[1:(n.vertices - 1)], ew, type = "l")

## The edge weights look prett wild. Note, though, how the smallest edge
## weights are in low density regions.

##-----------------------------------------------------------------
## Compute bootstrap edge weights

nboot <- 50
cwcbvw.out <- cwc.boot.nn.ve.weights(tg.vertices, tg.edges, X, nboot)
bvw <- cwcbvw.out$bootstrap.vertex.weights
bew <- cwcbvw.out$bootstrap.edge.weights
dim(bew)

i <- 30
qqnorm(bew[i,])

## Not suprisingly, bootstrap edge weights are highly discrete

bew.means <- apply(bew, 1, mean)

## Plot mean edge weights and bootstrap edge
ylim <- c(0, max(sqrt(bew.means)))
plot(0, 0, xlim = range(tg.vertices), ylim = ylim, type = "n")
for (i in 1:nrow(tg.edges)) {
  e <- tg.edges[i,]
  v1 <- tg.vertices[e[1]]
  v2 <- tg.vertices[e[2]]
  lines(c(v1, v2), sqrt(rep(bew.means[i], 2)), lwd = 3)
}
for (i in 1:nboot) points(tg.vertices[1:(n-1),], bew[,i], pch = ".")

## Plot upper and lower order stats of edge weights


midpoints <- matrix(0, nrow = n.edges, ncol = ncol(tg.vertices))
for (i in 1:ne) {
  midpoints[i,] <- 0.5 * (tg.vertices[tg.edges[i, 1],] +
                          tg.vertices[tg.edges[i, 2],])
}

trim <- 1
upper <- rep(0, n.edges)
lower <- rep(0, n.edges)
for (i in 1:n.edges) {
  sw <- sort(bew[i,])
  upper[i] <- sw[nboot  - trim]
  lower[i] <- sw[1 + trim]
}

tlower <- log(lower)
tupper <- log(upper)
ylim <- range(c(tlower, tupper))
plot(midpoints, tupper, ylim = ylim, type = "l", col = "red")
lines(midpoints, tlower, ylim = ylim, type = "l", col = "black")

##-----------------------------------------------------------------
## Try runt pruning
  
mdsfun <- gsl.make.nn.mdsfun(X)
gsl.cluster.out <- gsl.cluster(X, mdsfun, pruning.criterion = "size",
                        pruning.threshold = 0, assign.fluff = F)
gsl.runt.size(gsl.cluster.out)[1:10]
##  [1] 184 100  21  10   0   0   0   0   0  NA
gsl.cluster.out <- gsl.cluster(X, mdsfun, pruning.criterion = "size",
                        pruning.threshold = 19, assign.fluff = T)

cluster.id <- gsl.observation.labels(gsl.cluster.out)
unique(cluster.id)
in.core <- gsl.observation.in.core(gsl.cluster.out)

## -----------------------------------------------------------------
## Make plots just to be sure

SX <- as.matrix(X[order(X[,1]),])
Sid <- cluster.id[order(X[,1])]
Sic <- in.core[order(X[,1])]

plot (midpoints, sqrt(ew), ylim = c(0, max(sqrt(ew))), type = "l")
points(SX[(Sid == 3),],
       rep(0, sum(Sid == 3)), col = "black", pch = 20)
points(SX[((Sid == 3) & Sic),],
       rep(0, sum((Sid == 3) & Sic)), col = "black", pch = 19)
points(SX[(Sid == 4),],
       rep(0, sum(Sid == 4)), col = "red", pch = 20)
points(SX[((Sid == 4) & Sic),],
       rep(0, sum((Sid == 4) & Sic)), col = "red", pch = 19)
points(SX[(Sid == 5),],
       rep(0, sum(Sid == 5)), col = "green", pch = 20)
points(SX[((Sid == 5) & Sic),],
       rep(0, sum((Sid == 5) & Sic)), col = "green", pch = 19)

#################################################################
#################################################################
#################################################################

#################################################################

mdsmat.filename <- paste(data.name, "-mdsmat.R", sep = "")
mdsmat.pathname <- paste(data.dir, mdsmat.filename, sep = dir.sep)
if (file.exists(mdsmat.pathname))
  source(mdsmat.pathname, echo = T)

if (!file.exists(mdsmat.pathname)) {
  ## trial.h <- c(0.01, 0.05, 0.10, 0.15, 0.20, 0.30, 0.40, 0.60, 0.80, 1.00, 1.50,
  ## 2.00)
  cvs.out <- cv.search(X, trial.h)
  cvs.out
  h <- cvs.out$opt.smopar
  density <- make.gaussian.kernel.density.estimate(X, h)
  mdsmat <- gsl.mdsmat(X, density, kmax = kmax, interactive = T)
  dump("mdsmat", mdsmat.pathname)
}
mdsfun <- gsl.make.mdsfun.from.mdsmat(mdsmat)



  ##-----------------------------------------------------------------
## Maybe we should compute the "point estimate" tree from the mean of
## the resample weights, rather than from the weights for the original
## sample.
##
## Or maybe we should resample using bootstrap rather than half-sampling.
## 


## Some exploratory plots to better understand these pictures

nboot <- ncol(bew)
ne <- length(ew)
nv <- length(vw)
z.cov <- qnorm(1 - (1-confidence.level)/2)

nout <- ceiling(nboot * (1-confidence.level)/2)
lower <- nout + 1
upper <- nboot - nout
if (lower >= upper) {
  lower <- lower - 1
  upper <- upper + 1
}
upper.empirical.ew <- rep(0, ne)
lower.empirical.vw <- rep(0, nv)
for (i in 1:ne) {
  sew <- sort(bew[i,])
  upper.empirical.ew[i] <- sew[upper]
}
for (i in 1:nv) {
  svw <- sort(bvw[i,])
  lower.empirical.vw[i] <- svw[lower]
}
ew.mean <- apply(bew, 1, mean)
ew.sd <- apply(bew, 1, sd)
vw.mean <- apply(bvw, 1, mean)
vw.sd <- apply(bvw, 1, sd)
upper.gaussian.ew <- ew.mean + z.cov * ew.sd
lower.gaussian.vw <- vw.mean - z.cov * vw.sd

plot(upper.empirical.ew, upper.gaussian.ew, pch = ".")
cor(upper.empirical.ew, upper.gaussian.ew)

plot(ew.mean, upper.empirical.ew, pch = ".")


gamma.plus <- (sum(upper.empirical.ew * sqrt(ew.mean)) -
               sum(ew.mean^(3/2))) / sum(ew.mean)
gamma.minus <- (sum(lower.empirical.vw * sqrt(vw.mean)) -
                sum(vw.mean^(3/2))) / sum(vw.mean)
upper.poisson.ew <- ew.mean + gamma.plus * sqrt(ew.mean)
lower.poisson.vw <- vw.mean + gamma.minus * sqrt(vw.mean)

pch = "."
plot(ew.mean, upper.empirical.ew - ew.mean, pch = pch)
points(ew.mean, upper.gaussian.ew - ew.mean, pch = pch, col = "orange")
oew <- order(ew.mean)
lines(ew.mean[oew], upper.poisson.ew[oew] - ew.mean[oew], col = "red", lwd = 3)
title("Upper confidence bounds for edge weights")
plot(vw.mean, lower.empirical.vw - vw.mean, pch = 20)
points(vw.mean, lower.gaussian.vw - vw.mean, pch = 20, col = "orange")
ovw <- order(vw.mean)
lines(vw.mean[ovw], lower.poisson.vw[ovw] - vw.mean[ovw], col = "red", lwd = 3)
title("Lower confidence bounds for vertex weights")
par(mfrow = c(1,1))

plot(ew, ew.mean, pch = ".")
cor(ew, ew.mean)
plot(vw, vw.mean, pch = ".")
cor(vw, vw.mean)

nin.empirical <- nboot
nin.gaussian <- nboot
nin.poisson <- nboot
for (i in 1:nboot) {
  out.empirical <- sum(bew[,i] > upper.empirical.ew) + sum(bvw[,i] < lower.empiricalvw)
  out.gaussian <- sum(bew[,i] > upper.gaussian.ew) + sum(bvw[,i] < lower.gaussian.vw)
  out.poisson <- sum(bew[,i] > upper.poisson.ew) + sum(bvw[,i] < lower.poisson.vw)
  if (sum(out.empirical) > 0) nin.empirical <- nin.empirical - 1
  if (sum(out.gaussian) > 0) nin.gaussian <- nin.gaussian - 1
  if (sum(out.poisson) > 0) nin.poisson <- nin.poisson - 1
}
return(list(empirical.cov.prob = nin.empirical / nboot,
            simul.gaussian.cov.prob = nin.gaussian / nboot,
            simul.poisson.cov.prob = nin.poisson / nboot))



