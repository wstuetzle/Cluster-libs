## Figures for CWC paper (7-17-2010) (10-20-2017)
## ==============================================

## OS <- "Windows"
## OS <- "Linux"
OS <- "Mac"

if (OS == "Mac") {
  cwc.dir <- "/Users/wxs/Dropbox/Clustering/CWC-before-Germany-2011"
  dir.sep <- "/"
  quartz()
}


## if (OS == "Windows") {
##   cwc.dir <- "Z:\\Werner\\Clustering\\CWC-before-Germany-2011"
##   dir.sep <- "\\"
## }

## if (OS == "Linux") {
##   cwc.dir <- "/homes/wxs/Werner/Clustering/CWC-before-Germany-2011"
##   dir.sep <- "/"
## }

code.dir <- paste(cwc.dir, "Code", sep = dir.sep)
data.dir <- paste(cwc.dir, "Data", sep = dir.sep)

echo <- T

source(paste(code.dir, "gsl-functions-5-28-2010.R", sep = dir.sep), echo = echo)

source(paste(code.dir, "cwc-functions-8-28-2010-v2017.R", sep = dir.sep), echo = echo)
source(paste(code.dir, "cwc-functions-5-28-2010.R", sep = dir.sep), echo = echo)
source(paste(code.dir, "cwc-functions-7-13-2010.R", sep = dir.sep), echo = echo)

##-----------------------------------------------------------------
## ## figures.dir <- paste(cwc.dir, "Paper\\Figures", sep = "\\")

## cwc.dir <- "Z:\\Werner\\Clustering\\CWC"
## ## cwc.dir <- "C:\\Werner\\Clustering\\CWC"
## code.dir <- paste(cwc.dir, "Code", sep = "\\")
## data.dir <- paste(cwc.dir, "Data", sep = "\\")
## figures.dir <- paste(cwc.dir, "Paper\\Figures", sep = "\\")

## cwc.functions.filename <- "cwc-functions-7-13-2010.R"
## cwc.functions.pathname <- paste(code.dir, cwc.functions.filename, sep = "\\")

## echo = TRUE

## source(cwc.functions.pathname, echo = F)
##-----------------------------------------------------------------

##=================================================================
## Start with illustrating idea for univariate data using nice
## univariate trimodal.

##-----------------------------------------------------------------
## Read or generate data

data.name <- "nice-univariate-trimodal"
data.filename <- paste(data.name, ".R", sep = "")
data.pathname <- paste(data.dir, data.filename, sep = dir.sep)
data.pathname

if(file.exists(data.pathname)) source(data.pathname, echo = echo)
n <- nrow(X)

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
## Read or generate original and resample vertex and edge
## weights. Pictures look better if we choose equi-spaced vertices for
## the test graph

resample.type <- "half.sampling"
n.resamples <- 100 
kmax <- 1
ngrid <- 10
bandwidth <- "cv.mean"

rw.filename <- paste("equi-",
                     resample.weights.filename(data.name, resample.type,
                                               n.resamples, bandwidth,
                                               kmax, ngrid), sep = "")

rw.pathname <- paste(paste(data.dir, rw.filename, sep = dir.sep), ".R", sep = "")
rw.pathname

if (file.exists(rw.pathname)) {
  source(rw.pathname, echo = echo)
  tg.vertices <- resample.kernel.ve.weights$tg.vertices
}
if (!file.exists(rw.pathname)) {
  n.vertices <- 100
  nv <- n.vertices
  ne <- nv - 1
  tg.vertices <- as.matrix(seq(min(X), max(X), length = n.vertices), ncol = 1)
  tg.edges <- cbind((1:(n.vertices - 1)), (2:n.vertices))
  resample.kernel.ve.weights <-
    cwc.resample.kernel.ve.weights(X, tg.vertices = tg.vertices,
                                   tg.edges = tg.edges, kmax = kmax, 
                                   n.resamples = n.resamples, ngrid = ngrid,
                                   bandwidth = bandwidth,
                                   resample.type = resample.type,
                                   debug.filename = stdout())
  dump("resample.kernel.ve.weights", file = rw.pathname)
}
  

##-----------------------------------------------------------------
## Figure 1
## Show points, density estimate, halfsample estimates (with same
## bandwidth?).

original.vertex.weights <- resample.kernel.ve.weights$vertex.weights
resample.vertex.weights <- resample.kernel.ve.weights$resample.vertex.weights

## ylim <- range(cbind(original.vertex.weights, resample.vertex.weights))
ylim <- c(0, max(cbind(original.vertex.weights, resample.vertex.weights)))
plot(tg.vertices, original.vertex.weights, ylim = ylim, type = "l",
     lwd = 3, xaxt = "n", yaxt = "n", xlab = "", ylab = "")
for (i in 1:n.resamples) {
  points(tg.vertices, resample.vertex.weights[,i], pch = ".")
}
points(X, rep(0, n), pch = 20)

##-----------------------------------------------------------------
## Figure 2
## Plot sd of resample vertex weights vs mean

rvw.means <- apply(resample.vertex.weights, 1, mean)
rvw.sds <- apply(resample.vertex.weights, 1, sd)
coef <- lsfit(sqrt(rvw.means), rvw.sds, intercept = F)$coef
plot(rvw.means, rvw.sds, pch = 20, xlab = "mean of resample estimates",
     ylab = "standard deviation of resample estimates")
x <- sort(rvw.means)
lines(x, coef * sqrt(x), lwd = 3)

##-----------------------------------------------------------------
## Figure 3
## Plot poisson confidence bands

target.simul.cov.prob = 0.8
band.type = "poisson"

confidence.bounds <-
  cwc.confidence.bounds(resample.kernel.ve.weights, target.simul.cov.prob,
                        band.type = band.type)
confidence.bounds
cov.parameter <- confidence.bounds$cov.parameter
scp.out <- cwc.simul.cov.prob(resample.kernel.ve.weights, cov.parameter,
                              band.type = band.type)
scp.out
scp.out$simul.cov.prob
out.resamples <- scp.out$out.resamples
out.resamples
ylim <- c(0, max(resample.vertex.weights))
plot(tg.vertices, original.vertex.weights, ylim = ylim, 
     lwd = 3, xaxt = "n", yaxt = "n", xlab = "", ylab = "", type = "n")
for (i in 1:n.resamples) {
  points(tg.vertices, resample.vertex.weights[,i], pch = ".")
}
lines(tg.vertices, confidence.bounds$upper.vw, lwd= 3)
lines(tg.vertices, confidence.bounds$lower.vw, lwd= 3)
## for (i in out.resamples)
##   points(tg.vertices, resample.vertex.weights[,i], pch = ".", col = "red")

## We could add empirical bands but there is almost no difference.
band.type = "empirical"
confidence.bounds <-
  cwc.confidence.bounds(resample.kernel.ve.weights, target.simul.cov.prob,
                                  band.type = band.type)
lines(tg.vertices, confidence.bounds$upper.vw, lwd= 3, col = "red")
lines(tg.vertices, confidence.bounds$lower.vw, lwd= 3, col = "red")

## loc <- locator(1, type = "p")
dip1.x <- 1.956717
dip1.y <- 0.08683574
abline(v = dip1.x, lty = 2, lwd = 2)
abline(h = dip1.y, lty = 2, lwd = 2)

## loc <- locator(1, type = "p")
## dip2.x <- 0.5247248
## dip2.y <- 0.3961393
## abline(v = dip2.x, lty = 3, lwd = 2)
## abline(h = dip2.y, lty = 3, lwd = 2)


## Point to make: bands are very wide. High only becomes significant for
## fairly small simultanous coverage prob. Explain why level of valley
## is important
 

