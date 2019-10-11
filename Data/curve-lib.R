# The code in wxs\Pricur\Bellcore\curve-lib.S, modified for R
# (12-6-02)
#=================================================================

library("splines")

define.curve <- function(nseg = 100, loop = FALSE) {
  cat("\nPick locations by pointing and clicking left\n")
  cat("Finish by clicking right and choosing `stop'")
  plot(c(0,1), c(0,1), type = "n")
  points <- locator(type = "p", pch = 1)
  x <- cbind(points$x, points$y)
  if (loop) x <- rbind(x, x[1,])
  n <- nrow(x)
  diff <- x[2:n,] - x[1:(n-1),] 
  distances <- sqrt(apply(diff**2, 1, sum))
  arc.length <- cumsum(distances)
  arc.length <- c(0, arc.length)
  arc.length <- arc.length / max(arc.length)
  if (loop) {
    x1.spline <- interpSpline(arc.length, x[,1], period = 1)
    x2.spline <- interpSpline(arc.length, x[,2], period = 1)
  }
  else {
    x1.spline <- interpSpline(arc.length, x[,1])
    x2.spline <- interpSpline(arc.length, x[,2])
  }
  lambda <- seq(0, 1, length = nseg)
  lines(predict(x1.spline, lambda)$y, predict(x2.spline, lambda)$y)
  return(list(x1.spline = x1.spline, x2.spline = x2.spline))
}

#-----------------------------------------------------------------

create.points.around.curve <- function (curve, nobs, sd) {
  x1.spline <- curve$x1.spline
  x2.spline <- curve$x2.spline
  lambda <- runif(nobs)
  x1.c <- predict(x1.spline, lambda)$y
  x2.c <- predict(x2.spline, lambda)$y
  x1.obs <- x1.c + rnorm(nobs, 0, sd)
  x2.obs <- x2.c + rnorm(nobs, 0, sd)
  x <- cbind(x1.obs, x2.obs)
  return(x)
}

#-----------------------------------------------------------------

plot.curve <- function(curve, nseg = 100) {
  lambda <- seq(0, 1, length = nseg)
  lines(predict(curve$x1.spline, lambda)$y,
        predict(curve$x2.spline, lambda)$y)
}

