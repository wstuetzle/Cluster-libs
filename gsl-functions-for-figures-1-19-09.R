## R code for auxiliary function used in producing the figures for "A
## generalized single linkage method for estimating the cluster tree
## of a density" by Werner Stuetzle and Rebecca Nugent (1-19-09)
## =================================================================

gsl.define.curve <- function(nseg = 100, loop = FALSE) {
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

gsl.create.points.around.curve <- function (curve, nobs, sd) {
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

gsl.plot.curve <- function(curve, nseg = 100) {
  lambda <- seq(0, 1, length = nseg)
  lines(predict(curve$x1.spline, lambda)$y,
        predict(curve$x2.spline, lambda)$y)
}

##-----------------------------------------------------------------

gsl.plot.2d.density <- function(density, X, plot.title.prefix) {
  cex.main <- 0.8
  edog.out <- gsl.evaluate.density.on.grid(density, X)
  x1.grid <- edog.out$x1.grid
  x2.grid <- edog.out$x2.grid
  z <- edog.out$z
  persp(x1.grid, x2.grid, sqrt(z), box = F, theta = 45, phi = 30,
        expand = 0.5, shade = NA)
  tit <- paste(plot.title.prefix, ", sqrt(density)", sep = "")
  title(tit, cex.main = cex.main)
  contour(x = x1.grid, y = x2.grid, sqrt(z))
  points(X[,1:2], pch = 20)
  tit <- paste(plot.title.prefix, ", sqrt(density)", sep = "")
  title(tit, cex.main = cex.main)
  image(x = x1.grid, y = x2.grid, sqrt(z))
  points(X[,1:2], pch = 20)
  tit <- paste(plot.title.prefix, ", sqrt(density)", sep = "")
  title(tit, cex.main = cex.main)
}

gsl.evaluate.density.on.grid <- function(density, X, ngrid = 50) {
  x1.grid <- seq(from = min(X[,1]), to = max(X[,1]), length = ngrid)
  x2.grid <- seq(from = min(X[,2]), to = max(X[,2]), length = ngrid)
  Grid <- matrix(0, nrow = ngrid^2, ncol = 2)
  Grid[,1] <- rep(x1.grid, rep(ngrid, ngrid))
  Grid[,2] <- rep(x2.grid, ngrid)
  if (ncol(X) > 2) {
    Grid <- cbind(Grid, matrix(0, nrow = ngrid^2, ncol = ncol(X) - 2))
  }
  z <- matrix(density(Grid), nrow = ngrid, ncol = ngrid, byrow = T)
  return(list(x1.grid = x1.grid, x2.grid = x2.grid, z = z))
}

##=================================================================
## The mast of a minimum density similarity function will in general
## not be unique. This makes no difference for the clustering; any
## mast will do. To obtain a mast that does not look like spaghetti, we
## break ties in edge weights using Euclidean distance. Here is a
## (very slow) version of the function cwcint.mst that does the
## tiebreak. 

gsl.mst.with.tiebreak <- function(Dfun, X) {
  n <- attr(Dfun, "n")
  D <- as.matrix(dist(X))
  edges <- matrix(0, nrow = n-1, ncol = 2)
  edge.weights <- rep(0, n-1)
  out.points <- 1:(n-1)
  closest.in.points <- rep(n, (n-1))
  dist.to.closest.in.point <- Dfun(n, 1:(n-1))
  euclidean.dist.to.closest.in.point <- D[n, 1:(n-1)]
  for (i in (n-1):2) {
    jmin <- 1
    new.in.point <- out.points[1]
    closest.old.in.point <- closest.in.points[1]
    smallest.dist <- dist.to.closest.in.point[1]
    euclidean.dist <- euclidean.dist.to.closest.in.point[1]
    for (j in 2:i) {
      if ((dist.to.closest.in.point[j] < smallest.dist) ||
          ((dist.to.closest.in.point[j] == smallest.dist) &
           (euclidean.dist.to.closest.in.point[j] < euclidean.dist))) {
        jmin <- j
        new.in.point <- out.points[jmin]
        closest.old.in.point <- closest.in.points[jmin]
        smallest.dist <- dist.to.closest.in.point[jmin]
        euclidean.dist <- euclidean.dist.to.closest.in.point[jmin]
      }
    }
    edges[i, 1] <- new.in.point
    edges[i, 2] <- closest.old.in.point
    edge.weights[i] <- smallest.dist
    out.points[jmin] <- out.points[i]
    closest.in.points[jmin] <- closest.in.points[i]
    dist.to.closest.in.point[jmin] <- dist.to.closest.in.point[i]
    euclidean.dist.to.closest.in.point[jmin] <- euclidean.dist.to.closest.in.point[i]
    out.points <- out.points[1:(i-1)]
    closest.in.points <- closest.in.points[1:(i-1)]
    dist.to.closest.in.point <- dist.to.closest.in.point[1:(i-1)]
    euclidean.dist.to.closest.in.point <- euclidean.dist.to.closest.in.point[1:(i-1)]
    dist.to.new.in.point <- Dfun(new.in.point, out.points)
    euclidean.dist.to.new.in.point <- D[new.in.point, out.points]
    for (j in 1:(i-1)) {
      if((dist.to.new.in.point[j] < dist.to.closest.in.point[j]) ||
         ((dist.to.new.in.point[j] == dist.to.closest.in.point[j]) &
          (euclidean.dist.to.new.in.point[j] < euclidean.dist.to.closest.in.point[j]))) {
        closest.in.points[j] <- new.in.point
        dist.to.closest.in.point[j] <- dist.to.new.in.point[j]
        euclidean.dist.to.closest.in.point[j] <- euclidean.dist.to.new.in.point[j]
      }
    }
  }
  edges[1, 1] <- out.points[1]
  edges[1, 2] <- closest.in.points[1]
  edge.weights[1] <- dist.to.closest.in.point[1]
  return(list(edges = edges, edge.weights = edge.weights))
}

##-----------------------------------------------------------------

gsl.mast.with.tiebreak <- function(Sfun, X) {
  Dfun <- function(...) {
    return(- Sfun(...))
  }
  attr(Dfun, "n") <- attr(Sfun, "n")
  mst <- gsl.mst.with.tiebreak(Dfun, X)
  mast <- mst
  mast$edge.weights <- -mst$edge.weights
  return(mast)
}

##-----------------------------------------------------------------

gsl.cluster.with.tiebreak <- function(X, mdsfun, pruning.criterion = "size",
                                      pruning.threshold = 0) {
  n <- attr(mdsfun, "n")
  mast <- gsl.mast.with.tiebreak(mdsfun, X)
  mast.dendogram <- gsl.mast.dendogram(mast)
  obs.density <- rep(0, n)
  for (i in 1:n) obs.density[i] <- mdsfun(i, i)
  pc.values <- gsl.pruning.criteria(mast.dendogram, obs.density)
  cluster.tree <- gsl.compute.cluster.tree(mast.dendogram, pc.values,
                                      pruning.criterion = pruning.criterion,
                                      pruning.threshold = pruning.threshold)
  afo.out <- gsl.assign.fluff.oneshot(X, cluster.tree, mast.dendogram, obs.density)
  return(list(cluster.tree = cluster.tree, leaf.code = afo.out$leaf.code,
              cluster.core = afo.out$cluster.core, mast.dendogram = mast.dendogram,
              pc.values = pc.values, obs.density = obs.density))
}

##-----------------------------------------------------------------

algorithm.illustration.plot <- function(selected.node, X, gsl.cluster.out) {
  n <- nrow(X)
  cluster.tree <- gsl.cluster.out$cluster.tree
  mast.dendogram <- gsl.cluster.out$mast.dendogram
  leaf.code <- gsl.cluster.out$leaf.code
  cluster.core <- gsl.cluster.out$cluster.core
  obs.density <- gsl.cluster.out$obs.density
  leftson <- cluster.tree[selected.node, "leftson"]
  rightson <- cluster.tree[selected.node, "rightson"]
  node.leaf.code <- cluster.tree[selected.node, "leaf.code"]
  leftson.leaf.code <- cluster.tree[leftson, "leaf.code"]
  rightson.leaf.code <- cluster.tree[rightson, "leaf.code"]
  in.node.core <- gsl.select.obs.in.node(leaf.code, cluster.core,
                                            node.leaf.code)$in.core
  in.leftson.core <- gsl.select.obs.in.node(leaf.code, cluster.core,
                                               leftson.leaf.code)$in.core
  in.rightson.core <- gsl.select.obs.in.node(leaf.code, cluster.core,
                                                rightson.leaf.code)$in.core
  in.sons.core <- in.leftson.core | in.rightson.core
  plot(X, type = "n", xlab = "", ylab = "", xaxt = "n", yaxt = "n")
  points(X[in.node.core,], pch = 19, col = "grey")
  points(X[in.sons.core,], pch = 19, col = "black")
  for (i in 1:(n-1)) {
    edge <- mast.dendogram$edges[i,]
    v1 <- edge[1]
    v2 <- edge[2]
    endpoints <- rbind(X[v1,], X[v2,])
    if (in.node.core[v1] & in.node.core[v2]) lines(endpoints, lwd = 2, col = "grey")
    if (in.sons.core[v1] & in.sons.core[v2]) lines(endpoints, lwd = 2, col = "black")
  }
  broken.edge <- edges[cluster.tree[selected.node, "dendogram.node"],]
  v1 <- broken.edge[1]
  v2 <- broken.edge[2]
  endpoints <- rbind(X[v1,], X[v2,])
  lines(endpoints, col = "grey", lwd = 2)
  lines(endpoints, col = "black", lwd = 5, lty = "dotted")
}


#---------------------------------------------------------------------------
# Compute Hubert and Arabie adjusted Rand statistic 
# Equation (5) of Hubert and Arabie (1985)

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

#-----------------------------------------------------------------

choose <- function(n, k, order.matters = F)
{
	if(length(n) < length(k))
		n <- rep(n, length = length(k))
	if(length(n) > length(k))
		k <- rep(k, length = length(n))
	which <- !is.na(n) & !is.na(k)
	if(sum(which) == 0) {
		warning("All missing values in choose")
		return(rep(NA, length(n)))
	}
	if(any(!which)) {
		n <- n[which]
		k <- k[which]
	}
	if(any((n != round(n)) | k != round(k))) {
		warning("n and k are not integers, will be rounded")
		n <- round(n)
		k <- round(k)
	}
	result <- rep(0, length(n))
	zeros <- (k < 0) | (n < k)
	if(any(zeros)) {
		n <- n[!zeros]
		k <- k[!zeros]
	}
	result[!zeros] <- round(exp(lgamma(n + 1) - lgamma(n - k + 1) - {
		if(order.matters)
			0
		else lgamma(k + 1)
	}
	))
	if(any(!which)) {
		y <- rep(NA, length(which))
		y[which] <- result
		y
	}
	else result
}

##-----------------------------------------------------------------

diagonalize.table <- function(tab, shuffle = "columns") {
  if (nrow(tab) > ncol(tab)) shuffle <- "rows"
  if (shuffle == "rows") tab <- t(tab)
  colnames <- dimnames(tab)[[2]]
  for (i in 1:nrow(tab)) {
    jmax <- which.max(tab[i,])
    temp <- tab[,i]
    tab[,i] <- tab[,jmax]
    tab[,jmax] <- temp
    temp <- colnames[i]
    colnames[i] <- colnames[jmax]
    colnames[jmax] <- temp
  }
  dimnames(tab)[[2]] <- colnames
  if (shuffle == "rows") tab <- t(tab)
  return(tab)
}
