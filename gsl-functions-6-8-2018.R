## R code for Generalized Single Linkage Clustering (1-19-09)
##===========================================================
## Changes on May 20, 2009
## Fixed mistakes in gsl.assign.fluff.oneshot and gsl.knn.classifier:
## Did not work for univariate data
##
##
## The code is organized into 3 Chapters.
## - Chapter 1: Numerical routines for GSL
## - Chapter 2: Graphics and diagnostics
## - Chapter 3: Kernel density estimation and cross-validation
##
##
##=================================================================
## Chapter 1: Numerical routines
##=================================================================
## 
## Compute the cluster tree from a minimum density similarity function
## mdsfun. The value of mdsfun(i,i) is the density estimate for
## observation i. The value of mdsfun (i,j) for i != j is the minimum
## of the density estimate along the edge connecting observations i
## and j. The number n of observations is the attribute "n" of mdsfun.
## If computing all the edge weights is too expensive, then we can
## compute them for a smaller test graph, for example the Euclidean
## minimal spanning tree, or the union of several orthogonal msts, or
## the union of the Euclidean mst and a k-near-neighbor graph.
##
## Function returns the cluster tree and leaf.code and cluster.core
## for the observation.


## gsl.cluster <- function(X, mdsfun, pruning.criterion = "size",
##                         pruning.threshold = 0, gsl.cluster.out = NULL,
##                         assign.fluff = T) {
##   n <- attr(mdsfun, "n")
##   if (is.null(gsl.cluster.out)) {
##     mast <- gsl.mast(mdsfun)
##     mast.dendogram <- gsl.mast.dendogram(mast)
##     obs.density <- rep(0, n)
##     for (i in 1:n) obs.density[i] <- mdsfun(i, i)
##     pc.values <- gsl.pruning.criteria(mast.dendogram, obs.density)
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

## -----------------------------------------------------------------
## Modified versions of gsl functions to deal with additional argument
## obs.density.lower.confidence.bounds (5-28-2010)

gsl.cluster <- function(X, mdsfun, pruning.criterion = "size",
                        pruning.threshold = 0, gsl.cluster.out = NULL,
                        assign.fluff = T,
                        obs.density.lower.confidence.bounds = NULL) {
  n <- attr(mdsfun, "n")
  if (is.null(gsl.cluster.out)) {
    mast <- gsl.mast(mdsfun)
    mast.dendogram <- gsl.mast.dendogram(mast)
    obs.density <- rep(0, n)
    for (i in 1:n) obs.density[i] <- mdsfun(i, i)
    pc.values <- gsl.pruning.criteria(mast.dendogram, obs.density,
                                      obs.density.lower.confidence.bounds =
                                      obs.density.lower.confidence.bounds)
  }
  else {
    mast.dendogram <- gsl.cluster.out$mast.dendogram
    pc.values <- gsl.cluster.out$pc.values
    obs.density <- gsl.cluster.out$obs.density
  }
  cluster.tree <- gsl.compute.cluster.tree(mast.dendogram, pc.values,
                                   pruning.criterion = pruning.criterion,
                                   pruning.threshold = pruning.threshold)
  afo.out <- NULL
  if (assign.fluff) afo.out <- gsl.assign.fluff.oneshot(X, cluster.tree,
                                                        mast.dendogram, obs.density)
  
  return(list(cluster.tree = cluster.tree, leaf.code = afo.out$leaf.code,
              cluster.core = afo.out$cluster.core, mast.dendogram = mast.dendogram,
              pc.values = pc.values, obs.density = obs.density))
}

##-----------------------------------------------------------------

gsl.runt.size <- function(gsl.cluster.out) {
  return(rev(sort((gsl.cluster.out$cluster.tree[, "runt.size"]))))
}

gsl.runt.excess.mass <- function(gsl.cluster.out) {
  return(rev(sort((gsl.cluster.out$cluster.tree[, "runt.excess.mass"]))))
}

## 2-26-2011
gsl.runt.pruning.criterion <- function(gsl.cluster.out) {
  return(rev(sort((gsl.cluster.out$cluster.tree[, "runt.pruning.crit"]))))
}

gsl.observation.labels <- function(gsl.cluster.out) {
  return(gsl.cluster.out$leaf.code)
}

gsl.observation.in.core <- function(gsl.cluster.out) {
  return(gsl.cluster.out$cluster.core %% 2 == 1)
}

gsl.cluster.tree <- function(gsl.cluster.out) {
  cluster.tree <- gsl.cluster.out$cluster.tree
  columns <- c("leftson", "rightson", "runt.size", "runt.excess.mass",
               "leaf.code", "level")
  return(cluster.tree[,columns])
}

##-----------------------------------------------------------------

## gsl.mdsmat <- function(X, density, min.density.ut = NULL, kmin = 1,
##                        kmax = nrow(X) - 1, ngrid = 10, interactive = F) {
##   if (is.vector(X)) X <- matrix(X, ncol = 1)
##   n <- nrow(X)
##   p <- ncol(X)
##   all <- (kmin == 1) & (kmax == (n - 1))
##   if(is.null (min.density.ut)) min.density <- matrix(-1, nrow = n, ncol = n)
##   else {
##     min.density <- matrix(0, nrow = n, ncol = n)
##     min.density[upper.tri(min.density, diag = T)] <- min.density.ut
##     min.density <- min.density + t(min.density)
##     diag(min.density) <- diag(min.density) / 2
##     min.density[min.density == min(min.density)] <- -1
##   }
##   if (!all) {
##     Dist <- gsl.interpoint.distance.matrix(X, X)
##     NN <- matrix(0, nrow = n, ncol = n)
##     for (i in 1:n) {
##       NN[i,] <- order(Dist[i,])
##     }
##   }
##   for (i in 1:n) {
##     if (all) closest <- (1:n)[-i] else closest <- NN[i, (kmin + 1):(kmax + 1)]
##     n.closest <- kmax - kmin + 1
##     Xeval <- gsl.find.evaluation.points(X[i,], matrix(X[closest,], nrow = n.closest),
##                                         ngrid)
##     phat.eval <- matrix(density(Xeval), ncol = ngrid, byrow = T)
##     mins <- apply(phat.eval, 1, min) 
##     min.density[i, closest] <- mins
##     min.density[closest, i] <- mins
##     if (interactive) cat(" ", i)
##   }
##   diag(min.density) <- density(X)
##   actual.min <- min(min.density[min.density > 0])
##   min.density[min.density < 0] <- 0.5 * actual.min
##   up <- as.vector(min.density[upper.tri(min.density, diag = T)])
##   return(up)
## }

## 6-19-09: gsl.mdsmat modified to include Euclidean mst edges
gsl.mdsmat <- function(X, density, min.density.ut = NULL, kmin = 1,
                       kmax = nrow(X) - 1, ngrid = 10, interactive = F,
                       include.mst = T) {
  if (is.vector(X)) X <- matrix(X, ncol = 1)
  n <- nrow(X)
  p <- ncol(X)
  all <- (kmin == 1) & (kmax == (n - 1))
  if(is.null (min.density.ut)) min.density <- matrix(-1, nrow = n, ncol = n)
  else {
    min.density <- matrix(0, nrow = n, ncol = n)
    min.density[upper.tri(min.density, diag = T)] <- min.density.ut
    min.density <- min.density + t(min.density)
    diag(min.density) <- diag(min.density) / 2
    min.density[min.density == min(min.density)] <- -1
  }
  if (!all | include.mst) Dist <- gsl.interpoint.distance.matrix(X, X)
  if (!all) {
    ## Dist <- gsl.interpoint.distance.matrix(X, X)
    NN <- matrix(0, nrow = n, ncol = n)
    for (i in 1:n) {
      NN[i,] <- order(Dist[i,])
    }
  }
  if (!all & include.mst) {
    mst.mat.out <- gsl.mst.mat(Dist)
    mst.edges <- mst.mat.out$edges
    for (i in 1:(n-1)) {
      v1 <- mst.edges[i, 1]
      v2 <- mst.edges[i, 2]
      Xeval <- gsl.find.evaluation.points(X[v1,], matrix(X[v2,], nrow = 1),
                                          ngrid)
      phat.eval <- density(Xeval)
      min.density[v1, v2] <- min(phat.eval)
      min.density[v2, v1] <- min(phat.eval)
    }
  }
  for (i in 1:n) {
    if (all) closest <- (1:n)[-i] else closest <- NN[i, (kmin + 1):(kmax + 1)]
    n.closest <- kmax - kmin + 1
    Xeval <- gsl.find.evaluation.points(X[i,], matrix(X[closest,], nrow = n.closest),
                                        ngrid)
    phat.eval <- matrix(density(Xeval), ncol = ngrid, byrow = T)
    mins <- apply(phat.eval, 1, min) 
    min.density[i, closest] <- mins
    min.density[closest, i] <- mins
    if (interactive) cat(" ", i)
  }
  diag(min.density) <- density(X)
  actual.min <- min(min.density[min.density > 0])
  min.density[min.density < 0] <- 0.5 * actual.min
  up <- as.vector(min.density[upper.tri(min.density, diag = T)])
  return(up)
}

gsl.repeat.rows <- function(X, nrep) {
  n <- nrow(X)
  p <- ncol(X)
  vX <- as.vector(X)
  vX.rep <- rep(vX, rep(nrep, n*p))
  R <- matrix(vX.rep, ncol = p)
  return(R)
}

gsl.find.evaluation.points <- function(Xi, X.closest, ngrid) {
  n.closest <- nrow(X.closest)
  p <- length(Xi)
  tgrid <- seq(0, 1, length = ngrid)
  Xi.mat <- gsl.repeat.rows(matrix(Xi, nrow = 1), ngrid * n.closest)
  X.closest.mat <- gsl.repeat.rows(X.closest, ngrid)
  t <- rep(seq(0, 1, length = ngrid), n.closest)
  E <- (1 - t) * Xi.mat + t * X.closest.mat
  return(E)
}

## -----------------------------------------------------------------
## Make minimum density similarity function for the 1nn density estimate.
## Note: Generalized single linkage clustering  with this similarity matrix
## the same as regular single linkage clustering
## Using 1.e50 in place of Inf is a hack, but it works

gsl.make.nn.mdsfun <- function(X, on.the.fly = T) {
  if(on.the.fly) {
    X <- as.matrix(X)
    m <- ncol(X)
    mdsfun <- function(i, js) {
      Y <- matrix(X[js,], ncol = m)
      Z <- matrix(X[i,], nrow = length(js), ncol = m, byrow = T)
      d <- 1/apply((Y - Z)^2, 1, sum)
      d[d == Inf] <- 1.e50
      return(d)
    }
  }
  else {
    D <- 1/as.matrix(dist(X))
    D[D == Inf] <- 1.e50
    mdsfun <- function(i, js) {
      return(D[i, js])
    }
  }
  attr(mdsfun, "n") <- nrow(X)
  return(mdsfun)
}

##-----------------------------------------------------------------
## Make mdsfun from (upper triangle of) minimum density similarity
## matrix produced by "gsl.mdsmat"

gsl.make.mdsfun.from.mdsmat <- function(ut.mdsmat) {
  k <- length(ut.mdsmat)
  n <- (sqrt(1 + 8*k) - 1) / 2
  mdsmat <- matrix(0, nrow = n, ncol = n)
  mdsmat[upper.tri(mdsmat, diag = T)] <- ut.mdsmat
  mdsmat <- mdsmat + t(mdsmat)
  diag(mdsmat) <- diag(mdsmat) / 2
  mdsfun <- function(i, j) {
    return(mdsmat[i,j])
  }
  attr(mdsfun, "n") <- n
  return(mdsfun)
}

##-----------------------------------------------------------------

gsl.make.kernel.mdsfun <- function(X, bandwidth, kmax = 20, ngrid = 10,
                                   interactive = F) {
  density <- make.gaussian.kernel.density.estimate(X, bandwidth)
  mdsmat <- gsl.mdsmat(X, density, kmin = 1, kmax = kmax, ngrid = ngrid,
                       interactive = interactive)
  mdsfun <- gsl.make.mdsfun.from.mdsmat(mdsmat)
  return(mdsfun)
}

## gsl.make.kernel.mdsfun <- function(X, kmax = 20, ngrid = 10, interactive = F) {
##   cv.search.out <- cv.search(X)
##   h <- cv.search.out$opt.smopar
##   if (is.na(h)) stop("\nError in gsl.make.kernel.mdsfun: Cross-validation failed")
##   density <- make.gaussian.kernel.density.estimate(X, h)
##   kmax <- 20
##   kernel.mdsfun <- gsl.make.kernel.mdsfun(X, density, kmax = kmax, interactive = T)
##   return(kernel.mdsfun)
## }

##=================================================================
## Computes a minimal spanning tree from a distance matrix. Produces a
## 2 x (n-1) matrix of edges, and the corresponding edge weights.
## This version does not break ties using Euclidean distance.

gsl.mst.mat <- function(D) {
  n <- nrow(D)
  edges <- matrix(0, nrow = n-1, ncol = 2)
  edge.weights <- rep(0, n-1)
  ##D[diag(D)] <- max(D) + 1
  out.points <- 2:n
  closest.in.points <- rep(1, n-1)
  dist.to.closest.in.point <- D[1, (2:n)]
  for (i in (n-1):2) {
    o <- order(dist.to.closest.in.point)
    new.in.point <- out.points[o[1]]
    edges[i, 1] <- closest.in.points[o[1]]
    edges[i, 2] <- new.in.point
    edge.weights[i] <- dist.to.closest.in.point[o[1]]
    out.points[o[1]] <- out.points[i]
    closest.in.points[o[1]] <- closest.in.points[i]
    dist.to.closest.in.point[o[1]] <- dist.to.closest.in.point[i]
    out.points <- out.points[1:(i-1)]
    closest.in.points <- closest.in.points[1:(i-1)]
    dist.to.closest.in.point <- dist.to.closest.in.point[1:(i-1)]
    dist.to.new.in.point <- D[new.in.point, out.points]
    new.in.point.closer <- (dist.to.new.in.point < dist.to.closest.in.point)
    closest.in.points[new.in.point.closer] <- new.in.point
    dist.to.closest.in.point[new.in.point.closer] <-
      dist.to.new.in.point[new.in.point.closer]
  }
  edges[1, 1] <- out.points[1]
  edges[1, 2] <- closest.in.points[1]
  edge.weights[1] <- dist.to.closest.in.point[1]
  return(list(edges = edges, edge.weights = edge.weights))
}

##-----------------------------------------------------------------
## A more flexible version of the mst function. Distances are supplied
## by a function D. If we are computing the Euclidean mst, then distances
## can be computed on the fly. The size n of the data set is the attribute
## "n" of Dfun

gsl.mst <- function(Dfun) {
  n <- attr(Dfun, "n")
  edges <- matrix(0, nrow = n-1, ncol = 2)
  edge.weights <- rep(0, n-1)
  out.points <- 2:n
  closest.in.points <- rep(1, n-1)
  dist.to.closest.in.point <- Dfun(1, (2:n))
  for (i in (n-1):2) {
    o <- order(dist.to.closest.in.point)
    new.in.point <- out.points[o[1]]
    edges[i, 1] <- closest.in.points[o[1]]
    edges[i, 2] <- new.in.point
    edge.weights[i] <- dist.to.closest.in.point[o[1]]
    out.points[o[1]] <- out.points[i]
    closest.in.points[o[1]] <- closest.in.points[i]
    dist.to.closest.in.point[o[1]] <- dist.to.closest.in.point[i]
    out.points <- out.points[1:(i-1)]
    closest.in.points <- closest.in.points[1:(i-1)]
    dist.to.closest.in.point <- dist.to.closest.in.point[1:(i-1)]
    dist.to.new.in.point <- Dfun(new.in.point, out.points)
    new.in.point.closer <- (dist.to.new.in.point < dist.to.closest.in.point)
    closest.in.points[new.in.point.closer] <- new.in.point
    dist.to.closest.in.point[new.in.point.closer] <-
      dist.to.new.in.point[new.in.point.closer]
    ## cat(" ", i)
  }
  edges[1, 1] <- out.points[1]
  edges[1, 2] <- closest.in.points[1]
  edge.weights[1] <- dist.to.closest.in.point[1]
  return(list(edges = edges, edge.weights = edge.weights))
}

##-----------------------------------------------------------------
## Compute the maximal spanning tree for a similarity function

gsl.mast <- function(Sfun) {
  Dfun <- function(...) {
    return(- Sfun(...))
  }
  attr(Dfun, "n") <- attr(Sfun, "n")
  mst <- gsl.mst(Dfun)
  mast <- mst
  mast$edge.weights <- -mst$edge.weights
  return(mast)
}

##-----------------------------------------------------------------
## Compute mst dendogram. Root node corresponds to longest edge.  Use
## same algorithm as used in C code to compute runt sizes.
## 

gsl.mst.dendogram <- function(gsl.mst.out) {
  mst.edges <- gsl.mst.out$edges
  mst.edge.weights <- gsl.mst.out$edge.weights
  n <- nrow(mst.edges) + 1
  weight.order <- order(mst.edge.weights)
  group.members <- vector("list", n)
  for (i in 1:n) group.members[[i]] <- i
  group.size <- rep(1, n)
  group.index <- 1:n
  merge <- matrix(0, nrow = n-1, ncol = 2)
  edges <- matrix(0, nrow = n-1, ncol = 2)
  height <- rep(0, n-1)
  ## runt.size <- rep(0, n-1)
  group.row <- rep(0, n)
  for (i in 1:(n-1)) {
    e <- mst.edges[weight.order[i],]
    v1 <- e[1]
    v2 <- e[2]
    edges[i, 1] <- v1
    edges[i, 2] <- v2
    height[i] <- mst.edge.weights[weight.order[i]]
    gi1 <- group.index[v1]
    gi2 <- group.index[v2]
    s1 <- group.size[gi1]
    s2 <- group.size[gi2]
    ## runt.size[i] <- min(s1, s2)
    if (s1 == 1) merge[i, 1] <- -group.members[[gi1]]
    else merge[i, 1] <- group.row[gi1]
    if (s2 == 1) merge[i, 2] <- -group.members[[gi2]]
    else merge[i, 2] <- group.row[gi2]
    group.members[[gi1]] <- c(group.members[[gi1]], group.members[[gi2]])
    group.size[gi1] <- s1 + s2
    group.index[group.index == gi2] <- gi1
    group.row[gi1] <- i
  }
  return(list(merge = merge, height = height, edges = edges))
}

##-----------------------------------------------------------------
## Compute mast dendogram from mast. Root node corresponds to shortest
## edge

gsl.mast.dendogram <- function(gsl.mast.out) {
  gsl.mast.out$edge.weights <- - gsl.mast.out$edge.weights
  dend <- gsl.mst.dendogram(gsl.mast.out)
  dend$height <- - dend$height
  return(dend)
}

##-----------------------------------------------------------------
## Compute pruning criteria for each mast dendogram node. There are
## five of them:
## (a) runt size
## (b) runt excess mass
## (c) runt rise = runt (maximum density for leaves of (left resp right)
##     subtree - height of dendogram node
## (d) runt bootstrap vertex rise = runt (maximum lower confidence bound for
##     of leaves of subtree - level of dendogram node). Here the mast 
##     and the dendogram are computed using the upper confidence bounds 
##     for edge and vertex weights. Need to provide lower confidence bounds
##     for observation densities
## (e) runt.bootstrap.edge.rise = similar to (d), but look at edges in
##     subtree instead of leaves. Need to provide lower confidence limits
##     for mast edge weights.
##     
## For each interior dendogram node find the leaves in the left
## right subtrees for which the density is > the height of the
## dendogram node. These are the left and right sizes. The corresponding
## excess masses can be computed as 
## sum(1 - height.of.node / obs.density)/nobs, where the sum is over the
## leaves with density > height of dendogram node
##
## Code below does not implement pruning criterion (e)

gsl.pruning.criteria <- function(mast.dendogram, obs.density,
                                 obs.density.lower.confidence.bounds = NULL) {
  merge <- mast.dendogram$merge
  height <- mast.dendogram$height
  left.size <- rep(0, nrow(merge))
  right.size <- rep(0, nrow(merge))
  left.excess.mass <- rep(0, nrow(merge))
  right.excess.mass <- rep(0, nrow(merge))
  left.rise <- rep(0, nrow(merge))
  right.rise <- rep(0, nrow(merge))
  left.bootstrap.rise <- rep(0, nrow(merge))
  right.bootstrap.rise <- rep(0, nrow(merge))
  pc.int <- function(inode) {
    ileft <- merge[inode, 1]
    iright <- merge[inode, 2]
    if (ileft < 0) left.leaves <- -ileft
    if (ileft > 0) left.leaves <- pc.int(ileft)
    if (iright < 0) right.leaves <- -iright
    if (iright > 0) right.leaves <- pc.int(iright)
    left.high.density.cluster <- left.leaves[obs.density[left.leaves] >=
                                             height[inode]]
    right.high.density.cluster <- right.leaves[obs.density[right.leaves] >=
                                               height[inode]]
    ## if (!is.null(left.high.density.cluster)) {
    if (length(left.high.density.cluster) == 0) browser()
    if (length(left.high.density.cluster) > 0) {
      left.size[inode] <<- length(left.high.density.cluster)
      left.excess.mass[inode] <<- sum(1 - height[inode] /
                                      obs.density[left.high.density.cluster])/
                                        (nrow(merge) + 1)
      left.max <- max(obs.density[left.high.density.cluster])
      if (left.max == Inf) left.rise[inode] <<- Inf
      else left.rise[inode] <<- left.max - height[inode]
      ##left.rise[inode] <<- max(obs.density[left.high.density.cluster]) -
      ##  height[inode]
      ## if (left.rise[inode] == -Inf) browser()
      ##if (left.rise[inode] == Inf) browser()
      if (!is.null(obs.density.lower.confidence.bounds)) 
        left.bootstrap.rise[inode] <<-
          max(obs.density.lower.confidence.bounds[left.high.density.cluster]) -
            height[inode]
    }
    ## if (!is.null(right.high.density.cluster)) {    ## Changed 6-6-2018
    if (length(right.high.density.cluster) == 0) browser()
    if (length(right.high.density.cluster) > 0) {
      right.size[inode] <<- length(right.high.density.cluster)
      right.excess.mass[inode] <<- sum(1 - height[inode] /
                                       obs.density[right.high.density.cluster])/
                                         (nrow(merge) + 1)
      right.max <- max(obs.density[right.high.density.cluster])
      if (right.max == -Inf) browser()  ## 6-6-2018
      if (right.max == Inf) right.rise[inode] <<- Inf
      else right.rise[inode] <<- right.max - height[inode]
      ##right.rise[inode] <<- max(obs.density[right.high.density.cluster]) -
      ##  height[inode]
      ## if (right.rise[inode] == -Inf) browser()  
      ## if (right.rise[inode] == Inf) browser()
      if (!is.null(obs.density.lower.confidence.bounds))
        right.bootstrap.rise[inode] <<-
          max(obs.density.lower.confidence.bounds[right.high.density.cluster]) -
            height[inode]
    }
    return(c(left.leaves, right.leaves))
  }
  pc.int(nrow(merge))
  res <- cbind(left.size, right.size, left.excess.mass, right.excess.mass,
               left.rise, right.rise, left.bootstrap.rise, right.bootstrap.rise)
  colnames(res) <- c("left.size", "right.size", "left.excess.mass", "right.excess.mass",
                     "left.rise", "right.rise", "left.bootstrap.rise",
                     "right.bootstrap.rise")
  return(res)
}

##-----------------------------------------------------------------
## This is the key routine. It constructs the (pruned) cluster tree from
## the dendogram. It does not compute obs.leaf.code and obs.cluster.core.
## The arrangement of columns in the output matrix is to preserve
## compatibility with earlier versions of the code

gsl.compute.cluster.tree <- function(dendogram, gsl.pc.out, pruning.criterion = "size",
                             pruning.threshold = 0) {
  runt.size <- pmin(gsl.pc.out[,"left.size"], gsl.pc.out[,"right.size"])
  runt.excess.mass <- pmin(gsl.pc.out[,"left.excess.mass"],
                           gsl.pc.out[,"right.excess.mass"])
  runt.rise <- pmin(gsl.pc.out[,"left.rise"], gsl.pc.out[,"right.rise"])
  runt.bootstrap.rise <- pmin(gsl.pc.out[,"left.bootstrap.rise"],
                              gsl.pc.out[,"right.bootstrap.rise"])
  if (pruning.criterion == "size") runt.crit <- runt.size
  if (pruning.criterion == "excess.mass") runt.crit <- runt.excess.mass
  if (pruning.criterion == "rise") runt.crit <- runt.rise
  if (pruning.criterion == "bootstrap.rise") runt.crit <- runt.bootstrap.rise
  ##browser()
  merge <- dendogram$merge
  level <- dendogram$height
  nobs <- nrow(merge) + 1

  cluster.tree <- matrix(0, nrow = 2 * nobs - 1, ncol = 14)
  colnames(cluster.tree) <- c("leftson", "rightson", "runt.size", "runt.excess.mass",
                              "leaf.code", "level", "dendogram.node", "size",
                              "excess.mass", "rise", "bootstrap.rise",
                              "runt.rise", "runt.bootstrap.rise", "runt.pruning.crit")
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
      cluster.tree[parent, "runt.pruning.crit"] <<- runt.crit[large.runt.node]
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
    colnames(ct.out) <- c("leftson", "rightson", "runt.size", "runt.excess.mass",
                          "leaf.code", "level", "dendogram.node", "size",
                          "excess.mass", "rise", "bootstrap.rise",
                          "runt.rise", "runt.bootstrap.rise", "runt.pruning.crit")
    return(ct.out)
  }
  return(cluster.tree[1:nnodes,])
}

##-----------------------------------------------------------------
## One shot method for assigning fluff
## (1) Find observations in node using leaf code. (2) Find observations in
## cores of daughters. Set their leaf codes to the leaf codes of the daughters
## and set their cluster core to twice their previous value
## + 1. (3) Find the observations that are in node but not in the cores of the
## daughters. Set their leaf codes to the leaf codes of the daughters and their
## cluster core valaues to 2* the previous value. (4) Recurse
##
## Warning: works only if depth of cluster tree is < 63

gsl.assign.fluff.oneshot <- function(X, cluster.tree, dendogram, obs.density) {
  ## browser()
  merge <- dendogram$merge
  n <- nrow(merge) + 1
  m <- ncol(X)
  obs.leaf.code <- rep(1, n)
  obs.in.core <- rep(1, n)
  collect.obs.in.core <- function(dendogram.node, level) {
    is.in.core <- rep(F, n)
    collect.descendents.of.dendogram.node <- function(dendogram.node) {
      ## cat(" ", dendogram.node)
      if (dendogram.node < 0) {
        if (obs.density[-dendogram.node] > level)
          is.in.core[-dendogram.node] <<- T
        return()
      }
      leftson <- merge[dendogram.node, 1]
      rightson <- merge[dendogram.node, 2]
      collect.descendents.of.dendogram.node(leftson)
      collect.descendents.of.dendogram.node(rightson)
    }
    collect.descendents.of.dendogram.node(dendogram.node)
    return(is.in.core)
  }
  preserve.ncol <- function(x, m) {
    if (is.matrix(x)) return(x)
    else return(matrix(x, ncol = m, byrow = T))
  }
  afo.internal <- function(node.index) {
    ##browser()
    if (cluster.tree[node.index, "leftson"] == 0) return()
    dendogram.node <- cluster.tree[node.index, "dendogram.node"]
    node.leaf.code <- cluster.tree[node.index, "leaf.code"]
    leftson <- cluster.tree[node.index, "leftson"]
    rightson <- cluster.tree[node.index, "rightson"]
    ## cat("\n", node.index, leftson, rightson, "\n")
    leftson.leaf.code <- cluster.tree[leftson, "leaf.code"]
    rightson.leaf.code <- cluster.tree[rightson, "leaf.code"]
    in.node <- obs.leaf.code == node.leaf.code
    ## browser()
    in.left.core <- collect.obs.in.core(dendogram$merge[dendogram.node, 1],
                                        dendogram$height[dendogram.node])
    in.right.core <- collect.obs.in.core(dendogram$merge[dendogram.node, 2],
                                         dendogram$height[dendogram.node])
    obs.leaf.code[in.left.core] <<- leftson.leaf.code
    obs.in.core[in.left.core] <<- 2 * obs.in.core[in.left.core] + 1
    obs.leaf.code[in.right.core] <<- rightson.leaf.code
    obs.in.core[in.right.core] <<- 2 * obs.in.core[in.right.core] + 1
    ## browser()
    if (sum(in.left.core) + sum(in.right.core) < sum(in.node)) {
      X.train <- rbind(preserve.ncol(X[in.left.core,], m),
                       preserve.ncol(X[in.right.core,], m))
      y.train <- c(rep(1, sum(in.left.core)), rep(2, sum(in.right.core)))
      in.fluff <- in.node & !(in.left.core | in.right.core)
      obs.in.core[in.fluff] <<- 2 * obs.in.core[in.fluff]
      X.test <- preserve.ncol(X[in.fluff,], m)
      y.test <- as.vector(gsl.knn.classifier(X.train, y.train, X.test, kk = 1)$ypred)
      obs.leaf.code[in.fluff][y.test == 1] <<- leftson.leaf.code
      obs.leaf.code[in.fluff][y.test == 2] <<- rightson.leaf.code
    }
    afo.internal(leftson)
    afo.internal(rightson)
  }
  afo.internal(1)
  return(list(leaf.code = obs.leaf.code, cluster.core = obs.in.core))
}

## gsl.assign.fluff.oneshot <- function(X, cluster.tree, dendogram, obs.density, assign.flu) {
##   ## browser()
##   merge <- dendogram$merge
##   n <- nrow(merge) + 1
##   m <- ncol(X)
##   obs.leaf.code <- rep(1, n)
##   obs.in.core <- rep(1, n)
##   collect.obs.in.core <- function(dendogram.node, level) {
##     is.in.core <- rep(F, n)
##     collect.descendents.of.dendogram.node <- function(dendogram.node) {
##       ## cat(" ", dendogram.node)
##       if (dendogram.node < 0) {
##         if (obs.density[-dendogram.node] > level)
##           is.in.core[-dendogram.node] <<- T
##         return()
##       }
##       leftson <- merge[dendogram.node, 1]
##       rightson <- merge[dendogram.node, 2]
##       collect.descendents.of.dendogram.node(leftson)
##       collect.descendents.of.dendogram.node(rightson)
##     }
##     collect.descendents.of.dendogram.node(dendogram.node)
##     return(is.in.core)
##   }
##   preserve.ncol <- function(x, m) {
##     if (is.matrix(x)) return(x)
##     else return(matrix(x, ncol = m, byrow = T))
##   }
##   afo.internal <- function(node.index) {
##     ##browser()
##     if (cluster.tree[node.index, "leftson"] == 0) return()
##     dendogram.node <- cluster.tree[node.index, "dendogram.node"]
##     node.leaf.code <- cluster.tree[node.index, "leaf.code"]
##     leftson <- cluster.tree[node.index, "leftson"]
##     rightson <- cluster.tree[node.index, "rightson"]
##     ## cat("\n", node.index, leftson, rightson, "\n")
##     leftson.leaf.code <- cluster.tree[leftson, "leaf.code"]
##     rightson.leaf.code <- cluster.tree[rightson, "leaf.code"]
##     in.node <- obs.leaf.code == node.leaf.code
##     ## browser()
##     in.left.core <- collect.obs.in.core(dendogram$merge[dendogram.node, 1],
##                                         dendogram$height[dendogram.node])
##     in.right.core <- collect.obs.in.core(dendogram$merge[dendogram.node, 2],
##                                          dendogram$height[dendogram.node])
##     obs.leaf.code[in.left.core] <<- leftson.leaf.code
##     obs.in.core[in.left.core] <<- 2 * obs.in.core[in.left.core] + 1
##     obs.leaf.code[in.right.core] <<- rightson.leaf.code
##     obs.in.core[in.right.core] <<- 2 * obs.in.core[in.right.core] + 1
##     in.fluff <- in.node & !(in.left.core | in.right.core)
##     obs.in.core[in.fluff] <<- 2 * obs.in.core[in.fluff]

##     if (assign.fluff) {
##       if (sum(in.left.core) + sum(in.right.core) < sum(in.node)) {
##         X.train <- rbind(preserve.ncol(X[in.left.core,], m),
##                          preserve.ncol(X[in.right.core,], m))
##         y.train <- c(rep(1, sum(in.left.core)), rep(2, sum(in.right.core)))
##         ## in.fluff <- in.node & !(in.left.core | in.right.core)
##         ## obs.in.core[in.fluff] <<- 2 * obs.in.core[in.fluff]
##         X.test <- preserve.ncol(X[in.fluff,], m)
##         y.test <- as.vector(gsl.knn.classifier(X.train, y.train, X.test, kk = 1)$ypred)
##         obs.leaf.code[in.fluff][y.test == 1] <<- leftson.leaf.code
##         obs.leaf.code[in.fluff][y.test == 2] <<- rightson.leaf.code
##       }
##     }
##     afo.internal(leftson)
##     afo.internal(rightson)
##   }
##   afo.internal(1)
##   return(list(leaf.code = obs.leaf.code, cluster.core = obs.in.core))
## }


## -----------------------------------------------------------------
## Compute dendogram of cluster tree. Difficult issue is at which height
## to put the leaves. For the moment, offer three options:
## - "constant" puts them all at the same height, (1 + hang) * highest
##   split level
## - "mean" puts them at the mean density for the observations in the corresponding
##   high density cluster
## - "max" puts them at the maximum density for the corresponding high density
##   cluster

gsl.cluster.dendogram <- function(gsl.cluster.out, leaf.level.rule = "max", hang = 0.1) {
  cluster.tree <- gsl.cluster.out$cluster.tree
  mast.dendogram <- gsl.cluster.out$mast.dendogram
  obs.density <- gsl.cluster.out$obs.density
  if (sum(obs.density) == Inf) leaf.level.rule <- "constant"
  highest.level <- max(cluster.tree[, "level"])
  offset <- hang * highest.level
  obs.leaf.code <- gsl.cluster.out$leaf.code
  obs.in.core <- gsl.cluster.out$cluster.core
  cluster.dendogram <- cluster.tree
  convert.internal <- function(inode) {
    leftson <- cluster.tree[inode, "leftson"]
    rightson <- cluster.tree[inode, "rightson"]
    if (leftson == 0) {
      node.leaf.code <- cluster.tree[inode, "leaf.code"]
      obs.in.node <- gsl.select.obs.in.node(obs.leaf.code, obs.in.core, node.leaf.code)
      cluster.core.density <- obs.density[obs.in.node$in.core]
      if (leaf.level.rule == "mean") leaf.level <- mean(cluster.core.density)
      if (leaf.level.rule == "max") leaf.level <- max(cluster.core.density)
      if (leaf.level.rule == "constant")
        leaf.level <- (1 + hang) * highest.level
      cluster.dendogram[inode, "level"] <<- leaf.level
    }
    else {
      cluster.dendogram[inode, "level"] <<- cluster.tree[leftson, "level"]
      convert.internal(leftson)
      convert.internal(rightson)
    }
  }
  convert.internal(1)
  return(cluster.dendogram)
}

##-----------------------------------------------------------------

gsl.select.obs.in.node <- function(obs.leaf.code, obs.in.core, node.leaf.code) {
  n <- length(obs.leaf.code)
  node.leaf.code <- rep(node.leaf.code, n)
  delta.floor <- floor(log(obs.leaf.code, 2)) - floor(log(node.leaf.code, 2))
  in.node <- floor(obs.leaf.code / 2^delta.floor) == node.leaf.code
  in.core <- in.node & (floor(obs.in.core / 2^delta.floor) %% 2 == 1)
  return(list(in.node = in.node, in.core = in.core))
}

##-----------------------------------------------------------------

gsl.safe.lsfit <- function(X, y, tol = 1.e-6) {
  n <- nrow(X)
  p <- ncol(X)
  X1 <- cbind(rep(1,n), X)
  X1.svd <- svd(X1)
  U <- X1.svd$u
  V <- X1.svd$v
  d <- X1.svd$d
  dpos <- d[d/d[1] > tol]
  rank <- length(dpos)
  ytilde <- t(U[,1:rank]) %*% y
  z <- ytilde / dpos
  b <- V[,1:rank] %*% z
  return(list(coef = b, residuals = y - X1 %*% b))
}

## -----------------------------------------------------------------
##

gsl.knn.classifier <- function(X.train, y.train, X.test, kk = 1,
                               pi = rep(1/K, K), CV = F, block.size = 200) {
  ## browser()
  n.train <- nrow(X.train)
  n.test <- nrow(X.test)
  n.block <- n.test %/% block.size
  m <- ncol(X.train)
  if (n.block * block.size < n.test) n.block <- n.block + 1
  K <- max(y.train)
  ypred <- matrix(0, nrow = n.test, ncol = length(kk))
  weight <- vector("numeric", K)
  for (i in 1:K) {
    ni <- sum(y.train==i)
    if (ni > 0) {weight[i] <- pi[i] * n.train / ni}
    else {weight[i] <- 1}
  } Browser
  for (block in 1:n.block) {
    i.start <- (block - 1) * block.size + 1
    i.end <- min(block * block.size, n.test)
    ## cat("\ni.start = ", i.start)
    ## if (i.end > i.start) X.test.block <- X.test[i.start:i.end,]
    ## Changed 5-20-09
    if (i.end > i.start) {
        X.test.block <- matrix(X.test[i.start:i.end,], ncol = m)
    } else {
        X.test.block <- matrix(X.test[i.start,], ncol = m)
    }
    ## if (ncol(X.train) != ncol(X.test.block)) browser()
    dist <- gsl.interpoint.distance.matrix(X.train, X.test.block)
    if (CV) {dist[row(dist) == (col(dist) + (i.start - 1))] <- max(dist) + 1}
                                        # cat("Computed distance matrix \n")
    permut <- apply(dist, 2, order)
    ordered.class.ids <- matrix(y.train[permut], nrow = n.train)
    counts <- matrix(0, nrow = K, ncol = i.end - i.start + 1)
    ypred.block <- matrix(0, nrow = i.end - i.start + 1, ncol = length(kk))
                                        # cat("Sorted distances \n")
                                        # browser()
    for(j in 1:length(kk)) {
      k <- kk[j]
      for(i in 1:K){
                                        #browser()
        if(k==1) {counts[i,] <- weight[i] * (ordered.class.ids[1,] == i)}
        else {counts[i,] <- weight[i] * 
                apply((ordered.class.ids[1:k,] == i), 2, sum)}
      }
                                        #browser()
      ypred.block[,j] <- apply(counts, 2, which.max)
      counts <- sweep(counts, 2, apply(counts, 2, sum), "/")
    }
    ##browser()
    ypred[i.start:i.end,] <- ypred.block
  }
  return(list(ypred = ypred))
}
##
## -----------------------------------------------------------------
## Changed 3-20-2011
## Changed default for prior probs to assume no misrepresentation:
## pi[i] = n[i] / n.train

## gsl.knn.classifier <- function(X.train, y.train, X.test, kk = 1,
##                                pi = "empirical", CV = F, block.size = 200) {
##   n.train <- nrow(X.train)
##   n.test <- nrow(X.test)
##   n.block <- n.test %/% block.size
##   if (n.block * block.size < n.test) n.block <- n.block + 1
##   K <- max(y.train)
##   ypred <- matrix(0, nrow = n.test, ncol = length(kk))
##   weight <- vector("numeric", K)
##   for (i in 1:K) {
##     ni <- sum(y.train==i)
##     if (pi == "empirical") weight[i] <- 1 else {
##       if (ni > 0) {weight[i] <- pi[i] * n.train / ni}
##       else {weight[i] <- 1}
##     }
##   }
##   for (block in 1:n.block) {
##     i.start <- (block - 1) * block.size + 1
##     i.end <- min(block * block.size, n.test)
##     ## cat("\ni.start = ", i.start)
##     ## if (i.end > i.start) X.test.block <- X.test[i.start:i.end,]
##     ## Changed 5-20-09
##     if (i.end > i.start) X.test.block <- as.matrix(X.test[i.start:i.end,])
##     else X.test.block <- matrix(X.test[i.start,], ncol = ncol(X.train))
##     dist <- gsl.interpoint.distance.matrix(X.train, X.test.block)
##     if (CV) {dist[row(dist) == (col(dist) + (i.start - 1))] <- max(dist) + 1}
##                                         # cat("Computed distance matrix \n")
##     permut <- apply(dist, 2, order)
##     ordered.class.ids <- matrix(y.train[permut], nrow = n.train)
##     counts <- matrix(0, nrow = K, ncol = i.end - i.start + 1)
##     ypred.block <- matrix(0, nrow = i.end - i.start + 1, ncol = length(kk))
##                                         # cat("Sorted distances \n")
##                                         # browser()
##     for(j in 1:length(kk)) {
##       k <- kk[j]
##       for(i in 1:K){
##                                         #browser()
##         if(k==1) {counts[i,] <- weight[i] * (ordered.class.ids[1,] == i)}
##         else {counts[i,] <- weight[i] * 
##                 apply((ordered.class.ids[1:k,] == i), 2, sum)}
##       }
##                                         #browser()
##       ypred.block[,j] <- apply(counts, 2, which.max)
##       counts <- sweep(counts, 2, apply(counts, 2, sum), "/")
##     }
##     ypred[i.start:i.end,] <- ypred.block
##   }
##   return(list(ypred = ypred))
## }


## 6-8-2018
## This function has numerical problems.
## gsl.interpoint.distance.matrix(X, X)[i] is not the same as
## gsl.interpoint.distance.matrix(X, matrix(X[i, ], nrow = 1)
##
## gsl.interpoint.distance.matrix <- function (X, Y) {
##   p <- ncol(X)
##   nx <- nrow(X)
##   ny <- nrow(Y)
##   norm2.x <- apply(X^2, 1, sum)
##   norm2.y <- apply(Y^2, 1, sum)
##   dist <- matrix(norm2.x, nrow = nx, ncol = ny) +
##     matrix(norm2.y, nrow = nx, ncol = ny, byrow = T) -
##       2 * X %*% t(Y)
##   return(dist)
## }

## The next version does not have the problem but is a lot slower
## gsl.interpoint.distance.matrix <- function (X, Y) {
##   m <- ncol(X)
##   nx <- nrow(X)
##   ny <- nrow(Y)
##   Dist <- matrix(0, nrow = nx, ncol = ny)
##   for (i in 1:nx) {
##     Rowi.rep <- matrix(X[i,], nrow = ny, ncol = m, byrow = T)
##     Diff <- (Rowi.rep - Y)^2
##     Dist[i, ] <- apply(Diff, 1, sum)
##   }
##   return(Dist)
## }

## This version also does not have the problem
gsl.interpoint.distance.matrix <- function(X, Y) {
  XY <- rbind(X, Y)
  XY.dist <- as.matrix(dist(XY))
  row.ind <- 1:nrow(X)
  col.ind <- (nrow(X) + 1):(nrow(X) + nrow(Y))
  Dist <- XY.dist[row.ind, col.ind]
  if (nrow(X) == 1) Dist <- matrix(Dist, nrow = 1)
  if (nrow(Y) == 1) Dist <- matrix(Dist, ncol = 1)
  return(Dist^2)
}











##=================================================================
## Chapter 2: Graphics and Diagnostics
##=================================================================
##
## Draw the cluster tree

gsl.draw.cluster.tree <- function(gsl.cluster.out, draw.cluster.numbers = F,
                          draw.cluster.leaf.codes = F,
                          draw.runt.excess.mass = F,
                          draw.runt.size = F) {
  ## browser()
  dct.out <- gsl.draw.tree(gsl.cluster.out$cluster.tree,
                          draw.cluster.numbers = draw.cluster.numbers,
                          draw.cluster.leaf.codes = draw.cluster.leaf.codes,
                          draw.runt.excess.mass = draw.runt.excess.mass,
                          draw.runt.size = draw.runt.size)
  return(dct.out)
}

##-----------------------------------------------------------------
## Draw the cluster dendogram

gsl.draw.cluster.dendogram <- function(gsl.cluster.out, leaf.level.rule = "max",
                                       hang = 0.1, ...) {
  cd.out <- gsl.cluster.dendogram(gsl.cluster.out,
                                  leaf.level.rule = leaf.level.rule,
                                  hang = hang)
  dct.out <- gsl.draw.tree(cd.out, ...)
  return(dct.out)
}

##-----------------------------------------------------------------
## Gnanadesikan, Kettenring, Landwehr diagnostic plot. Find Fisher direction
## discriminating between daughters of an internal dendogram or cluster tree
## node. Project observation in node on Fisher direction and show histogram.

gsl.gkl.diagnostic.plot <- function(X, gsl.cluster.out, tree.type = "cluster.tree",
                                    n.breaks = 20, ...){
  par(mfrow = c(3,1))
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
                                        #title(sub = "(c)")
  selected.node <- identify(node.xy$x, node.xy$y, labels = rep(" ", length(node.xy$x)),
                            n = 1)
  points(node.xy$x[selected.node], node.xy$y[selected.node], pch = 18, cex = 2, 
         col = "blue1")
  left.son <- cluster.tree[selected.node, "leftson"]
  right.son <- cluster.tree[selected.node, "rightson"]
  if (left.son == 0) {
    cat("\nYou should select an interior node\n")
    return()
  }
  ## browser()
  m <- ncol(X)
  selected.leaf.code <- cluster.tree[selected.node, "leaf.code"]
  left.son.leaf.code <- cluster.tree[left.son, "leaf.code"]
  right.son.leaf.code <- cluster.tree[right.son, "leaf.code"]
  selected <- gsl.select.obs.in.node(obs.leaf.code, obs.cluster.core, selected.leaf.code)
  left <- gsl.select.obs.in.node(obs.leaf.code, obs.cluster.core, left.son.leaf.code)
  right <- gsl.select.obs.in.node(obs.leaf.code, obs.cluster.core, right.son.leaf.code)
  is.left.in.core <- left$in.core
  is.right.in.core <- right$in.core
  is.selected <- selected$in.node
  is.in.core <- selected$in.core
  if (sum(is.left.in.core) == 1) left.X <- matrix(X[is.left.in.core,], nrow = 1) else
    ## left.X <- X[is.left.in.core,]
    left.X <- matrix(X[is.left.in.core,], ncol = m)
  if (sum(is.right.in.core) == 1) right.X <- matrix(X[is.right.in.core,], nrow = 1) else
    ## right.X <- X[is.right.in.core,]
    right.X <- matrix(X[is.right.in.core,], ncol = m)
  X.train <- rbind(left.X, right.X)
  y.train <- c(rep(0, nrow(left.X)), rep(1, nrow(right.X)))
  lsmod <- gsl.safe.lsfit(X.train, y.train)
  yhat.selected <- matrix(X[is.selected,], ncol = m) %*% lsmod$coef[2:(ncol(X)+1)]
  yhat.in.core <- matrix(X[is.in.core,], ncol = m) %*% lsmod$coef[2:(ncol(X)+1)]
  breaks <- seq(from = min(yhat.selected), to = max(yhat.selected),
                length = n.breaks)
  hist(yhat.selected, breaks = breaks,  xaxt = "n",
       xlab = "", col = "blue1", main = "Observations in node")
  hist(yhat.in.core, breaks = breaks,  xaxt = "n",
       xlab = "", col = "blue1", main = "Observations in core")
  par(mfrow = c(1,1))
  return()
}

##-----------------------------------------------------------------

gsl.all.gkl.diagnostic.plots <- function(X, gsl.cluster.out, tree.type = "cluster.tree",
                                         n.breaks = 20,...){
  par(mfrow = c(3,1))
  cluster.tree <- gsl.cluster.out$cluster.tree
  obs.leaf.code <- gsl.cluster.out$leaf.code
  obs.cluster.core <- gsl.cluster.out$cluster.core
  agdp.internal <- function(selected.node) {
    left.son <- cluster.tree[selected.node, "leftson"]
    right.son <- cluster.tree[selected.node, "rightson"]
    if (left.son == 0) return()
    ##browser()
    if (tree.type == "cluster.tree") {
      node.xy <- gsl.draw.cluster.tree(gsl.cluster.out, ...)
      title("Cluster tree")
    }
    else {
      node.xy <- gsl.draw.cluster.dendogram(gsl.cluster.out, ...)
      title("Cluster dendogram")
    }
    points(node.xy$x[selected.node], node.xy$y[selected.node], pch = 18, cex = 2, 
           col = "blue1")
    selected.leaf.code <- cluster.tree[selected.node, "leaf.code"]
    left.son.leaf.code <- cluster.tree[left.son, "leaf.code"]
    right.son.leaf.code <- cluster.tree[right.son, "leaf.code"]
    selected <- gsl.select.obs.in.node(obs.leaf.code, obs.cluster.core, selected.leaf.code)
    left <- gsl.select.obs.in.node(obs.leaf.code, obs.cluster.core, left.son.leaf.code)
    right <- gsl.select.obs.in.node(obs.leaf.code, obs.cluster.core, right.son.leaf.code)
    is.left.in.core <- left$in.core
    is.right.in.core <- right$in.core
    is.selected <- selected$in.node
    is.in.core <- selected$in.core
    if (sum(is.left.in.core) == 1) left.X <- matrix(X[is.left.in.core,], nrow = 1)
    else left.X <- X[is.left.in.core,]
    if (sum(is.right.in.core) == 1) right.X <- matrix(X[is.right.in.core,], nrow = 1)
    else right.X <- X[is.right.in.core,]
    X.train <- rbind(left.X, right.X)
    y.train <- c(rep(0, nrow(left.X)), rep(1, nrow(right.X)))
    lsmod <- gsl.safe.lsfit(X.train, y.train)
    yhat.selected <- X[is.selected,] %*% lsmod$coef[2:(ncol(X)+1)]
    yhat.in.core <- X[is.in.core,] %*% lsmod$coef[2:(ncol(X)+1)]
    breaks <- seq(from = min(yhat.selected), to = max(yhat.selected),
                  length = n.breaks)
    hist(yhat.selected, breaks = breaks,  xaxt = "n",
         xlab = "", col = "blue1", main = "Observations in node")
    hist(yhat.in.core, breaks = breaks,  xaxt = "n",
         xlab = "", col = "blue1", main = "Observations in core")
    bla <- readline("\nHit return to see next plot or any character and return to exit ")
    if (bla != "") stop
    agdp.internal(left.son)
    agdp.internal(right.son)
  }
  agdp.internal(1)
  par(mfrow = c(1,1))
}

##-----------------------------------------------------------------
## 

gsl.two.dim.diagnostic.plot <- function(X, gsl.cluster.out, tree.type = "cluster.tree",
                                        ...){
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
                            n = 1)[1]
  points(node.xy$x[selected.node], node.xy$y[selected.node], pch = 18, cex = 2, 
         col = "black")
  left.son <- cluster.tree[selected.node, "leftson"]
  right.son <- cluster.tree[selected.node, "rightson"]
  if (left.son == 0) {
    cat("\nYou should select an interior node\n")
    return()
  }
  points(node.xy$x[left.son], node.xy$y[left.son], pch = 18, cex = 2, 
         col = "red1")
  points(node.xy$x[right.son], node.xy$y[right.son], pch = 18, cex = 2, 
         col = "blue1")
  left.son.leaf.code <- cluster.tree[left.son, "leaf.code"]
  right.son.leaf.code <- cluster.tree[right.son, "leaf.code"]
  left <- gsl.select.obs.in.node(obs.leaf.code, obs.cluster.core, 
                                 left.son.leaf.code)
  right <- gsl.select.obs.in.node(obs.leaf.code, obs.cluster.core, 
                                  right.son.leaf.code)
  is.left.not.core <- left$in.node & !left$in.core
  is.left.core <- left$in.node & left$in.core
  is.right.not.core <- right$in.node & !right$in.core
  is.right.core <- right$in.node & right$in.core
  plot(X[,1], X[,2], type = "n", xlab = "", ylab = "")
  points(X[left$in.node, 1], X[left$in.node, 2], pch = 20, col = "red1")
  points(X[left$in.core, 1], X[left$in.core, 2], pch = 19, col = "red1")
  points(X[right$in.node, 1], X[right$in.node, 2], pch = 20, col = "blue1")
  points(X[right$in.core, 1], X[right$in.core, 2], pch = 19, col = "blue1")
  points(X[!left$in.node & !right$in.node, 1],
         X[!left$in.node & !right$in.node, 2], pch = 18, col = "black")
  return(selected.node)
}

##-----------------------------------------------------------------

gsl.visually.select.observations.in.node <- function(gsl.cluster.out) {
  cluster.tree <- gsl.cluster.out$cluster.tree
  obs.leaf.code <- gsl.cluster.out$leaf.code
  obs.cluster.core <- gsl.cluster.out$cluster.core
  node.xy <- gsl.draw.cluster.tree(gsl.cluster.out, 1)
  selected.node <- identify(node.xy$x, node.xy$y,
                            labels = rep(" ", length(node.xy$x)), n = 1)
  points(node.xy$x[selected.node], node.xy$y[selected.node], pch = 18, cex = 2, 
         col = "black")
  selected.leaf.code <- cluster.tree[selected.node, "leaf.code"]
  selected <- gsl.select.obs.in.node(obs.leaf.code, obs.cluster.core, selected.leaf.code)
  return(list(in.node = selected$in.node, in.core = selected$in.core))
}

##=================================================================

gsl.subtree.gkl.diagnostic.plot <- function(X, gsl.cluster.out, ...){
  n <- nrow(X)
  m <- ncol(X)
  cluster.tree <- gsl.cluster.out$cluster.tree
  obs.leaf.code <- gsl.cluster.out$leaf.code
  obs.cluster.core <- gsl.cluster.out$cluster.core
  par(mfrow = c(2, 1))
  node.xy <- gsl.draw.cluster.tree(gsl.cluster.out, ...)
  title("Cluster tree")
  selected.node <- identify(node.xy$x, node.xy$y, labels = rep(" ", length(node.xy$x)),
                            n = 1)[1]
  points(node.xy$x[selected.node], node.xy$y[selected.node], pch = 18, cex = 2, 
         col = "blue1")
  left.daughter <- cluster.tree[selected.node, 1]
  if (left.daughter == 0) {
    cat("\nYou should select an interior node\n")
    return()
  }
  find.leaf.inodes.for.subtree <- function(inode, ifs = NULL) {
    leftson <- cluster.tree[inode, 1]
    rightson <- cluster.tree[inode, 2]
    if (leftson == 0) ifs <- c(ifs, inode)
    else {
      ifs <- find.leaf.inodes.for.subtree(leftson, ifs)
      ifs <- find.leaf.inodes.for.subtree(rightson, ifs)
    }
    return(ifs)
  }
  leaf.inodes <- find.leaf.inodes.for.subtree(selected.node) 
  leaf.codes <- cluster.tree[leaf.inodes, "leaf.code"]
  n.leaves <- length(leaf.codes)
  X.in.leaves <- matrix(0, nrow = n, ncol = m)
  y.in.leaves <- rep(0, n)
  core.in.leaves <- rep(0, n)
  n.in.leaves <- 0
  for (i in 1:n.leaves) {
    lc <- leaf.codes[i]
    selected <- gsl.select.obs.in.node(obs.leaf.code, obs.cluster.core, lc)
    in.leaf <- selected$in.node
    in.leaf.core <- selected$in.core
    n.in.leaf <- sum(in.leaf)
    ind.range <- (n.in.leaves + 1):(n.in.leaves + n.in.leaf)
    X.in.leaves[ind.range,] <- X[in.leaf,]
    y.in.leaves[ind.range] <- i
    core.in.leaves[ind.range] <- in.leaf.core[in.leaf]
    n.in.leaves <- n.in.leaves + n.in.leaf
    ## browser()
  }
  X.in.leaves <- X.in.leaves[1:n.in.leaves,]
  y.in.leaves <- y.in.leaves[1:n.in.leaves]
  core.in.leaves <- core.in.leaves[1:n.in.leaves]
  ## browser()
  lda.out <- lda(X.in.leaves[(core.in.leaves == 1), ],
                 y.in.leaves[core.in.leaves == 1])
  if (n.leaves == 2) {
    par(mfrow = c(3, 1))
    node.xy <- gsl.draw.cluster.tree(gsl.cluster.out, ...)
    title("Cluster tree")
    points(node.xy$x[selected.node], node.xy$y[selected.node], pch = 18, cex = 2, 
           col = "blue1")
    X.proj <- X.in.leaves %*% lda.out$scaling
    hist(X.proj, xlim = range(X.proj), nclass = 20,  xaxt = "n", xlab = "",
         col = "blue1", main = "Observations in node")
    X.proj.core <- X.proj[core.in.leaves == 1]
    hist(X.proj.core, xlim = range(X.proj), nclass = 20,  xaxt = "n",
         xlab = "", col = "blue1", main = "Observations in core")
    par(mfrow = c(1, 1))
    return()
  }
  par(mfrow = c(2, 1))
  node.xy <- gsl.draw.cluster.tree(gsl.cluster.out, ...)
  title("Cluster tree")
  points(node.xy$x[selected.node], node.xy$y[selected.node], pch = 18, cex = 2, 
         col = "blue1")
  X.proj <- X.in.leaves %*% lda.out$scaling[,1:2]
  plot(X.proj, col = "grey", pch = 20, xaxt = "n", yaxt = "n")
  points(X.proj[(core.in.leaves == 1),], col = "black", pch = 20)
  readline("\nHit return to see colored plot ")

  node.xy <- gsl.draw.cluster.tree(gsl.cluster.out, ...)
  title("Cluster tree")
  points(node.xy$x[selected.node], node.xy$y[selected.node], pch = 18, cex = 2, 
         col = "black")
  points(node.xy$x[leaf.inodes], node.xy$y[leaf.inodes], pch = 19,
         col = rainbow(n.leaves))
  col <- rainbow(n.leaves)[y.in.leaves]     
  plot(X.proj, pch = ".", col = col, xaxt = "n", yaxt = "n")
  points(X.proj[(core.in.leaves == 1),], col = col[core.in.leaves == 1],
         pch = 20)
  return()
}


##-----------------------------------------------------------------
## Draw possibly annotated cluster tree

gsl.draw.tree <- function(cluster.tree, draw.cluster.numbers = F,
                          draw.cluster.leaf.codes = F,
                          draw.runt.excess.mass = F,
                          draw.runt.size = F) {
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
  
  draw.edges(1)
  return(list(x = node.x, y = node.y))
}


##=================================================================
## Chapter 3: Density estimation and cross-validation
##=================================================================

## Error noticed 7-25-2013
## sphere <- function(X) {
##   Centered <- sweep(X, 2, apply(X, 2, mean))
##   Sigma  <- var(Centered)
##   R <- chol(Sigma)
##   Sphered <- t(solve(t(R), t(X)))
##   return(Sphered)
## }

sphere <- function(X) {
  Centered <- sweep(X, 2, apply(X, 2, mean))
  Sigma  <- var(Centered)
  R <- chol(Sigma)
  Sphered <- t(solve(t(R), t(Centered)))
  return(Sphered)
}

standardize <- function(X) {
  centered <- sweep(X, 2, apply(X, 2, mean))
  sds <- sqrt(apply(centered, 2, var))
  scaled <- sweep(centered, 2, sds, FUN = "/")
  return(scaled)
}

##-----------------------------------------------------------------

make.gaussian.kernel.density.estimate <- function(X.train, h) {
  density.estimate <- function(X.eval) {
    phat <- gaussian.kernel.density.estimate(X.train, X.eval, h, cv = F)
    return(phat)
  }
  return(density.estimate)
}

##-----------------------------------------------------------------
## This estimate does not normalize the density - unnecessary in our application

## gaussian.kernel.density.estimate <- function(X.obs, X.eval, h, cv = F) {
##   K <- function(tsquared, h) {
##     return(exp(-tsquared / (2 * h^2)) / h^p)
##   }
##   if (is.vector(X.obs)) {
##     X.obs <- matrix(X.obs, ncol = 1)
##     X.eval <- matrix(X.eval, ncol = 1)
##   }
##   n.eval <- nrow(X.eval)
##   n.obs <- nrow(X.obs)
##   p <- ncol(X.obs)
##   dens <- rep(0, n.eval)
##   squared.dist <- gsl.interpoint.distance.matrix(X.eval, X.obs) 
##   for (i in 1:n.eval) {
##     kernel.values <- K(squared.dist[i,], h)
##     dens[i] <- sum(kernel.values) / n.obs
##     if (cv) {
##       dens[i] <- (sum(kernel.values) - kernel.values[i]) / (n.obs - 1)
##     }
##   }
##   return(dens)
## }

gaussian.kernel.density.estimate <- function(X.obs, X.eval, h, cv = F) {
  ## K <- function(tsquared, h) {
  ##     return(exp(-tsquared / (2 * h^2)) / h^p)
  ##   }
  K <- function(tsquared, h) {
    return(exp(-tsquared / (2 * h^2)) / (((2 * pi)^(p/2)) * h^p))
  }
  if (is.vector(X.obs)) {
    X.obs <- matrix(X.obs, ncol = 1)
    X.eval <- matrix(X.eval, ncol = 1)
  }
  n.eval <- nrow(X.eval)
  n.obs <- nrow(X.obs)
  p <- ncol(X.obs)
  dens <- rep(0, n.eval)
  squared.dist <- gsl.interpoint.distance.matrix(X.eval, X.obs) 
  for (i in 1:n.eval) {
    kernel.values <- K(squared.dist[i,], h)
    dens[i] <- sum(kernel.values) / n.obs
    if (cv) {
      dens[i] <- (sum(kernel.values) - kernel.values[i]) / (n.obs - 1)
    }
  }
  return(dens)
}

##-----------------------------------------------------------------

gaussian.kernel.least.squares.cv <- function(X, h) {
  if (is.vector(X)) X <- matrix(X, ncol = 1)
  n <- nrow(X)
  p <- ncol(X)
  nh <- length(h)
  cv <- rep(0, nh)
  K <- function(tsquared, h) {
    return(exp(-tsquared / (2 * h^2)) / h^p)
  }
  Kstar <- function(tsquared, h) {
    return(K(tsquared, h*sqrt(2)))
  }
  D <- gsl.interpoint.distance.matrix(X, X)
  for (i in 1:nh) {
    cv[i] <- sum(Kstar(D, h[i]))/n^2 - 2 * sum(K(D, h[i]))/(n * (n-1)) +
      2 * K(0, h[i])/(n-1)
  }
  return(cv)
}

##-----------------------------------------------------------------

golden.section <- function(fun, xleft = 0.05, xright = 1, ngrid = 10, maxit = 8){
  xgrid <- seq(xleft, xright, length = ngrid)
  fgrid <- rep(0, ngrid)
  for (i in 1:ngrid) fgrid[i] <- fun(xgrid[i])
  iopt <- which.min(fgrid)
  if ((iopt == 1)|(iopt == ngrid)) return(-1)
  xleft <- xgrid[iopt-1]; fleft = fgrid[iopt-1]
  xmiddle <- xgrid[iopt]; fmiddle = fgrid[iopt]
  xright <- xgrid[iopt+1]; fright = fgrid[iopt+1]
  ## We have initial bracket
  for (i in 1:maxit) {
    if ((xmiddle - xleft) > (xright - xmiddle)) {
      xnew <- xmiddle - 0.38 * (xmiddle - xleft)
      fnew <- fun(xnew)
      if (fnew < fmiddle) {
        xright <- xmiddle; fright <- fmiddle; xmiddle <- xnew; fmiddle <- fnew
      }
      else {
        xleft <- xnew; fleft = fnew
      }
    }
    else {
      xnew <- xright - 0.38 * (xright - xmiddle)
      fnew <- fun(xnew)
      if (fnew < fmiddle) {
        xleft <- xmiddle; fleft <- fmiddle; xmiddle <- xnew; fmiddle <- fnew
      }
      else{
        xright <- xnew; fright <- fnew
      }
    }
  }
  return(xmiddle)
}

##-----------------------------------------------------------------
## cv.search either takes a vector of k's or a vector of h's as input.
## If it's a vector of k's then it's an adaptive estimate and we won't
## search for h.
## cv.search <- function(X, cv.function = gaussian.kernel.least.squares.cv,
##                       trial.par = c(0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.6, 0.8,
##                                     1.0, 1.5, 2.0),
##                       search.for.h = T) {
cv.search <- function(X, trial.par = c(0.001, 0.01, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4,
                           0.6, 0.8, 1.0, 1.5, 2.0, 4, 8)) {
  cv.function <- gaussian.kernel.least.squares.cv
  search.for.h <- T
  cv <- cv.function(X, trial.par)
  finite <- is.finite(cv)
  opt.smopar <- NA
  finite.par <- NA
  finite.cv <- NA
  if (sum(finite) > 3) {
    finite.cv <- cv[finite]
    finite.par <- trial.par[finite]
    imin <- which.min(finite.cv)
    if((imin != 1) & (imin != length(finite.cv))) {
      opt.smopar <- finite.par[imin]
      if (search.for.h) {
        fun <- function(h) {
          return(cv.function(X, h))
        }
        opt.smopar <- golden.section(fun, xleft = finite.par[imin - 1],
                                     xright = finite.par[imin + 1], maxit = 8)
      }
    }
  }
  return(list(opt.smopar = opt.smopar, smopars = finite.par, cv = finite.cv))
}

##-----------------------------------------------------------------
## Compute Hubert and Arabie adjusted Rand statistic 
## Equation (5) of Hubert and Arabie (1985)

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



