## Diagnostic plots and testing for unimodality (11-2-2010)
## ========================================================

## Requires
## library("MASS")
## library("fdrtool")

## This file contains functions written by Jeremy Tantrum in
## Winter 2003 and functions written by WXS in November 2011.

## Here are Jeremy's functions

## plot.ucdf(x) - plots the cdf of x and the closest unimodal cdf of x.

## plot.silverman(x) - plots the unimodal Gaussian smoother closest to the
##                     x and the closest bimodal Gaussian smoother.

## calcdip(x) - calculates the dip test statistic for x, using the mode found
##              by the closest unimodal Gaussian smoother.

## unisample(cd.out,n) - generates a sample from the unimodal distribution
##                       returned by the output of calcdip.

##=================================================================

## Jeremy's UCDF code (assuming x is ordered
## (1) Find smallest bandwidth that gives unimodal kernel estimate
##     and the observation x_m closest to the corresponding mode
## (2) Find greatest convex minorant of F_n on [x_1...x_m] and smallest
##     convex minorant on [x_m ... x_n + 1]

plot.ucdf <- function(x,plot.it=T,COL1=1,COL2=2,LWD1=3,LWD2=2)
{
  x.cdf <- (1:length(x))/length(x)
  plot(c(min(x),sort(x)),c(0,x.cdf),type='n',xlab="data",ylab="cdf")
  lines(c(par('usr')[1],sort(c(x,x)),par('usr')[2]),sort(c(0,0,x.cdf,x.cdf)),
        lwd=LWD1,col=COL1)
  h.0 <- critwidth(x,4)
  x.f0 <- density(x,width=h.0$low)
  x.mode <- x.f0$x[order(x.f0$y)[length(x.f0$x)]]
  x.mode2 <- x[order(abs(x-x.mode))[1]]
  if(x.mode2 < sort(x)[4])
    x.mode2 <- sort(x)[4]
  if(x.mode2 > sort(x)[length(x)-4])
    x.mode2 <- sort(x)[length(x)-4]
  x.ord <- sort(x)

  x.split <- x.ord < x.mode2
  x.split2 <- x.ord >= x.mode2
  hull1 <- chull(c(x.ord[x.split],x.mode2),c(0,x.cdf[x.split]))
  n.1 <- sum(x.split)+1
  hlist <- c(hull1,hull1,hull1)
  start <- (1:length(hlist))[hlist==n.1][2]
  if(hlist[start+1]==1) hlist <- rev(hlist)
  start <- (1:length(hlist))[hlist==n.1][2]
  x.hull1 <- n.1
  i <- start
  if(findslope(x.ord,x.cdf,1,hlist[i]) > findslope(x.ord,x.cdf,1,hlist[i+1]))
    while(hlist[i] > 1)
      {
        i <- i + 1
        x.hull1 <- c(x.hull1, hlist[i])
      }
  else
    while(hlist[i] > 1)
      {
        i <- i - 1
        x.hull1 <- c(x.hull1, hlist[i])
      }
  x.hull1 <- sort(x.hull1)
  hull2 <- chull(x.ord[x.split2],x.cdf[x.split2])
  n.2 <- sum(x.split2)
  hlist <- c(hull2,hull2,hull2)
  start <- (1:length(hlist))[hlist==n.2][2]
  if(hlist[start+1]==1) hlist <- rev(hlist)
  start <- (1:length(hlist))[hlist==n.2][2]
  x.hull2 <- n.2
  i <- start
  if(findslope(x.ord[x.split2],x.cdf[x.split2],1,hlist[i]) <
     findslope(x.ord[x.split2],x.cdf[x.split2],1,hlist[i+1]))
    while(hlist[i] > 1)
      {
        i <- i +1
        x.hull2 <- c(x.hull2, hlist[i])
      }
  else
    while(hlist[i] > 1)
      {
        i <- i - 1
        x.hull2 <- c(x.hull2, hlist[i])
      }
  x.hull2 <- sort(x.hull2)
  lines(c(x.ord[x.split],x.mode2)[x.hull1],c(0,x.cdf[x.split])[x.hull1],
        col=COL2,lwd=LWD2)
  lines(x.ord[x.split2][x.hull2],x.cdf[x.split2][x.hull2],col=COL2,lwd=LWD2)

  hull.out <- list(hull1=x.hull1, hull2=x.hull2+n.1-1,
                   x = x.ord[sort(unique(c(x.hull1,x.hull2+n.1-1)))],
                   y = x.cdf[sort(unique(c(x.hull1,x.hull2+n.1-1)))])
  delta <- rep(0,length(x))
  for(i in 1:length(x))
    delta[i] <- abs(x.cdf[i] - fofx(sort(x)[i],hull.out))

  return(invisible(list(cdf=x.cdf,unicurve=hull.out, mode=x.mode2,
                        dip = max(delta), delta=delta,
                        dipwhere = order(delta)[length(delta)] )))
}

##-----------------------------------------------------------------

plot.silverman <- function(x,COL1=1,COL2="grey",...)
{
  h.0 <- critwidth(x,4,tol=0.0001)$high
  h.1 <- critwidth2(x,h.0,tol=0.001)$high
  f.c <- density(x,window="g",width=h.0,n=100)
  f.n <- density(x,width=h.1,n=100)
  plot(c(f.c$x,f.n$x),c(f.c$y,f.n$y),type='n',xlab="",ylab="density",...)
  lines(f.n,lwd=3,col=COL1)
  lines(f.c,lwd=3,col=COL2)
  points(x,rep(0,length(x)),pch='|')
}

##-----------------------------------------------------------------

critwidth <- function(g.data,start,tol=0.001,n.points=200)
{
  if(is.unimodal(density(g.data,window="g",width=start,n=n.points)))
    {
      high <- start
      low <- start/2
      while(is.unimodal(density(g.data,window="g",width=low,n=n.points)))
        low <- low/2
    }
  else
    {
      low <- start
      high <- start*2
      while(!is.unimodal(density(g.data,window="g",width=high,n=n.points)))
        high <- high*2
    }
                                        #is.unimodal(low)=F and is.unimodal(high)=T
  while(high-low>tol)
    {
      wdth <- 0.5 * (high+low)
      if(is.unimodal(density(g.data,window="g",width=wdth,n=n.points)))
        high <- wdth
      else
        low <- wdth
    }
  return(list(low=low,high=high))
}

##-----------------------------------------------------------------

critwidth2 <- function(g.data,h.0,tol=0.001,n.points=200)
{
                                        #h.0 is the critical width for a is.unimodal
  start <- h.0 + 2 * tol
  if(is.bimodal(density(g.data,window="g",width=start,n=n.points)))
    {
      high <- start
      low <- start/2
      while(is.bimodal(density(g.data,window="g",width=low,n=n.points)))
        low <- low/2
    }
  else
    {
      low <- start
      high <- start*2
      while(!is.bimodal(density(g.data,window="g",width=high,n=n.points)))
        high <- high*2
    }
  ##is.unimodal(low)=F and is.unimodal(high)=T
  while(high-low>tol)
    {
      wdth <- 0.5 * (high+low)
      if(is.bimodal(density(g.data,window="g",width=wdth,n=n.points)))
        high <- wdth
      else
        low <- wdth
    }
  return(list(low=low,high=high))
}

##-----------------------------------------------------------------

is.unimodal <- function(dens)
{
  ##dens is a list of dens$x and dens$y
  cdf <- cumsum(dens$y)
  n <- length(cdf)
  cdf.diff1 <- cdf[-1] - cdf[-n]
  cdf.diff2 <- cdf.diff1[-1] - cdf.diff1[-(n-1)]
  return(!any(order(-sign(cdf.diff2)) - 1:(n-2) >0))
}

##-----------------------------------------------------------------

is.bimodal <- function(dens)
{
  ##dens is a list of dens$x and dens$y
  cdf <- cumsum(dens$y)
  n <- length(cdf)
  cdf.diff1 <- cdf[-1] - cdf[-n]
  cdf.diff2 <- cdf.diff1[-1] - cdf.diff1[-(n-1)]
  return(sum(sign(cdf.diff2)[-1] - sign(cdf.diff2)[-(n-2)] < 0)<=2)
}

##-----------------------------------------------------------------

calcdip <- function(x,plot.it=T,calc.it=T)
{
  h.0 <- critwidth(x,4)
  x.f0 <- density(x,width=h.0$low)
  x.mode <- x.f0$x[order(x.f0$y)[length(x.f0$x)]]
  x.mode2 <- x[order(abs(x-x.mode))[1]]
  if(x.mode2 < sort(x)[4])
    x.mode2 <- sort(x)[4]
  if(x.mode2 > sort(x)[length(x)-4])
    x.mode2 <- sort(x)[length(x)-4]
  x.cdf <- (1:length(x))/length(x)
  hull.out <- findhulls(sort(x),x.cdf,x.mode2,plot.it=plot.it,xlab="",ylab="CDF")
  delta <- rep(0,length(x))
  if(calc.it)
    for(i in 1:length(x))
      delta[i] <- abs(x.cdf[i] - fofx(sort(x)[i],hull.out))
  return(list(dip = max(delta),unicurve = hull.out))
}

##-----------------------------------------------------------------

unisample <- function(hull.out,size)
{
  n <- length(hull.out$x)
  min.x <- hull.out$x[1] - hull.out$y[1]/findslope(hull.out$x,hull.out$y,1,2)
  probs <- hull.out$y[-1] - c(0,hull.out$y[-c(1,n)])
  where <- sample(1:(n-1),size,replace =T,prob=probs)
  out <- numeric(0)
  x <- c(min.x,hull.out$x[-1])
  for(i in 2:n)
    {
      x.s <- sum(where==i-1)
      if(x.s>0)
        out <- c(out,runif(x.s,x[i-1],x[i]))
    }
  return(out)
}

##-----------------------------------------------------------------

findhulls <- function(x.ord,x.cdf,x.mode,plot.it=T,...)
{
  x.split <- x.ord <= x.mode
  x.split2 <- x.ord >= x.mode
  hull1 <- chull(x.ord[x.split],x.cdf[x.split])
  n.1 <- sum(x.split)
  hlist <- c(hull1,hull1,hull1)
  start <- (1:length(hlist))[hlist==n.1][2]
  if(hlist[start+1]==1) hlist <- rev(hlist)
  start <- (1:length(hlist))[hlist==n.1][2]
  x.hull1 <- n.1
  i <- start
  if(findslope(x.ord,x.cdf,1,hlist[i]) > findslope(x.ord,x.cdf,1,hlist[i+1]))
    while(hlist[i] > 1)
      {
        i <- i + 1
        x.hull1 <- c(x.hull1, hlist[i])
      }
  else
    while(hlist[i] > 1)
      {
        i <- i - 1
        x.hull1 <- c(x.hull1, hlist[i])
      }
  x.hull1 <- sort(x.hull1)  
  hull2 <- chull(x.ord[x.split2],x.cdf[x.split2])
  n.2 <- sum(x.split2)
  hlist <- c(hull2,hull2,hull2)
  start <- (1:length(hlist))[hlist==n.2][2]
  if(hlist[start+1]==1) hlist <- rev(hlist)
  start <- (1:length(hlist))[hlist==n.2][2]
  x.hull2 <- n.2
  i <- start
  if(findslope(x.ord[x.split2],x.cdf[x.split2],1,hlist[i]) <
     findslope(x.ord[x.split2],x.cdf[x.split2],1,hlist[i+1]))
    while(hlist[i] > 1)
      {
        i <- i +1
        x.hull2 <- c(x.hull2, hlist[i])
      }
  else
    while(hlist[i] > 1)
      {
        i <- i - 1
        x.hull2 <- c(x.hull2, hlist[i])
      }
  x.hull2 <- sort(x.hull2)
  if(plot.it)
    {
      plot(x.ord,x.cdf,...)
      lines(x.ord[x.split][x.hull1],x.cdf[x.split][x.hull1])
      lines(x.ord[x.split2][x.hull2],x.cdf[x.split2][x.hull2])
    }
  return(list(hull1=x.hull1, hull2=x.hull2+n.1 -1,
              x = x.ord[sort(unique(c(x.hull1,x.hull2+n.1-1)))],
              y = x.cdf[sort(unique(c(x.hull1,x.hull2+n.1-1)))]))
}

##-----------------------------------------------------------------

findslope <- function(x,y,i,j)  return((y[j] - y[i])/(x[j]-x[i]))

##-----------------------------------------------------------------

fofx <- function(x,hull.out)
{
  n <- length(hull.out$x)
  if(x <= hull.out$x[1])
    return(0)
  if(x >= hull.out$x[n])
    return(1)
  where <- (1:n)[order(c(x,hull.out$x))==1]
  return( (x-hull.out$x[where-1])/(hull.out$x[where]-hull.out$x[where-1]) *
         (hull.out$y[where]-hull.out$y[where-1]) + hull.out$y[where-1])
}

##=================================================================
##=================================================================
## Here are WXS functions, implementing some of the same functionality

birge.find.closest.unimodal.cdf <- function (x, mode = NULL) {
  x <- sort(x)
  n <- length(x)
  if (is.null(mode)) {
    h.start <- 4.0
    h.crit <- critwidth(x, h.start)$high
    density.out <- density(x, width = h.crit)
    dx <- density.out$x
    dy <- density.out$y
    mode <- dx[which(dy == max(dy))[1]]
  }
  i.mode <- (1:n)[order(abs(x - mode))[1]]
  cdf.x <- c(x, x[n] + 1)
  cdf.y <- c(0, (1:n)/n)
  gcm <- gcmlcm(cdf.x[1:i.mode], cdf.y[1:i.mode], type = "gcm")
  lcm <- gcmlcm(cdf.x[i.mode:(n+1)], cdf.y[i.mode:(n+1)], type = "lcm")
  knots.x <- c(gcm$x.knots, lcm$x.knots[2:length(lcm$x.knots)])
  knots.y <- c(gcm$y.knots, lcm$y.knots[2:length(lcm$x.knots)])
  return(list(knots.x = knots.x, knots.y = knots.y, mode = mode))
}

##-----------------------------------------------------------------
## Piecewise linear interpolation

pl.inter <- function(knots.x, knots.y, query.x) {
  nknots <- length(knots.x)
  nquery <- length(query.x)
  query.y <- rep(0, nquery)
  for (i in 1:nquery) {
    x <- query.x[i]
    if (x <= knots.x[1]) {
      query.y[i] <- knots.y[1]
      next
    }
    if (x >= knots.x[nknots]) {
      query.y[i] <- knots.y[nknots]
      next
    }
    j <- ((1:nknots)[x <= knots.x])[1]
    j1 <- j -1
    query.y[i] <- knots.y[j1] + (knots.y[j] - knots.y[j1]) *
      (x - knots.x[j1]) / (knots.x[j] - knots.x[j1])
  }
  return(query.y)
}

##-----------------------------------------------------------------

plot.pl.density <- function(knots.x, knots.y) {
  nknots <- length(knots.x)
  levels <- rep(0, nknots - 1)
  for (i in 1:(nknots - 1)) {
    levels[i] <- (knots.y[i+1] - knots.y[i]) / (knots.x[i+1] - knots.x[i])
  }
  plot(knots.x[1:(nknots-1)], levels, type= "n")
  for (i in 1:(nknots - 1)) {
    lines(knots.x[i:(i+1)], rep(levels[i], 2), lwd = 3)
  }
}

##=================================================================
## Estimate a monotonically decreasing density on [lb, infty]
## Assumes no ties
## Randomly merge inversions

mdd.random <- function(x, lb = 0) {
  merge.blocks <- function(i) {
    blocks[i, "upper.index"] <<- blocks[i+1, "upper.index"]
    blocks[i, "count"] <<- blocks[i, "count"] + blocks[i+1, "count"]
    blocks[i, "length"] <<- blocks[i, "length"] + blocks[i+1, "length"]
    blocks[i, "phat"] <<- blocks[i, "count"] / (n * blocks[i, "length"])
    if((i+2) <= n) 
      for (j in (i+2):n) blocks[(j-1),] <<- blocks[j,]
    nblocks <<- nblocks - 1
  }
  pick.block <- function() {
    candidates <- NULL
    for (i in 1:(nblocks - 1)) 
      if (blocks[i, "phat"] < blocks[(i+1), "phat"])
        candidates <- c(candidates, i)
    if (length(candidates) <= 1) return(candidates)
    else return(sample(candidates, 1))
  }
  x <- sort(x)
  n <- length(x)
  nblocks <- n
  blocks <- matrix(0, nrow = n, ncol = 5)
  colnames(blocks) <- c("lower.index", "upper.index", "count", "length", "phat")
  blocks[,"lower.index"] <- 0:(n-1)
  blocks[,"upper.index"] <- 1:n
  blocks[,"count"] <- rep(1, n)
  blocks[,"length"] <- x - c(lb, x[-n])
  blocks[,"phat"] <- 1 / (n * blocks[,"length"])

  picked.block <- pick.block()
  while(!is.null(picked.block)) {
    merge.blocks(picked.block)
    picked.block <- pick.block()
  }
  return(blocks[1:nblocks,])
}


##=================================================================
## Isotonic regression.
## Assumes that y is sorted oin increasing order of x.
## Assumes no ties in the x-es. Randomly merge inversions

isotonic.regression.random <- function(y) {
  merge.blocks <- function(i) {
    blocks[i, "upper.index"] <<- blocks[i+1, "upper.index"]
    blocks[i, "m"] <<- (blocks[i, "count"] * blocks[i, "m"] +
                        blocks[i+1, "count"] * blocks[i+1, "m"]) /
                          (blocks[i, "count"] + blocks[i+1, "count"])
    blocks[i, "count"] <<- blocks[i, "count"] + blocks[i+1, "count"]
    if((i+2) <= nblocks)
      blocks[(i+1):(nblocks-1),] <<- blocks[(i+2):nblocks,]
    ## for (j in (i+2):nblocks) blocks[(j-1),] <<- blocks[j,]
    nblocks <<- nblocks - 1
  }
  pick.block <- function() {
    candidates <- NULL
    if (nblocks == 1) return(candidates)
    for (i in 1:(nblocks - 1))
      if (blocks[i, "m"] > blocks[(i+1), "m"])
        candidates <- c(candidates, i)
    if (length(candidates) <= 1) return(candidates)
    else return(sample(candidates, 1))
  }
  n <- length(y)
  nblocks <- n
  blocks <- matrix(0, nrow = n, ncol = 4)
  colnames(blocks) <- c("lower.index", "upper.index", "count", "m")
  blocks[,"lower.index"] <- 1:n
  blocks[,"upper.index"] <- 1:n
  blocks[,"count"] <- rep(1, n)
  blocks[,"m"] <- y

  picked.block <- pick.block()
  while(!is.null(picked.block)) {
    merge.blocks(picked.block)
    picked.block <- pick.block()
  }
  yhat <- rep(0, n)
  for (i in 1:nblocks)
    yhat[blocks[i, "lower.index"]:blocks[i, "upper.index"]] <- blocks[i, "m"]
  blocks <- matrix(blocks[1:nblocks,], nrow = nblocks, ncol = 4)
  colnames(blocks) <- c("lower.index", "upper.index", "count", "m")
  return(list(blocks = blocks, yhat = yhat))
}

##-----------------------------------------------------------------
## Assumes y is sorted in increasing order of x

isotonic.regression <- function(y) {
  merge.blocks <- function(i) {
    blocks[i, "upper.index"] <<- blocks[i+1, "upper.index"]
    blocks[i, "m"] <<- (blocks[i, "count"] * blocks[i, "m"] +
                        blocks[i+1, "count"] * blocks[i+1, "m"]) /
                          (blocks[i, "count"] + blocks[i+1, "count"])
    blocks[i, "count"] <<- blocks[i, "count"] + blocks[i+1, "count"]
    if((i+2) <= nblocks) 
      ## for (j in (i+2):nblocks) blocks[(j-1),] <<- blocks[j,]
      blocks[(i+1):(nblocks-1),] <<- blocks[(i+2):nblocks,]
    nblocks <<- nblocks - 1
  }
  
  n <- length(y)
  nblocks <- n
  blocks <- matrix(0, nrow = n, ncol = 4)
  colnames(blocks) <- c("lower.index", "upper.index", "count", "m")
  blocks[,"lower.index"] <- 1:n
  blocks[,"upper.index"] <- 1:n
  blocks[,"count"] <- rep(1, n)
  blocks[,"m"] <- y

  i <- 1
  while (i < nblocks) {
    if (blocks[i, "m"] <= blocks[(i+1), "m"]) {i <- i+1; next}
    merge.blocks(i)
    if (i > 1) {
      for (j in (i-1):1) {
        if (blocks[j, "m"] <= blocks[(j+1), "m"]) break
        merge.blocks(j)
        i <- i-1
      }
    }
  }
  yhat <- rep(0, n)
  for (i in 1:nblocks)
    yhat[blocks[i, "lower.index"]:blocks[i, "upper.index"]] <- blocks[i, "m"]
  blocks <- matrix(blocks[1:nblocks,], nrow = nblocks, ncol = 4)
  colnames(blocks) <- c("lower.index", "upper.index", "count", "m")
  return(list(blocks = blocks, yhat = yhat))
}

##-----------------------------------------------------------------

anti.isotonic.regression <- function(y) {
  ir.out <- isotonic.regression(-y)
  ir.out$blocks[, "m"] <- -ir.out$blocks[, "m"]
  ir.out$yhat <- - ir.out$yhat
  return(ir.out)
}

##-----------------------------------------------------------------
## PAV does not work for unimodal regression.
## Alternative: Pick index k. Find monotonically increasing fit
## on [x[1]...x[k+1]) and monotonically decreasing fit on
## [x[k+1]...x[n]]. Then either x[k] or x[k+1] is a mode. Optimize
## k over 0..n

unimodal.regression <- function(y) {
  n <- length(y)
  rss <- rep(0, (n+1))
  ur.out <- anti.isotonic.regression(y)
  rss.opt <- sum((ur.out$yhat - y)^2)
  rss[1] <- rss.opt
  for (k in 1:(n-1)) {
    ir.out <- isotonic.regression(y[1:k])
    air.out <- anti.isotonic.regression(y[(k+1):n])
    rss.k <- sum((ir.out$yhat - y[1:k])^2) +
             sum((air.out$yhat - y[(k+1):n])^2)
    rss[k+1] <- rss.k
    if (rss.k < rss.opt) {
      rss.opt <- rss.k
      air.out$blocks[, "lower.index"] <- air.out$blocks[, "lower.index"] + k
      air.out$blocks[, "upper.index"] <- air.out$blocks[, "upper.index"] + k
      ur.out <- list(blocks = rbind(ir.out$blocks, air.out$blocks),
                     yhat = c(ir.out$yhat, air.out$yhat))
    }
  }
  ir.out <- isotonic.regression(y)
  rss.k <- sum((ir.out$yhat - y)^2)
  rss[n+1] <- rss.k
  if (rss.k <= rss.opt) ur.out <- ir.out
  return(list(blocks = ur.out$blocks, yhat = ur.out$yhat, rss = rss))
}

##-----------------------------------------------------------------
## Assumes that density = 0 outside of [grid[1], grid[ngrid]

unimodalize.density <- function(density, grid) {
spacing <- grid[2] - grid[1]
ngrid <- length(grid)
phat <- density(grid)
mass <- spacing * sum(phat[1:(ngrid - 1)])
phat <- phat / mass
ur.out <- unimodal.regression(phat[1:(ngrid - 1)])
phat <- c(ur.out$yhat, ur.out$yhat[ngrid - 1])
return(list(blocks = ur.out$blocks, phat = phat, grid = grid))
}

##-----------------------------------------------------------------
## cv.search assumes that data have sd = 1

## unimodal.gkde.fit <- function(x, ngrid = 100) {
## cv.out <- cv.search(x)
## h <- cv.out$opt.smopar
## if (is.na(h)) stop("unimodal.gkde: CV failed")
## density <- make.gaussian.kernel.density.estimate(x, h)
## grid <- seq(from = min(x) - 2*h, to = max(x) + 2*h, length.out = ngrid)
## return(unimodalize.density(density, grid))
## }

unimodal.gkde.fit <- function(x, ngrid = 100) {
cv.out <- cv.search(x / sd(x))
h <- cv.out$opt.smopar
if (is.na(h)) stop("unimodal.gkde: CV failed")
h <- h * sd(x)
density <- make.gaussian.kernel.density.estimate(x, h)
grid <- seq(from = min(x) - 2*h, to = max(x) + 2*h, length.out = ngrid)
return(unimodalize.density(density, grid))
}



##-----------------------------------------------------------------

unimodal.gkde.sample <- function(unimodal.gkde.fit, n) {
grid <- unimodal.gkde.fit$grid
ngrid1 <- length(grid) - 1
phat <- unimodal.gkde.fit$phat
spacing <- grid[2] - grid[1]
u <- runif(n, min = 0, max = spacing)
igrid <- sample(1:ngrid1, n, replace = T, prob = spacing * phat[1:ngrid1])
return(u + grid[igrid])
}

##-----------------------------------------------------------------
## Should generate unimodal density with same mean and covariance as
## the data.

fit.unimodals <- function(X) {
  n <- nrow(X)
  m <- ncol(X)
  unimodal.fits <- vector("list", m)
  for (i in 1:m) {
    ugf.out <- unimodal.gkde.fit(X[,i], ngrid = 100)
    unimodal.fits[[i]] <- ugf.out
    cat(" ",i)
  }
  return(unimodal.fits)
}

##-----------------------------------------------------------------
## X is supposed to be in eigen-coordinates, so that the uniform reference
## has the same "shape" as the data

runt.size.reference <- function(X, nrep, nrs = 20, reference.type =
                                "unimodal") {
  n <- nrow(X)
  m <- ncol(X)
  rs.unimodal <- matrix(nrow = nrep, ncol = nrs)
  X.unimodal <- matrix(0, nrow = n, ncol = m)
  if (reference.type == "unimodal") unimodal.fits <- fit.unimodals(X)
  for (i in 1:nrep) {
    for (j in 1:m) {
      if (reference.type == "uniform")
        X.unimodal[,j] <- sd(X[,j]) * runif(n)
      if(reference.type == "unimodal")
        X.unimodal[,j] <- unimodal.gkde.sample(unimodal.fits[[j]], n)
    }
    ## browser()
    mdsfun <- gsl.make.nn.mdsfun(X.unimodal)
    gsl.cluster.out <- gsl.cluster(X.unimodal, mdsfun, assign.fluff = F)
    rs.unimodal[i,] <- gsl.runt.size(gsl.cluster.out)[1:nrs]
    cat(" ", i)
  }
  return(rs.unimodal)
}

##-----------------------------------------------------------------

assess.runt.sizes <- function(X, scenario, reference = reference,
                              nrep = nrep, nrs = nrs) {
  XX <- prcomp(X)$x
  rs.reference <- runt.size.reference(XX, nrep, reference = reference)
  mdsfun <- gsl.make.nn.mdsfun(X)
  gsl.cluster.out <- gsl.cluster(X, mdsfun, assign.fluff = F)
  rs.original <- gsl.runt.size(gsl.cluster.out)[1:nrs]
  rs.significance <- rep(0, nrs)
  for (i in 1:nrs) {
    rs.significance[i] <- sum(rs.reference[,i] >= rs.original[i]) / nrep
  }
##   if (make.ps) {
##     figure.filename <- paste(scenario, ".ps", sep = "")
##     figure.pathname <- paste(figure.dir, figure.filename, sep = dir.sep)
##     postscript(file = figure.pathname, onefile = T, horizontal = F,
##     width = 6, height = 6)
##   }
  plot(rs.original, pch = 20, col = "red")
  bplot(rs.reference, add = T)
  title(paste(scenario, ", rs"))

  inc.original <- rs.original[1:(nrs-1)] - rs.original[2:nrs]
  inc.reference <- rs.reference[, 1:(nrs-1)] - rs.reference[, 2:nrs]
  inc.significance <- rep(0, (nrs-1))
  for (i in 1:(nrs-1)) {
    inc.significance[i] <- sum(inc.reference[,i] >= inc.original[i]) / nrep
  }
  plot(inc.original, pch = 20, col = "red")
  bplot(inc.reference, add = T)
  title(paste(scenario, ", rs inc"))
 
  relinc.original <- (rs.original[1:(nrs-1)] - rs.original[2:nrs]) /
                      rs.original[2:nrs]
  relinc.reference <- (rs.reference[, 1:(nrs-1)] - rs.reference[, 2:nrs]) /
                      rs.reference[, 2:nrs]
  relinc.significance <- rep(0, (nrs-1))
  for (i in 1:(nrs-1)) {
    relinc.significance[i] <- sum(relinc.reference[,i] >= relinc.original[i]) / nrep
  }
  plot(relinc.original, pch = 20, col = "red")
  bplot(relinc.reference, add = T)
  title(paste(scenario, ", rs rel inc"))
  return(list(rs.significance = rs.significance,
              inc.significance = inc.significance,
              relinc.significance = relinc.significance))
}



#################################################################
## Stuff below does not work. Density won't be unimodal
## Make unimodal data with approximately the same covariance as X but
## unimodalized marginals

## plot(X, pch = 20)
## var(X)
## n <- nrow(X)
## m <- ncol(X)
## hist(X[,1], col = "green")
## hist(X[,2], col = "green")

## U <- sphere(matrix(runif(n*m), nrow = n))
## R <- chol(var(X))
## Y <- U %*% R
## plot(Y, pch = 20)
## points(X, pch = 20, col = "green")
## var(Y)



## fu.out <- fit.unimodals(X)
## z1 <- unimodal.gkde.sample(fu.out[[1]], n)
## z2 <- unimodal.gkde.sample(fu.out[[2]], n)

## hist(z1, col = "green")
## hist(X[,1], col = "green")

## hist(z2, col = "green")
## hist(X[,2], col = "green")

## Z <- X
## Z[order(Y[,1]), 1] <- sort(z1)
## Z[order(Y[,2]), 2] <- sort(z2)

## var(Z)
## var(X)

## cor(X, method = "spearman")
## cor(Y, method = "spearman")
## cor(Z, method = "spearman")

## plot(Z, pch = 20)
## hist(Z[,1], col = "green")
## hist(Z[,2], col = "green")


## ## Transform X to be more highly correlated

## Sigma <- matrix(c(1, 0.8, 0.8, 1), nrow = 2)
## XX <- sphere(X) %*% t(chol(Sigma))
## plot(XX, pch = 20)

## X <- XX


#################################################################
## Stuff below is deprecated

## unimodal.regression.random <- function(x, y, mode.x) {
##   merge.blocks <- function(i) {
##     blocks[i, "upper.index"] <<- blocks[i+1, "upper.index"]
##     blocks[i, "m"] <<- (blocks[i, "count"] * blocks[i, "m"] +
##                         blocks[i+1, "count"] * blocks[i+1, "m"]) /
##                           (blocks[i, "count"] + blocks[i+1, "count"])
##     blocks[i, "count"] <<- blocks[i, "count"] + blocks[i+1, "count"]
##     if((i+2) <= n) 
##       for (j in (i+2):n) blocks[(j-1),] <<- blocks[j,]
##     nblocks <<- nblocks - 1
##   }
##   pick.block <- function() {
##     candidates <- NULL
##     for (i in 1:(nblocks - 1))
##       if (is.violation(blocks, i, mode.x)) candidates <- c(candidates, i)
##     if (length(candidates) <= 1) return(candidates)
##     else return(sample(candidates, 1))
##   }
##   y <- y[order(x)]
##   n <- length(x)
##   nblocks <- n
##   blocks <- matrix(0, nrow = n, ncol = 4)
##   colnames(blocks) <- c("lower.index", "upper.index", "count", "m")
##   blocks[,"lower.index"] <- 1:n
##   blocks[,"upper.index"] <- 1:n
##   blocks[,"count"] <- rep(1, n)
##   blocks[,"m"] <- y

##   picked.block <- pick.block()
##   while(!is.null(picked.block)) {
##     merge.blocks(picked.block)
##     picked.block <- pick.block()
##   }
##   yhat <- rep(0, n)
##   for (i in 1:nblocks)
##     yhat[blocks[i, "lower.index"]:blocks[i, "upper.index"]] <- blocks[i, "m"]
##   return(list(blocks = blocks[1:nblocks,], yhat = yhat))
## }

## ##-----------------------------------------------------------------

## unimodal.regression <- function(x, y, mode.x) {
##   merge.blocks <- function(i) {
##     blocks[i, "upper.index"] <<- blocks[i+1, "upper.index"]
##     blocks[i, "m"] <<- (blocks[i, "count"] * blocks[i, "m"] +
##                         blocks[i+1, "count"] * blocks[i+1, "m"]) /
##                           (blocks[i, "count"] + blocks[i+1, "count"])
##     blocks[i, "count"] <<- blocks[i, "count"] + blocks[i+1, "count"]
##     if((i+2) <= nblocks) 
##       blocks[(i+1):(nblocks-1),] <<- blocks[(i+2):nblocks,]
##     nblocks <<- nblocks - 1
##   }

##   ox <- order(x)
##   y <- y[ox]
##   n <- length(x)
##   nblocks <- n
##   blocks <- matrix(0, nrow = n, ncol = 4)
##   colnames(blocks) <- c("lower.index", "upper.index", "count", "m")
##   blocks[,"lower.index"] <- 1:n
##   blocks[,"upper.index"] <- 1:n
##   blocks[,"count"] <- rep(1, n)
##   blocks[,"m"] <- y

##   i <- 1
##   while (i < nblocks) {
##     if (!is.violation(blocks, i, mode.x)) {i <- i+1; next}
##     merge.blocks(i)
##     if (i > 1) {
##       for (j in (i-1):1) {
##         if (!is.violation(blocks, j, mode.x)) break
##         merge.blocks(j)
##         i <- i-1
##       }
##     }
##   }
##   yhat <- rep(0, n)
##   for (i in 1:nblocks)
##     yhat[blocks[i, "lower.index"]:blocks[i, "upper.index"]] <- blocks[i, "m"]
##   yhat[ox] <- yhat
##   return(list(blocks = blocks[1:nblocks,], yhat = yhat))
## }


##-----------------------------------------------------------------
## Let u[i], v[i] be the lower and upper indices for block i.
## We assume that the regression function is constant in
## [x[u[i]], x[u[i+1]]) for i = 1...(nblocks - 1)
## and [x[u[nblocks-1]], x[u[nblocks]]]

## is.violation <- function(blocks, i, mode.x) {
##   left.of.mode <- function(i) {
##     return(x[blocks[i+1, "lower.index"]] <= mode.x)
##   }
##   right.of.mode <- function(i) {
##     return(x[blocks[i, "lower.index"]] > mode.x)
##   }
##   violation <- FALSE
##   if (left.of.mode(i) & (blocks[i, "m"] > blocks[i+1, "m"]))
##     violation <- TRUE
##   if (right.of.mode(i+1) &  (blocks[i, "m"] < blocks[i+1, "m"]))
##     violation <- TRUE
##   return(violation)
## }

##-----------------------------------------------------------------











