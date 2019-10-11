## Generalized single linkage clustering: an illustration
## ======================================================
## Werner Stuetzle (wxs@u.washington.edu), 3-19-2009

## This document assumes that you have read the paper

## W. Stuetzle and R. Nugent. A generalized single linkage method for
## estimating the cluster tree of a density. Journal of Computational
## and Graphical Statistics, 2009 (to appear).

##-----------------------------------------------------------------

library(MASS)

quartz()
echo <- T

cluster.dir <- "/Users/wxs/Dropbox/Clustering"
libs.dir <- paste(cluster.dir, "A-Git-repos", "Cluster-libs", sep = "/")
data.dir <- paste(cluster.dir, "Data", "Test-collection", sep = "/")

source(paste(libs.dir, "gsl-functions-6-8-2018.R", sep = "/"), echo = echo)
source(paste(libs.dir, "gsl-functions-for-figures-1-19-09.R", sep = "/"), echo = echo)

##-----------------------------------------------------------------
## Read the olive oil data. These are measurements of the
## concentrations of 8 chemicals in 572 olive oil samples from 9
## different areas of Italy. The areas are grouped into
## regions. Region 1 consists of areas 1..4, region 2 consists of
## areas 5 and 6, and region 3 consists of areas 7..9.
##
## We will try to recover the areas and regions frm the measurements

filename <- "olive.R"
source(paste(data.dir, filename, sep = "/"), echo = echo)

## Creates a predictor array "X" and a label vector "group.id" with
## values 1, 2..ngroup

print(n <- nrow(X))
print(m <- ncol(X))
print(ngroup <- max(group.id))

table(group.id)

##-----------------------------------------------------------------
## First we must choose a density estimate. The package provides two
## options: the nearest neighbor estimate, and a kernel estimate with
## spherical Gaussian kernel.

##-----------------------------------------------------------------
## Let's first try the nearest neighbor estimate. We use the function
## gsl.make.nn.mdsfun to generate a "minimum density similarity
## function" nn.mdsfun(i, j), which returns the minimum of the nearest
## neighbor density estimate over the edge (line segment) connecting
## observations i and j. (Actually, that's not really true, but don't
## worry about it - the graph cluster tree comes out right.

nn.mdsfun <- gsl.make.nn.mdsfun(X)

## Next we use gsl.cluster to compute the graph cluster tree, the runt
## sizes and excess masses, and assorted other stuff. The output of
## gsl.cluster is not intended for your consumption. 

gsl.cluster.out <- gsl.cluster(X, nn.mdsfun)

## Let's have a look at the runt sizes

gsl.runt.size(gsl.cluster.out)[1:10]

## Here's what you should see
##  [1] 168  97  59  51  42  42  33  13  13  12

## We scan the runt sizes form small to large and look for the first
## gap (see Section 5 and Section 7 of the paper). It seems to be
## between 13 and 33. So let's choose 33 as the runt size
## threshold. This means the pruned graph cluster tree will have 7
## internal nodes and 8 leaves. Let's compute the pruned tree.

gsl.cluster.out <- gsl.cluster(X, nn.mdsfun, pruning.criterion =
                               "size", pruning.threshold = 33)

## Draw the pruned graph cluster tree. The numbers below the interior
## nodes are the runt sizes. The numbers above the leaves are the leaf
## codes 

bla <- gsl.draw.cluster.tree(gsl.cluster.out, draw.runt.size = T,
                             draw.cluster.leaf.codes = T)

## Let's do some diagnostics to assess cluster separation. Before
## proceeding, have a look at the help files for the functions
## gsl.gkl.diagnostic.plot and gsl.subtree.gkl.diagnostic.plot

gsl.gkl.diagnostic.plot(X, gsl.cluster.out, draw.runt.size = T,
                        draw.cluster.leaf.codes = T)

## Click on the root node. You will see two clearly bimodal
## histograms. So the daughters of the root node are are clearly
## separated in feature space.

## Do the same thing again, but click on the left daughter of the root
## node (the node with runt size 42)

gsl.gkl.diagnostic.plot(X, gsl.cluster.out, draw.runt.size = T,
                        draw.cluster.leaf.codes = T)

## That also looks clearly bimodal. We could proceed in this way, but
## instead look at the three clusters which are left descendents of
## the root node.

gsl.subtree.gkl.diagnostic.plot(X, gsl.cluster.out)

## Click on the left daughter of the root node. The projection of the
## observations in the left daughter on the Fisher plane clearly shows
## three groups. Now hit "Enter" and you see that the three groups
## really do correspond to the three clusters identified by the
## procedure

## Next, let's check if the three clusters really correspond to
## different areas

obs.in.left.daughter.of.root <-
  gsl.visually.select.observations.in.node(gsl.cluster.out)$in.node

## Click on left daughter of root

cluster.labels <- gsl.observation.labels(gsl.cluster.out)

table(group.id[obs.in.left.daughter.of.root],
      cluster.labels[obs.in.left.daughter.of.root])

## Most of the left children are from areas 7, 8, 9, with some
## smattering of areas 1, 2, and 4

## Now let's look at the right daughter of the root node

gsl.gkl.diagnostic.plot(X, gsl.cluster.out, draw.runt.size = T,
                        draw.cluster.leaf.codes = T)

## Clearly bimodal. We can keep playing this game and see which splits
## really separate distinct groups in the data and which ones appear
## spurious.

## Let's look at the correspondence between areas and cluster labels

diagonalize.table(table(group.id, cluster.labels), shuffle = "columns")

## Not perfect: The areas in region 1 seems to be hard to separate.
## Let's see if we can do better if we restrict ourselves to
## observations in the cluster cores

in.core <- gsl.observation.in.core(gsl.cluster.out)
diagonalize.table(table(cluster.labels[in.core], group.id[in.core]))

## Clusters 5, 8, 9, 12, 14, 15 mostly correspong to areas 8, 7, 9, 2,
## 6, 5, respcectively.

##=================================================================
## Now let's try a kernel density estimate

## Sphere the data - a good idea if we want to use a kernel density
## estimate and determine the bandwidth by least squares
## cross-validation

X <- sphere(X)
var(X)

##=================================================================
## Find optimal bandwidth for kernel density estimate using least
## squares CV
hopt <- cv.search(X)$opt.smopar
hopt
## 0.2287662

## Create minimum density similarity function. This takes a few seconds.
kernel.mdsfun <- gsl.make.kernel.mdsfun(X, hopt)

gsl.cluster.out <- gsl.cluster(X, kernel.mdsfun)

gsl.runt.size(gsl.cluster.out)[1:10]

## Here's what you should see
## 129  89  47  33  25  24  24  20  11   9

## We scan the runt sizes form small to large and look for the first
## gap (see Section 5 and Section 7 of the paper). It seems to be
## between 11 and 20. So let's choose 20 as the runt size
## threshold. This means the pruned graph cluster tree will have 8
## internal nodes and 9 leaves. Let's compute the pruned tree.

gsl.cluster.out <- gsl.cluster(X, kernel.mdsfun, pruning.criterion =
                               "size", pruning.threshold = 20)

## Now we just re-run the code above.

#################################################################
#################################################################
## Code below edited 8-9-2019
## The code below seems to be for the Iris data

filename <- "iris.R"
source(paste(data.dir, filename, sep = "/"), echo = echo)

## Creates a predictor array "X" and a label vector "group.id" with
## values 1, 2..ngroup

print(n <- nrow(X))
print(m <- ncol(X))
print(ngroup <- max(group.id))

table(group.id)


##-----------------------------------------------------------------
## First we must choose a density estimate. The package provides two
## options: the nearest neighbor estimate, and a kernel estimate with
## spherical Gaussian kernel.

##-----------------------------------------------------------------
## Let's first try the nearest neighbor estimate. We use the function
## gsl.make.nn.mdsfun to generate a "minimum density similarity
## function" nn.mdsfun(i, j), which returns the minimum of the nearest
## neighbor density estimate over the edge (line segment) connecting
## observations i and j. (Actually, that's not really true, but don't
## worry about it - the graph cluster tree comes out right.

nn.mdsfun <- gsl.make.nn.mdsfun(X)

## Next we use gsl.cluster to compute the graph cluster tree, the runt
## sizes and excess masses, and assorted other stuff. The output of
## gsl.cluster is not intended for your consumption. 

gsl.cluster.out <- gsl.cluster(X, nn.mdsfun)

## Let's have a look at the runt sizes

gsl.runt.size(gsl.cluster.out)[1:10]

## Here's what you should see
##  [1] 50 38 16 14 14  7  5  4  4  4

## We scan the runt sizes form small to large and look for the first
## gap (see Section 5 and Section 7 of the paper). It seems to be
## between 38 and 16. So let's choose 38 as the runt size
## threshold. This means the pruned graph cluster tree will have 2
## internal nodes and 3 leaves. Let's compute the pruned tree.

gsl.cluster.out <- gsl.cluster(X, nn.mdsfun, pruning.criterion =
                               "size", pruning.threshold = 38)

gsl.cluster.out <- gsl.cluster(X, nn.mdsfun, pruning.criterion =
                               "size", pruning.threshold = 14)

## Draw the pruned graph cluster tree. The numbers below the interior
## nodes are the runt sizes. The numbers above the leaves are the leaf
## codes 

bla <- gsl.draw.cluster.tree(gsl.cluster.out, draw.runt.size = T,
                             draw.cluster.leaf.codes = T)

## Let's do some diagnostics to assess cluster separation. Before
## proceeding, have a look at the help files for the functions
## gsl.gkl.diagnostic.plot and gsl.subtree.gkl.diagnostic.plot

gsl.gkl.diagnostic.plot(X, gsl.cluster.out, draw.runt.size = T,
                        draw.cluster.leaf.codes = T)

## Click on the root node. You will see two clearly bimodal
## histograms. So the daughters of the root node are are clearly
## separated in feature space.

## Do the same thing again, but click on the right daughter of the root
## node (the node with runt size 38)

gsl.gkl.diagnostic.plot(X, gsl.cluster.out, draw.runt.size = T,
                        draw.cluster.leaf.codes = T)
## Clearly bimodal

## Now the left daughter of the root node

gsl.gkl.diagnostic.plot(X, gsl.cluster.out, draw.runt.size = T,
                        draw.cluster.leaf.codes = T)

## Also looks bimodal. Let's see how many observations there are

obs.in.left.daughter.of.root <-
  gsl.visually.select.observations.in.node(gsl.cluster.out)$in.node

## Click on the left daughter of the root node

sum(obs.in.left.daughter.of.root)
## [1] 50

## Only 50 observations. So the default of 20 breaks in the histograms
## may be too many. Let's try 10

gsl.gkl.diagnostic.plot(X, gsl.cluster.out, draw.runt.size = T,
                        draw.cluster.leaf.codes = T, n.breaks = 10)

## Not much bimodality left. Now check out the right daughter of the
## root node

gsl.gkl.diagnostic.plot(X, gsl.cluster.out, draw.runt.size = T,
                        draw.cluster.leaf.codes = T)

## Click on the right daughter of the root node. Pretty clear
## bimodality. I know, I know - it helps to know what one is looking
## for. Now check out the right grand-daughters of the root node - you
## get the idea

## We can also assess the separation between the leaves of a subtree

gsl.subtree.gkl.diagnostic.plot(X, gsl.cluster.out)

## Click on the right daughter of the root node. Focus on the black
## dots in the scatterplot, which are the observations in the leaf
## cores. There seem to be two clusters, not four. Then hit "return",
## and the observations in the leaves will be shown in different colors.

## Maybe a better runt size threshold would be 38, which would give a
## pruned cluster tree with three leaves

gsl.cluster.out <- gsl.cluster(X, nn.mdsfun, pruning.criterion =
                           "size", pruning.threshold = 38)

## Let's get the cluster labels and see how they align with the
## species

observation.labels <- gsl.observation.labels(gsl.cluster.out)
table(observation.labels, group.id)
##                   Species
## observation.labels setosa versicolor virginica
##                  2     50          0         0
##                  6      0         45         1
##                  7      0          5        49

## Well how about it? The labels reflect the species almost perfectly.
## Let's check if the agreement gets better if we only look at
## observations in the leaf cores?

in.core <- gsl.observation.in.core(gsl.cluster.out)
table(observation.labels[in.core], group.id[in.core]
##     Setosa versicolor virginica
##   2     50          0         0
##   6      0         39         0
##   7      0          3        35

## There are still some Virginica that end up in the Versicolor
## cluster.

## Let's look how well the clusters obtained with runt size threshold
## 14 match up with the species

gsl.cluster.out <- gsl.cluster(X, nn.mdsfun, pruning.criterion =
                           "size", pruning.threshold = 14)
observation.labels <- gsl.observation.labels(gsl.cluster.out)
table(observation.labels, group.id)
##                   Species
## observation.labels setosa versicolor virginica
##                 4      30          0         0
##                 5      20          0         0
##                 12      0         30         1
##                 13      0         15         0
##                 14      0          0        36
##                 15      0          5        13

## Each of the species is broken up into two clusters.

##-----------------------------------------------------------------

## Now try a kernel density estimate. We use a standard Gaussian
## kernel on the sphered data, which is the same as using a Gaussian
## kernel with the same covariance as the data.

XX <- sphere(X)

## Use least squares cross-validation to determine the bandwidth

cv.search.out <- cv.search(XX)
h <- cv.search.out$opt.smopar
h
##[1] 0.3510157

## Now make a minimum density similarity function based on the kernel
## density estimate

kernel.mdsfun <- gsl.make.kernel.mdsfun(XX, h, interactive = T)

## Setting "interactive = T" makes the function print out a "progress
## report" for the impatient.

gsl.cluster.out <- gsl.cluster(XX, kernel.mdsfun)
gsl.runt.size(gsl.cluster.out)[1:10]
## [1] 46  7  6  5  5  3  3  3  2  2

round(nrow(XX) * gsl.runt.excess.mass(gsl.cluster.out), 0)[1:10]
##  [1] 27  2  2  2  2  2  1  1  1  1

## Runt excess mass is easier to scan if we convert it from the
## probability scale to the "equivalent number of observations" scale.

## Both runt size and runt excess mass suggest two clusters.

gsl.cluster.out <- gsl.cluster(XX, nn.mdsfun, pruning.criterion =
                               "size", pruning.threshold = 46)

## Let's get the cluster labels and see how they align with the
## species

observation.labels <- gsl.observation.labels(gsl.cluster.out)
table(observation.labels, group.id)
##                    Species
## observation.labels setosa versicolor virginica
##                  2     50          0         0
##                  3      0         50        50


## Seems like we are oversmoothing? We have observed more generally
## that nearest neighbor estimation tends to do not much worse, and
## often does better, than kernel density estimation.
