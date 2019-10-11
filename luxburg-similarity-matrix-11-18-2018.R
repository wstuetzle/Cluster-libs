## Compute Luxburg similarity graph
## ================================
##
## Vertex weights w_i = knn density estimates
##
## Two versions for edges:
## * Symmetric knn graph: connect vertices i and j if i is among k
##   nearest neighbors of j AND vice versa
## * Asymetric knn graph: connect vertices i and j if i is among k
##   nearest  neighbors of j OR vice versa
##
## Edge weights w_ij = min(w_i, w_j)
##
## GSL computes mast of similarity graph but assumes the graph is
## connected. However, a knn graph does not have to be connected.
##
## Idea: Compute Euclidean mst. Suppose (i,j) is a Euclidan mst
## edge. If (i,j) is not an edge of the knn graph then add (i,j) with
## weight 0 to the edges of the knn graph. The resulting augmented knn
## graph will be connected. The connected components of the threshold
## graph for threshold 0 are the connected components of the knn graph.
## 
## An edge with weight 0 will only be added to the mast of the
## augmented knn graph if there is no edge with weight > 0 connecting
## any out vertex to any in vertex.  => A Euclidean mst edge will only
## be added to the mast of the knn graph if the current in vertices and
## the current out vertices are in different connected components of
## the knn graph.
##
## Need a function luxsfun(i, js) that computes the similarity between
## a vertex i and a vector of vertices.


