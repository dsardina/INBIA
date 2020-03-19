
mds <- function(g, dg = NULL){
  # Algorithm for finding Minimum Dominating Set.
  # Uses a greedy approach with approximation error less than ln(delta) -2
  # 
  # Args:
  #   g: The graph used to compute approximate Minimum Dominating Set
  #   
  # Returns:
  #   A list with a vector of node ids belonging to the MDS and
  #   the visited node ids
  nv       <- length(V(g))
  mds      <- rep(FALSE, nv)
  visited  <- rep(FALSE, nv)
  consider <- rep(TRUE, nv)
  
  if(is.null(dg)) dg <- degree(g)
  
  while ( !all(visited) ) {
    dcons <- dg[consider]
    mids  <- which(dg == max(dcons))
    mid   <- mids[consider[mids]][1]
    
    mds[mid]      <- TRUE
    visited[mid]  <- TRUE
    nids          <- neighbors(g, mid, "all")
    visited[nids] <- TRUE
    
    consider[visited] <- FALSE
  }
  
  return(list(mds=mds, visited=visited))
}

mds.pgraph <- function(g, visited, mds){
  # First line description
  # 
  # Args:
  #   var1: Description
  #   var2: Description
  #   
  # Returns:
  #   Output description
  V(g)$color <- ifelse(visited, "orange", "blue")
  clr         <- V(g)$color
  clr[mds]    <- "red"
  V(g)$color <- clr
  
  l <- layout_with_fr(g)
  l <- norm_coords(l, ymin=-1, ymax=1, xmin=-1, xmax=1)
  plot(g, rescale = FALSE, layout = l*1.2, vertex.size = 5, vertex.label.cex = 0.8)
}

# nv  <- 153
# g   <- erdos.renyi.game(nv, 0.2, "gnp")
# g   <- sample_pa(nv, directed = FALSE)
# ci <- collective.influence(g)
# res <- mds(g, dg = ci$ci)
# mds.pgraph(g, res$visited, res$mds)
 
 
# netm <- get.adjacency(g, sparse=F)
# palf <- colorRampPalette(c("gold", "dark orange")) 
# heatmap(netm, Rowv = NA, Colv = NA, col = palf(100))
# 
# 
# ## Solve the binary integer linear programming problem
# library(lpSolve)
# 
# vertex.num <- 100
# obj.coef <- rep(1, vertex.num)
# A <- rep(1, constraints)
# 
