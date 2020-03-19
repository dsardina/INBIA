
collective.influence <- function(g, ball.size = 2){
  # First line description
  # 
  # Args:
  #   var1: Description
  #   var2: Description
  #   
  # Returns:
  #   Output description
  diam <- diameter(graph = g, directed = "FALSE", is.connected(g))
  
  if(ball.size > diam) stop("Current ball size is greater than the input graph's diameter.")
  
  ci           <- vector(mode = "list", length = length(V(g)))
  names(ci)    <- V(g)$name
  balls        <- vector(mode = "list", length = length(V(g)))
  names(balls) <- V(g)$name
                  
  for(i in 1:length(V(g))){
    if(i %% 10 == 0) cat("Nodes processed:", i, "\n")
    ball          <- c()
    current.nodes <- i
    
    for (j in 1:ball.size){
      
      current.nodes <- unique( unlist( get.neighbors(g, current.nodes) ) )
      ball          <- c(ball, current.nodes)
      
    }
    
    ball        <- setdiff(ball, i)
    ball.degree <- get.ball.degree(g, ball)
    vdegree     <- degree(g, i) -1
    
    balls[[i]] <- ball
    ci[[i]]    <- unname(vdegree) * sum( ball.degree -1 )
  }
  
  return( list(ci=unlist(ci), balls=balls) )
}

get.neighbors <- function(g, vs){
  # First line description
  # 
  # Args:
  #   var1: Description
  #   var2: Description
  #   
  # Returns:
  #   Output description
  lapply(vs, function(v){
    return(neighbors(graph = g, v = v, mode = "all"))
  })
}

get.ball.degree <- function(g, ball){
  # Computes the degree of the nodes inside the ball in the graph.
  # 
  # Args:
  #   g: The graph used for collective influence
  #   ball: the ids of the nodes present in the ball
  #   
  # Returns:
  #   A numeric vector of degrees for the node in the ball
  degrees <- lapply(ball, function(n){
    degree(g, V(g)[n])
  })
  
  return(unlist(degrees))
}


ci.pgraph <- function(g, ball, v){
  # Plots the graph and color the nodes of the ball of the vertex v
  #
  # Args:
  #   g: The graph used for collective influence
  #   ball: A numeric vector with vertices ids within the ball
  #   v: An igraph vertex within the graph used for collective influence
  #
  # Returns:
  #   The plot of the graph with normalized coordinates and colored nodes
  #   of the ball of vertex v
  
  l <- layout.lgl(g)
  l <- norm_coords(l, ymin=-1, ymax=1, xmin=-1, xmax=1)

  ball.color <- rep(FALSE, length(V(g)))
  ball.color[ball] <- TRUE

  V(g)$color <- ifelse(ball.color, "red", "orange")
  V(g)$color[v] <- "palegreen3"

  plot(g, rescale = FALSE, layout = l*1.2, vertex.size = 5)
}

# nv <- 40
# g <- erdos.renyi.game(nv, 0.2, "gnp")
# g <- sample_pa(nv, directed = FALSE)
# 
# ci <- collective.influence(g)
# ci.pgraph(g, ci$balls[[4]], V(g)[4])
# 
# source("scripts/method.functions.R")
# plot.graph(g, ci$ci)
