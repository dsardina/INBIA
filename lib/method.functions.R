
# FUNCTION SECTION ================================================================

# saveTOMPlot - Save the TOM plot from WGCNA
# 
#
saveTOMPlot <- function(TOM, geneTree, moduleColors, tname, fname){
  png(fname, width = 6, height = 5, unit = "in", res = 300)
  on.exit(dev.off())
  TOMplot(TOM, geneTree, moduleColors, main = paste0("Network heatmap for ", tname))
}

#
#
#
saveDendogram <- function(dendo, mcolors, fname){
  png(fname, width = 6, height = 5, unit = "in", res = 300)
  on.exit(dev.off())
  
  # Plot the dendrogram and the module colors underneath
  plotDendroAndColors(dendo, mcolors,
                      "Module colors",
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05)
}

#
#
#
saveScalefreePlot <- function(sft, fname, powers){
  png(fname, width = 9, height = 5, unit = "in", res = 300)
  on.exit(dev.off())
  par(mfrow = c(1,2))
  cex1 <- 0.9;
  
  # Scale-free topology fit index as a function of the soft-thresholding power
  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
       main = paste("Scale independence"))
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,cex=cex1,col="red")
  
  # this line corresponds to using an R^2 cut-off of h
  abline(h=0.90, col="red")
  
  # Mean connectivity as a function of the soft-thresholding power
  plot(sft$fitIndices[,1], sft$fitIndices[,5],
       xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
       main = paste("Mean connectivity"))
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
}

#
#
#
wgcna <- function(dataset, tname, dynMods = FALSE, blockMods = FALSE){
  powers <- c(c(1:10), seq(from = 12, to=20, by=2))
  
  # Compute soft threshold
  cat("Pick Soft Threshold Table:\n")
  sft <- pickSoftThreshold(dataset, powerVector = powers)
  
  # Save the plot for further visual analysis
  dir.create("results/wgcna", recursive = TRUE, showWarnings = FALSE)
  saveScalefreePlot(sft, paste0("results/wgcna/", tname, ".scalefree.png"), powers)
  
  # Retrieve the best value for power
  best.power <- which.max(-sign(sft$fitIndices[,3])*sft$fitIndices[,2])
  cat("Best power:", best.power, "\n")
  
  
  if(blockMods){
    # Construct the network and identifying modules
    # NOTE: read the help before execute this function on new dataset
    # see pamRespectsDendro
    net <- blockwiseModules(dataset, power = best.power,
                            TOMType = "unsigned", minModuleSize = 30,
                            reassignThreshold = 0, mergeCutHeight = 0.25,
                            numericLabels = TRUE, pamRespectsDendro = FALSE)
    
    # Convert modules to colors
    moduleColors <- labels2colors(net$colors)
    
    # Create and save the dendograms
    saveDendogram(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]], paste0("results/wgcna/", tname, ".moduleDendogram.png")) 
  }
  
  TOM <- TOMsimilarityFromExpr(dataset, networkType = "unsigned", TOMType = "unsigned", power = best.power)
  
  if (dynMods){
    
    # Compute dissimilarity matrix
    dissTOM   <- 1-TOM
    geneTree  <- hclust(as.dist(dissTOM), method = "average")
    
    # Use dynamicTreeCutting
    minModuleSize <- 10
    
    # Module identification using dynamic tree cut:
    dynamicMods <- cutreeDynamic(dendro = geneTree,
                                 distM = dissTOM,
                                 method = "tree",
                                 minClusterSize = minModuleSize)
    
    dynamicColors <- labels2colors(dynamicMods)
    
    saveDendogram(geneTree, dynamicColors, paste0("results/wgcna/", tname, ".dynamicModuleDendogram.png"))
    
    # Transform dissTOM with a power to make moderately strong connections more visible in the heatmap
    plotTOM <- dissTOM^7
    
    # Set diagonal to NA for a nicer plot
    diag(plotTOM) <- NA
    
    # Call the plot function
    cat("Saving TOM plots...\n")
    # saveTOMPlot((1-adjacency)^7, adjGeneTree, moduleColors, tname, paste0("wgcna/", tname, ".Adjplot.png"))
    saveTOMPlot(plotTOM, geneTree, dynamicColors, tname, paste0("results/wgcna/", tname, ".TOMplot.png"))
  }
  
  # create matrix and edges from predictions
  return(TOM)
}

# invcov2parcor - Invert the partial correlation matrix
#
# From Multi-Method approach source. Convert inverse covariance matrix to partial correlation matrix 
invcov2parcor <- function(mat){
  vars <- diag(mat)
  len <- length(vars)
  tmp <- pcor <- matrix(NA, nrow=len, ncol=len)
  for(m in 1:len){
    tmp[m, ] <- mat[m, ]/sqrt(vars[m])   
  }#end for m
  for(m in 1:len){
    pcor[, m] <- tmp[, m]/sqrt(vars[m])   
  }#end for m
  return(pcor)
}

#
#
#
elasticNetGridSearch = function(x, y, nfold=10, alpha=0.5){
  cv.alpha <- lambda.opt <- matrix(NA, ncol = length(alpha), nrow = length(nfold))
  
  # grid search involving cross-validation splits and alpha values
  x <- as.matrix(x)
  y <- as.matrix(y)
  
  for(i in 1:length(nfold)){
    for(j in 1:length(alpha)){
      res <- cv.glmnet(x, y, alpha = alpha[j])
      # cross-validation error when lambda is the largest such that it's within 1se of predictive power
      cv.alpha[i,j] <- res$cvm[which(res$lambda == res$lambda.1se)]
      lambda.opt[i,j] <- res$lambda.1se
    }#end for j
  }#end for i
  opt.ind <- which(cv.alpha == min(cv.alpha), arr.ind=TRUE)
  k.opt <- nfold[opt.ind[1]]
  a.opt <- alpha[opt.ind[2]]
  lambda <- lambda.opt[opt.ind[1], opt.ind[2]]
  model <- glmnet(x, y, alpha = a.opt, intercept = FALSE)
  pred <- predict(model, newx = x, type = "coefficients", s = lambda)
  intercept <- pred[1]
  coefficients <- pred[-1]
  names(coefficients) <- rownames(pred)[-1]
  return(list(lambda.opt = lambda, cv.grid = cv.alpha, nfold.opt = k.opt,
              alpha.opt = a.opt, intercept = intercept, coefficients = coefficients))
}

#
#
#
elasticNetwork <- function(X, nfold=10, alpha=0.5, verbose=FALSE){
  
  X <- scale(X)
  B <- matrix(0, nrow = ncol(X), ncol = ncol(X))
  colnames(B) <- 1:ncol(X)
  pcor <-  NULL
  n.opt <- alpha.opt <- rep(NA, ncol(X))
  
  if (verbose == TRUE) {
    cat(paste("Performing local elastic net regressions\n"))
    cat(paste("Vertex no "))
  }
  for(i in 1:ncol(X)){
    if (verbose == TRUE) {
      if ((i/10) == floor(i/10)) {
        cat(paste(i, "..."))
      }
    }
    M <- X[,-i]
    p <- X[,i]
    fit <- elasticNetGridSearch(M, p, nfold = nfold, alpha = alpha)
    coefi <- fit$coefficients
    B[i,-i] <- coefi
    
    n.opt[i] <- fit$nfold.opt
    alpha.opt[i] <- fit$alpha.opt
  }#end for i
  pcor <- Beta2parcor(B, verbose = verbose)
  cat(paste("\n"))
  
  return(list(pcor = pcor, n.opt = n.opt, alpha.opt = alpha.opt))
  
}

#
#
#
get.edges.from.path <- function(path){
  # Extract all the proteins in the path
  pp <- get.proteins(path)
  
  # Create the data frame of edges
  edges <- lapply(1:(length(pp)-1), function(i){
    return(c(pp[i], pp[i+1]))
  })
  
  return(do.call(rbind, edges))
}

#
#
#
get.prediction.graph <- function(paths){
  res <- lapply(paths, get.edges.from.path)
  
  return(graph.edge(do.call(rbind, res)))
}

#
#
#
load.tcpa.map <- function(){
  # Load function for tcpa map table
  # convTab <- read.delim("conversion/conversion.table.txt", header = FALSE)
  convTab <- read.delim("conversion/conversion.table_new.txt", header = FALSE)
  rownames(convTab) <- convTab[,1]
  
  return(convTab)
}

# Find shared gene symbols between interactions and TCPA
#
#
commonGenes <- function(intr, data){
  # take all gene namems from interactions set
  intr_genes <- unique(c(intr[,1], intr[,2]))
  
  # take all gene names from given TCPA dataset removing first 2 columns
  data_genes <- colnames(data)
  
  return(sum(data_genes %in% intr_genes))
}

# Given a vector of boolean values return the percentage of TRUE values
#
#
percentage <- function(values) return(sum(values)*100/length(values))

#
#
#
write.file <- function(data, fname, col.names = TRUE) write.table(x = data, file = fname, quote = FALSE, sep = "\t", row.names = FALSE, col.names = col.names)

#
#
#
subgraph <- function(interactions, gene.symbols){
  ids <- interactions[,1] %in% gene.symbols & interactions[,2] %in% gene.symbols
  return(interactions[ids,])
}

#
#
#
getEdges <- function(m){
  # Function devised for mutual information algorithms like aracne, clr, mrnet
  # m is symmetric adjacency matrix where the cells show the strength of the relationship  
  m[lower.tri(m, diag=TRUE)] <- 0
  
  res   <- cbind(which(m > 0, arr.ind=TRUE), value = m[m > 0])
  edges <- as.data.frame(res[order(abs(res[,3]), decreasing=TRUE),])
  rownames(edges) <- 1:nrow(edges)
  
  return(edges)
}

# Convert a vector of protein names
#
#
convert <- function(x, convt) return(convt[as.character(x),2])

#
#
#
is.phospho <- function(v) grepl("_p", v, fixed = TRUE)

#
#
#
gplot <- function(g){
  plot(g, layout = layout.fruchterman.reingold,
       vertex.size = 5,
       vertex.label.cex = .5, 
       edge.color="black") 
}

#
#
#
plot.graph <- function(g, size, vcolor=NULL){
  # l <- layout.lgl(g)
  l <- layout_with_fr(g)
  l <- norm_coords(l, ymin=-1, ymax=1, xmin=-1, xmax=1)
  lsize <- ((size - min(size)) / max(size - min(size)) * 10) + 2
  lsize[is.nan(lsize)] <- 0
  if(!is.null(vcolor) && length(vcolor) == length(V(g))) V(g)$color <- vcolor
  plot(g, rescale = FALSE, layout = l*1.2, vertex.size = lsize, vertex.label.cex = 0.8)
}

#
#
#
methods.sp.plot <- function(ds, stacked = TRUE){
  m <- ggplot(ds, aes(length, fill=method))
  if(stacked) m <- m + geom_bar()
  else m <- m + geom_bar(position="dodge")
  return(m)
}

#
#
#
graph.edge <- function(edges) return(graph.data.frame(as.data.frame(edges), directed=FALSE))

#
#
#
graph.matrix <- function(m) return(graph.adjacency(data.matrix(m), mode = "undirected", diag = F, weighted = TRUE))

#
#
#
get.spath <- function(g, v1, v2){
  if(all(c(v1, v2) %in% V(g)$name)) return(get.shortest.paths(g, v1, v2)$vpath)
  else return(NA)
}

#
#
#
get.path.length.distr <- function(paths){
  if(length(paths) == 0) return(numeric(0))
  distr <- lapply(paths, function(mp){
    if(is.na(mp)){
      return(NA)
    }
    else{
      lapply(unlist(strsplit(mp, split = "|", fixed = TRUE)), function(p){
        return(length(unlist(strsplit(p, split = ",", fixed = TRUE))))
      })
    }
  })
  return(unlist(distr))
}

#
#
#
save.plot <- function(distr, title, fname, type="jpg"){
  m <- ggplot(NULL, aes(x=distr))
  
  ggsave(fname, m + geom_histogram(binwidth=0.5) + ggtitle(title))
}

#
#
#
unique.genes <- function(edges) return(unique(c(edges[,1], edges[,2])))

#
#
#
compute.all.interactions <- function(geneset, with.repetition = FALSE){
  do.call(rbind, lapply(1:(length(geneset)-1), function(i){
    do.call(rbind, lapply((i+1):length(geneset), function(j){
      return(c(geneset[i], geneset[j]))
      #         cat(i,j,"\n")
    }))
  }))
}

#
#
#
which.edge <- function(edge, edge.list){
  res1 <- (as.character(edge.list[, 1]) == as.character(edge[1])) & (edge.list[, 3] == as.logical(edge[3])) & (as.character(edge.list[, 2]) == as.character(edge[2])) & (edge.list[, 4] == as.logical(edge[4]))
  res2 <- (as.character(edge.list[, 1]) == as.character(edge[2])) & (edge.list[, 3] == as.logical(edge[4])) & (as.character(edge.list[, 2]) == as.character(edge[1])) & (edge.list[, 4] == as.logical(edge[3]))
  
  if(sum(res1) > 1 || sum(res2) > 1) stop("Found more than one occurence!")
  
  if(any(res1)) return(which(res1))
  if(any(res2)) return(which(res2))
  
  stop("Error: match not found!")
}

#
#
#
get.not.predicted <- function(method.ds, genes){
  # cat("Computing not predicted\n")
  adj <- create.adj.mtx(method.ds, genes)
  
  not.predicted <- lapply(1:(nrow(adj)-1), function(i){
    res <- lapply((i+1):ncol(adj), function(j){
      if(!adj[i, j]) return(c(rownames(adj)[i], colnames(adj)[j]))
    })
    
    return(do.call(rbind, res))
  })
  return(do.call(rbind, not.predicted))
}

#
#
#
get.fn <- function(final.paths, method.ds, genes, cl){
  method.not.predicted <- get.not.predicted(method.ds, genes)
  if(is.null(method.not.predicted)) return(nrow(final.paths)) 
  
  clusterExport(cl, varlist = c("method.not.predicted", "final.paths"), envir = environment())
  
  res <- parLapply(cl, 1:nrow(method.not.predicted), function(index){
    tmp <- final.paths[final.paths[,1] == method.not.predicted[index, 1], ]
    return(method.not.predicted[index, 2] %in% tmp[,2])
  })
  
  return(sum(unlist(res)))
}

#
#
#
get.tn <- function(final.paths, method.ds, genes, cl){
#   cat("Creating not predicted ...\n")
  method.not.predicted <- get.not.predicted(method.ds, genes)
  if(is.null(method.not.predicted)) return(0)
  
  clusterExport(cl, varlist = c("method.not.predicted", "final.paths"), envir = environment())
  
  res <- parLapply(cl, 1:nrow(method.not.predicted), function(index){
    tmp <- final.paths[final.paths[,1] == method.not.predicted[index, 1], ]
    
    return(!(method.not.predicted[index, 2] %in% tmp[,2]))
  })
  
  return(sum(unlist(res)))
}

#
#
#
par.compute.shortest.paths <- function(edges, g){
  require(parallel)
  
  cat("Init parallel execution...\n", file = "par.res.txt", append = TRUE)
  cl <- makeCluster(detectCores() - 1)
  on.exit(stopCluster(cl))
  clusterExport(cl, list("edges", "g", "extract.genes","get.proteins",
                         "get.spath", "V","get.shortest.paths"), envir = environment())
  
  cat("Computation started...\n", file = "par.res.txt", append = TRUE)
  res <- parLapply(cl, 1:nrow(edges), function(i){
    if(i %% 1000 == 0) cat("Computed paths:", i, "\n", file = "par.res.txt", append = TRUE)
    # get all genes separated by comma
    gs1 <- extract.genes(as.character(edges[i, 1]))
    gs2 <- extract.genes(as.character(edges[i, 2]))
    
    # compute shortest paths between genes
    pres <- sapply(gs1, function(g1){
      return(sapply(gs2, function(g2){
        return(get.spath(g, g1, g2))
      }))
    })
    
    if(all(is.na(unlist(pres)))) return(NA)
    else{
      final <- lapply(pres, function(path){
        return(paste(path$name, collapse =","))
      })
      return(paste(unlist(final), collapse = "|"))
    }
  })
  
  return(res)
}

#
#
#
compute.shortest.paths <- function(edges, g){
  if(is.null(edges) || nrow(edges) == 0) return(numeric(0))
  
  res <- lapply(1:nrow(edges), function(i){
    # if(i %% 1000 == 0) cat("Computed paths:", i, "\n")
    # get all genes separated by comma
    gs1 <- extract.genes(as.character(edges[i,1]))
    gs2 <- extract.genes(as.character(edges[i,2]))

    # compute shortest paths between genes
    pres <- sapply(gs1, function(g1){
      return(sapply(gs2, function(g2){
        return(get.spath(g, g1, g2))
      }))
    })

    if(all(is.na(unlist(pres)))) return(NA)
    else{
      final <- lapply(pres, function(path){
        return(paste(path$name, collapse =","))
      })
      return(paste(unlist(final), collapse = "|"))
    }
  })
  return(unlist(res))
}

#
#
#
get.phospho.converted.interactions <- function(res, convTab){
  intr.ds <- data.frame(v1=convert(res[,1], convTab), v2=convert(res[,2], convTab),
                        v1.phospho = is.phospho(res[,1]), v2.phospho = is.phospho(res[,2]))
  return(intr.ds)
}

#
#
#
net.from.index <- function(tissue.net, gene.names){
  tissue.net[,1] <- gene.names[tissue.net[,1]]
  tissue.net[,2] <- gene.names[tissue.net[,2]]
  rownames(tissue.net) <- NULL
  return(tissue.net)
}

#
#
#
significant.edges <- function(tissue.net, gene.names){
  tissue.net$node1 <- gene.names[tissue.net$node1]
  tissue.net$node2 <- gene.names[tissue.net$node2]
  return(tissue.net[tissue.net$pval < 0.05,])
}

#
#
#
create.adj.mtx <- function(edges, ugenes = unique.genes(edges)){
  # cat("Creating adj matrix...\n")
  # ugenes <- unique.genes(edges)
  m <- matrix(FALSE, nrow = length(ugenes), ncol = length(ugenes), dimnames = list(ugenes, ugenes))
  if(nrow(edges) == 0) return(m)
    
  for(i in 1:nrow(edges)){
#     if(i %% 1000 == 0) cat(" ...",i)
    
    c1 <- as.character(edges[i,1])
    c2 <- as.character(edges[i,2])
    # cat(i, c1, c2, "\n")
    m[c1,c2] <- TRUE
    m[c2,c1] <- TRUE
  }
  return(m)
}

#
#
#
create.kedge.adj.mtx <- function(kedges){
  # cat("\nCreating matrix...")
  if(is.null(kedges) || is.na(kedges) || ncol(kedges) < 2){
    stop("Error: Edges must be at least a couple. Empty matrix.")
  }
  
  ugenes <- unique(as.vector(kedges))
  m      <- matrix(FALSE, nrow = length(ugenes), ncol = length(ugenes), dimnames = list(ugenes, ugenes))
  
  apply(kedges, 1, function(ke){
    for(i in 1:(length(ke)-1)){
      m[as.character(ke[i]), as.character(ke[i+1])] <<- TRUE
      # cat(as.character(ke[i]), "-", as.character(ke[i+1]), "\n")
    }
  })
  
  return(m)
}

# Return only those interactions whose proteins are in
# the set of genes given in input
#
#
filter.predictions <- function(preds, genes){
  if(length(preds) == 0) return(c())
  if(ncol(preds) < 2) stop("Wrong predictions format")
  
  good.ids <- apply(preds, 1, function(pred){
    f1 <- as.character(pred[1]) %in% genes & as.character(pred[2]) %in% genes # both names are in the set of genes
    f2 <- as.character(pred[1]) != as.character(pred[2]) # the interaction is not a loop
    return(f1 && f2)
  })
  
  return(good.ids)
}

# Return the protein names contained in a comma separated list
#
#
get.proteins <- function(s) unlist(strsplit(s, ",", fixed = TRUE))

#
#
#
extract.genes <- function(genes) return(unlist(sapply(genes, get.proteins)))

# Functions that report the interaction type: P-P, p-P or P-p, p-p
#
#
both.phospho       <- function(ds) return(ds$v1.phospho & ds$v2.phospho)
atleastone.phospho <- function(ds) return(xor(ds$v1.phospho, ds$v2.phospho))
none.phospho       <- function(ds) return(!ds$v1.phospho & !ds$v2.phospho)

#
#
#
print.percentage <- function(value) cat("Direct edges in iref: ", paste(format(value, digits=2),"%,"), "\nIndirect edges: ", paste(format(100-value, digits=3),"%,"), "\n")

# Compute the edges intersection between pairs and dataset
#
#
edges.intersection <- function(method.edges, ds){
  if(nrow(method.edges) == 0) return(c())
  
  esm <- create.adj.mtx(ds)
  
#   cat("\nStarting the comparison...")
  tmp <- sapply(1:nrow(method.edges), function(i){
#     if(i %% 1000 == 0) cat(" ...",i)
    gc1 <- extract.genes(as.character(method.edges[i,1]))
    gc2 <- extract.genes(as.character(method.edges[i,2]))
    
    for(j in 1:length(gc1)){
      for(k in 1:length(gc2)){
        if(gc1[j] %in% rownames(esm) && gc2[k] %in% rownames(esm) && esm[gc1[j], gc2[k]])
          # print(paste(gc1[j], gc2[k]))
          return(TRUE)
      }
    }
    return(FALSE)
  })
#   cat("\nend!")
  return(tmp)
}

#
#
#
load.dataset <- function(dname, ncol){
  ds <- read.csv(dname, check.names = FALSE)
#   rownames(ds) <- ds[,1]
  ds <- ds[,-c(1:ncol)]
  return(ds)
}

#
#
#
mds.plot <- function(expr){
  d <- dist(expr) # euclidean distances between the rows
  fit <- cmdscale(d, eig=TRUE, k=2) # k is the number of dim
  
  x <- fit$points[,1]
  y <- fit$points[,2]
  plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2", main="Metric  MDS",  type="n")
  text(x, y, labels = row.names(expr), cex=.7)
}

#
#
#
get.phospho.protein <- function(ns){
  return(grep("_p", ns))
}

#
#
#
unique.predictions <- function(tissue.predictions){
  temp.pred <- unique(tissue.predictions)
  
  spred     <- paste0(temp.pred[,1], "-", temp.pred[,2])
  rev.spred <- paste0(temp.pred[,2], "-", temp.pred[,1])
  
  fpred <- c()
  
  for( i in 1:length(spred)){
    if( !(spred[i] %in% fpred) && !(rev.spred[i] %in% fpred) ) fpred <- c(fpred, spred[i])
  }
  
  return(do.call(rbind, strsplit(fpred, "-", fixed = TRUE)))
}

#
#
#
get.best.predictions <- function(methods.predictions, tissue.name, method.name){
  if( !exists("methods.predictions") ) stop("Variable methods.predictions not found. Exit.\n")
  
  temp  <- methods.predictions[[ tissue.name ]]
  mname <- paste0(method.name, ".ds")
  
  if( !mname %in% names(temp)) cat("Warning: method not found in the predictions, return NULL.\n")
  return(temp[[ mname ]])
}

# Specific function to load TissueNet datasets
#
#
load.PPI <- function(f){
  ds <- read.delim(f, header = FALSE)
  colnames(ds) <- c("proteinA", "proteinB", "typeA", "typeB", "descr")
  return(ds)
}

# Replace symbol with string specified in the dataframe's 2nd column
#
#
convert <- function(x, map, type=NULL){
  if(any(duplicated(map[,1]))) stop("Conversion table is not a valid map. Exit.\n")
  
  if(is.null(type)){
    type <- 2
    rownames(map) <- map[,1]
  } else if(type == 1){
    rownames(map) <- map[,2]
  } else if(type == 2){
    rownames(map) <- map[,1]
  } else{
    warning("Invalid type for conversion. Set to default: 2nd column.\n")
    type <- 2
    rownames(map) <- map[,1]
  }
  
  return(map[as.character(x), type])
}

# Return a pair of elements of position (id) and (id+1) from vector v
#
#
#
elbind <- function(id, v){ return( c(v[id], v[id+1]) ) }

# Returns a matrix of interaction pairs from a path of interactions
#
#
#
get.predictions <- function(path){
  # get the elements of the path
  elems <- unlist(strsplit(path, ",", fixed = TRUE))
  
  # create interaction pairs
  do.call(rbind, lapply(1:(length(elems)-1), elbind, elems))
}

#
#
#
extract.gene.pair <- function(pairs){
  temp <- unlist(pairs)
  ids.rm <- temp == ""
  fpairs <- temp[!ids.rm]
  return(sort(unique(unlist(strsplit(fpairs, split = "-", fixed = TRUE)))))
}

#
#
#
get.matrix.edges <- function(m){
  
  res   <- cbind(which(m == TRUE, arr.ind=TRUE))
  edges <- as.data.frame(res)
  rownames(edges) <- 1:nrow(edges)
  colnames(edges) <- c("v1", "v2")
  
  return(net.from.index(tissue.net = edges, gene.names = colnames(m)))
}

#
#
#
percent <- function(first, last) return(last*100 / first)

# Compute n/p distribution
#
#
np.ratios <- function(ress){
  ratios <- lapply(ress, function(res){
    tissue.data <- res$data
    n <- nrow(tissue.data)
    p <- ncol(tissue.data)
    return(n/p)
  })
  
  plot(unlist(ratios), main = "n over p distribution", xlab = "tissues", ylab = "ratio") 
  abline(h=1)
}
