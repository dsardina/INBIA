source("lib/algorithms/mds.R")
source("lib/algorithms/collective.influence.R")

# Best method results (see classification results)
# bmr <- read.delim("results/classification/induced/best.methods_induced.txt")
# rownames(bmr) <- bmr[,1]

get.best.networks <- function(bmr){
  cat("### SELECTING BEST NETWORKS ###\n")
  dir.create("results/bestnetworks", recursive = TRUE, showWarnings = FALSE)
  
  cat("k=", k)
  for(k in 1:4){
    cat(",", k)
    tissue.nets <- lapply(bmr$TissueType, function(tissue){
      
      cat("..Tissue:", tissue, "\n")
      mname      <- bmr[tissue, k+1] # extract best method for current tissue and k
      best.preds <- get.best.predictions(methods.predictions, tissue, mname)
      path.ids   <- best.preds$path.length == (k+1)
      kpath      <- best.preds$gs.shortest.path[path.ids]
      
      kpath.matrix <- do.call(rbind, strsplit(kpath, ",", fixed = TRUE))
      
      if(is.null(kpath.matrix)) kgraph <- NULL
      else kgraph <- graph.matrix(create.kedge.adj.mtx(kpath.matrix))
      
      if(is.null(kgraph)){
        kmds <- NULL
        ci   <- NULL
      }
      else{
        kmds <- mds(kgraph)
        ci   <- collective.influence(kgraph)
      }
      
      return(list(k=k, tissue.name=tissue, tissue.mds=kmds,
                  tissue.ci=ci, tissue.kgraph=kgraph))
    })
    save(tissue.nets, file=paste0("results/bestnetworks/tissuenet_k", k, ".RData"))
    cat("\n")
  }
}


## MDS WITH COLLECTIVE INFLUENCE
for(k in 1:4){
  cat("k=", k, "\n")
  tissue.nets <- lapply(bmr$CancerType, function(tissue){
    
    cat("  Tissue:", tissue, "\n")
    mname         <- bmr[tissue, k+1] # extract best method for current tissue and k
    mname_induced <- bmr_induced[tissue, k+1] # extract best method for current tissue and k
    
    best.preds_induced <- get.best.predictions(tissue, mname_induced)
    best.preds_human   <- get.best.predictions(tissue, mname)
    
    path.ids       <- best.preds_induced$path.length == (k+1)
    human.path.ids <- best.preds_human$path.length.human == (k+1)
    
    kpath       <- best.preds_induced$iref.shortest.path[path.ids]
    human.kpath <- best.preds_human$human.iref.shortest.path[human.path.ids]
    
    kpath.matrix       <- do.call(rbind, strsplit(kpath, ",", fixed = TRUE))
    human.kpath.matrix <- do.call(rbind, strsplit(human.kpath, ",", fixed = TRUE))
    
    if(is.null(kpath.matrix)) kgraph <- NULL
    else kgraph <- graph.matrix(create.kedge.adj.mtx(kpath.matrix))
    
    if(is.null(human.kpath.matrix)) human.kgraph <- NULL
    else human.kgraph <- graph.matrix(create.kedge.adj.mtx(human.kpath.matrix))
    
    if(is.null(kgraph)){
      kmds <- NULL
      ci   <- NULL
    }
    else{
      ci   <- collective.influence(kgraph)
      kmds <- mds(kgraph, ci$ci)
    }
    
    if(is.null(human.kgraph)){
      human.kmds <- NULL
      human.ci   <- NULL
    }
    else{
      human.ci   <- collective.influence(human.kgraph)
      human.kmds <- mds(human.kgraph, human.ci$ci)
    }
    
    return(list(k=k, tissue.name=tissue, tissue.mds=kmds, human.tissue.mds=human.kmds, tissue.ci=ci,
                tissue.human.ci=human.ci, tissue.kgraph=kgraph, human.tissue.kgraph=human.kgraph))
  })
  save(tissue.nets, file=paste0("bestanalysis/tissuenet_CI_k", k, ".RData"))
}


load("bestanalysis/tissuenet_k1.RData")

# Extract MDS proteins from tissue
tmds.tissue.proteins <- lapply(tissue.nets, function(tnet){
  g <- tnet$tissue.kgraph
  mds.proteins <- V(g)[tnet$tissue.mds$mds]
  return(names(mds.proteins))
})

# Extract MDS proteins from tissue within the whole human protein interaction network
human.tmds.tissue.proteins <- lapply(tissue.nets, function(tnet){
  g <- tnet$human.tissue.kgraph
  mds.proteins <- V(g)[tnet$tissue.mds$mds]
  return(names(mds.proteins))
})

# Count the MDS proteins for each tissue
tmds.proteins       <- table(unlist(tmds.tissue.proteins))
human.tmds.proteins <- table(unlist(human.tmds.tissue.proteins))

# Compute the tissue-specific and house-keeping proteins
ts.mds <- tmds.proteins[tmds.proteins <= 3]
hk.mds <- tmds.proteins[tmds.proteins >= 11]

lapply(tmds.tissue.proteins, function(tpr){
  tpr[tpr %in% names(ts.mds)]
})

################################################################################
## Example                                                                    ##
################################################################################
names(tissue.nets) <- tissue.names
tname <- tissue.names[2]
g       <- tissue.nets[[tname]]$tissue.kgraph
vmds    <- rep("orange", length(V(g)))
mds.ids <- tissue.nets[[tname]]$tissue.mds$mds

vmds[mds.ids]                       <- "red"
# vmds[which(names(V(g)) == "YBX1")] <- "green"

plot.graph(g, tissue.nets[[tname]]$tissue.ci$ci, vmds)
################################################################################

## Clustering of collective influence
tree   <- hclust(dist(tissue.nets[[tname]]$tissue.ci$ci), method = "average")
groups <- cutree(tree, k = 2)
plot(tree)

# Plotting
mds.pgraph(kgraph, kres$visited, kres$mds)
plot.graph(kgraph, ci$ci)
ci.pgraph(kgraph, ci$balls$TP53, V(kgraph)["TP53"])

lmatrix <- laplacian_matrix(kgraph)
eigen.laplacian <- eigen(data.matrix(lmatrix), symmetric = TRUE)
