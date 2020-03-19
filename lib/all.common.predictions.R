setwd("D:/Dottorato/material/magi")
options(stringsAsFactors = FALSE)

load("rdata/methods.predictions_wgcna_all.RData") # load methods.predictions
inbia.methods.predictions <- methods.predictions
names(inbia.methods.predictions) <- tissue.names
rm(methods.predictions)

load("pera/methods.predictions.pera_all.RData")
pera.methods.predictions <- methods.predictions
names(pera.methods.predictions) <- tissue.names
rm(methods.predictions)

load("rdata/tcpa.proteins.RData") # load tcpa.genes
convTab <- load.tcpa.map()
cgenes <- extract.genes(unique(convTab[tcpa.genes, 2]))
rm(tcpa.genes)

# human.iref <- read.delim("IRefNet_New_HGNC.txt", header = FALSE)
# induced.human.iref <- subgraph(human.iref, cgenes)
# induced.iref.graph <- graph.edge(edges = as.data.frame(induced.human.iref))


# adj.matrices <- lapply(methods.predictions, function(tissue.preds){
#   cat("Tissue predictions\n")
#   lapply(tissue.preds, function(mpreds){
#     cat("...creating adjacency matrix\n")
#     create.adj.mtx(mpreds, ugenes = cgenes)
#   })
# })

# lapply(methods.predictions, function(tissue.preds){
#   cat("Tissue predictions\n")
#   # bmethods <- c("clr.ds", "glasso.ds", "pls.ds", "mrnet.ds")
#   bmethods <- c("elasticnet.ds", "pls.ds", "aracnea.ds", "spearman.ds", "wgcna.ds", "clr.ds")
#   lapply(bmethods, function(mname){
#     preds <- get.phosphorilated(tissue.preds[[mname]])
#     
#     cat("...creating adjacency matrix for", mname, "\n")
#     create.adj.mtx(preds, ugenes = agenes)
#   })
# }) 

get.phosphorilated <- function(ds){
  ph <- ds$v1.phospho
  ds$v1[ph] <- paste0(ds$v1[ph], "_P")
  ph <- ds$v2.phospho
  ds$v2[ph] <- paste0(ds$v2[ph], "_P")
  return(ds)
}

genes <- convTab$V2
# ph <- is.phospho(convTab$V1)
# genes[ph] <- paste0(genes[ph], "_P")
# genes <- unique(genes)

adj.matrices <- function(methods.predictions, genes, filter=NULL, best.methods=NULL){
  lapply(1:length(methods.predictions), function(i){
    tissue.preds <- methods.predictions[[i]]
    cat("Tissue predictions for", tissue.names[i], "\n")
    
    # if no filter, select all methods
    if(is.null(filter)) filter <- names(tissue.preds)
    if(!is.null(best.methods)) filter <- paste0(best.methods[i,2], ".ds")
    
    if(!is.null(best.methods)){
      cat("...creating adjacency matrix for", filter, "\n")
      return(create.adj.mtx(tissue.preds[[filter]], ugenes = genes))
    }
    
    every <- lapply(filter, function(mname){
      cat("  creating adjacency matrix for", mname, "\n")
      create.adj.mtx(tissue.preds[[mname]], ugenes = genes)
    })
    
  }) 
}

common.predictions <- function(adj.matrices){
  lapply(adj.matrices, function(adj){
    res <- adj[[1]]
    
    for(i in 2:length(adj)){
      res <- res & adj[[i]]
    }
    
    return(unique.predictions(get.matrix.edges(res)))
  }) 
}

tissue.common.predictions <- function(adj.matrices){
  res <- adj.matrices[[1]]
  
  for(i in 2:length(adj.matrices)){
    res <- res & adj.matrices[[i]]
  }
  
  unique.predictions(get.matrix.edges(res))
}


#################################################
## COMPUTE ALL METHODS INTERSECTION BY TISSUE
#################################################

inbia.adj.matrices <- adj.matrices(inbia.methods.predictions, genes) # no filter, no best methods
pera.adj.matrices <- adj.matrices(pera.methods.predictions, genes) # no filter, no best methods

inbia.tc <- common.predictions(inbia.adj.matrices)
pera.tc <- common.predictions(pera.adj.matrices)

sapply(inbia.tc, length)
sapply(pera.tc, length)

# FUNCTIONAL ANALYSIS
"functional/methods.common/inbia"
"functional/methods.common/pera/"
save.tissue.genes <- function(methods.common, dir.name){
  lapply(1:length(tissue.names), function(i){
    tc <- methods.common[[i]]
    ug <- unique.genes(tc)
    write.table(ug, file = paste0(dir.name, "tc.", tissue.names[i], ".txt"), quote=FALSE,
                row.names = FALSE, col.names = FALSE)
  })
}

save.tissue.genes(inbia.tc, "functional/methods.common/inbia/")
save.tissue.genes(pera.tc, "functional/methods.common/pera/")

#################################################
## COMPUTE ALL TISSUE BEST METHODS INTERSECTION
#################################################

inbia.bmr <- read.delim("bestanalysis/results.best.methods_induced.txt")
pera.bmr <- read.delim("pera/bestmethods/results.best.methods.txt")

inbia.best.adj.matrices <- adj.matrices(inbia.methods.predictions, genes, best.methods = inbia.bmr) # no filter
pera.best.adj.matrices <- adj.matrices(pera.methods.predictions, genes, best.methods = pera.bmr) # no filter

inbia.best.tc <- tissue.common.predictions(inbia.best.adj.matrices)
pera.best.tc <- tissue.common.predictions(pera.best.adj.matrices)

# write.table(inbia.best.tc, "comparison/allcomp/inbia.allcancers.edges.txt", quote=FALSE, row.names = FALSE, col.names = FALSE, sep="\t")
# write.table(pera.best.tc, "comparison/allcomp/pera.allcancers.edges.txt", quote=FALSE, row.names = FALSE, col.names = FALSE, sep="\t")

## Extract gene for functional analysis
inbia.allcancer.genes <- unique.genes(inbia.best.tc)
pera.allcancer.genes <- unique.genes(pera.best.tc)

write.table(inbia.allcancer.genes, "functional/allcancer.common/inbia.allcancer.genes.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(pera.allcancer.genes, "functional/allcancer.common/pera.allcancer.genes.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

source("scripts/mds.R")
source("scripts/collective.influence.R")

graph.commons <- function(tcp){
  cat("Computing graph statistics...\n")
  g <- graph.edge(tcp)
  cmds <- mds(g)
  cci   <- collective.influence(g)
  return(list(common.graph=g, common.mds=cmds, common.ci=cci))
}

inbia.common.graph <- graph.commons(inbia.best.tc)
pera.common.graph <- graph.commons(pera.best.tc)

plot.common.graph <- function(g.commons){
  g       <- g.commons$common.graph
  vmds    <- rep("orange", length(V(g)))
  names(vmds) <- names(V(g))
  
  mds.ids <- g.commons$common.mds$mds
  vmds[mds.ids] <- "red"
  
  plot.graph(g, g.commons$common.ci$ci, vcolor = vmds)
}

png(filename = "comparison/allcomp/all.cancer.plot/inbia.graph.png", width = 8, height = 8, unit = "in", res = 600)
plot.common.graph(inbia.common.graph)
dev.off()

png(filename = "comparison/allcomp/all.cancer.plot/pera.graph.png", width = 8, height = 8, unit = "in", res = 600)
plot.common.graph(pera.common.graph)
dev.off()

#################################################################
## Compute the quality of common predictions by using IRefIndex
common.intersection <- lapply(tissue.common.predictions, function(preds){
  cat("Computing dataset for tissue\n")
  common.ds <- as.data.frame(preds)
  
  common.ds$present.iref       <- edges.intersection(preds, induced.human.iref)
  common.ds$iref.shortest.path <- compute.shortest.paths(preds, induced.iref.graph)
  common.ds$path.length        <- get.path.length.distr(common.ds$iref.shortest.path)
  
  return(common.ds)
})


#########################
## Compute True Positive
#########################

tp.values <- lapply(common.intersection, function(predictions){
  lapply(2:5, function(k){
    sum(predictions$path.length <= k)
  })
})

tp <- lapply(tp.values, function(vals){
  do.call(cbind, vals)
})

tp.complete <- as.data.frame(do.call(cbind, tp))


##########################
## Compute False Positive
##########################

total.predictions <- lapply(tissue.common.predictions, function(prs){ nrow(prs) })
total.prediction.ds <- do.call(cbind, total.predictions)

fp.values <- vector(mode = "list", length = length(common.intersection))

for(i in 1:length(common.intersection)){
  tprs <- total.prediction.ds[,i]
  ktp  <- tp[[i]]
  fp.values[[i]] <- tprs - ktp
}

fp.complete <- as.data.frame(do.call(cbind, fp.values))


##########################
## Compute True Negative
##########################

iref.genes_ind     <- unique.genes(induced.human.iref)
all.iref.edges_ind <- compute.all.interactions(iref.genes_ind)
iref.allsp_ind     <- compute.shortest.paths(all.iref.edges_ind, induced.iref.graph)

all.paths.iref_ind           <- as.data.frame(cbind(all.iref.edges_ind, iref.allsp_ind, get.path.length.distr(iref.allsp_ind)))
colnames(all.paths.iref_ind) <- c("node1", "node2", "shortest.path", "path.length")

library(parallel)

cl <- makeCluster(detectCores()-1)

## Compute TN in all tissue for all methods
tn.values <- vector(mode = "list", length = length(common.intersection))

for(i in 1:length(common.intersection)){
  cat(i,"...\n")
  predictions <- common.intersection[[i]]
  
  tn.values[[i]] <- lapply(2:5, function(k){
    cat("   k =",k, "...\n")
    iref.paths <- all.paths.iref_ind[as.numeric(all.paths.iref_ind$path.length) <= k, ]
    return(get.tn(iref.paths, predictions, iref.genes_ind, cl))
  })
}

stopCluster(cl)
rm(cl)

tn <- lapply(tn.values, function(vals){
  do.call(cbind, vals)
})

tn.complete <- as.data.frame(do.call(cbind, tn))

##########################
## Compute False Negative
##########################

cl <- makeCluster(detectCores()-1)

## Compute FN in all tissue for all methods
fn.values <- vector(mode = "list", length = length(common.intersection))

for(i in 1:length(common.intersection)){
  cat(i,"...\n")
  predictions <- common.intersection[[i]]
  
  fn.values[[i]] <- lapply(2:5, function(k){
    cat("   k =",k, "...\n")
    iref.paths <- all.paths.iref_ind[as.numeric(all.paths.iref_ind$path.length) <= k, ]
    
    return(get.fn(iref.paths, predictions, iref.genes_ind, cl))
  })
}

stopCluster(cl)
rm(cl)

fn <- lapply(fn.values, function(vals){
  do.call(cbind, vals)
})

fn.complete <- as.data.frame(do.call(cbind, fn))

save(tp.complete, fp.complete,
     tn.complete, fn.complete, file="comparison/allcomp/best.common.predictions.classification.results_induced.RData")

#############################
## Classification results
#############################
accuracy    <- function(tp, fp, tn, fn) return((tp + tn)/(tp + fp + tn + fn))
# error       <- function(accs) return(1- accs)
precision   <- function(tp, fp) return(tp/(tp + fp))
sensitivity <- function(tp, fn) return(tp / (tp + fn)) # also called recall
specificity <- function(tn, fp) return(tn/(tn + fp))
fmeasure    <- function(precision, recall) return(2*precision*recall/(precision+recall))

get.best.accuracy <- function(id){
  tmp <- (id*4)
  ids <- (tmp-3):tmp
  return( t(acc)[ids,-1] )
}

get.best.precision <- function(id){
  tmp <- (id*4)
  ids <- (tmp-3):tmp
  return( t(prec)[ids,-1] )
}

get.best.sensitivity <- function(id){
  tmp <- (id*4)
  ids <- (tmp-3):tmp
  return( t(sens)[ids,-1] )
}

get.best.fmeasure <- function(id){
  tmp <- (id*4)
  ids <- (tmp-3):tmp
  return( t(fmeas)[ids,-1] )
}

## Compute classification results

kcol  <- rep(c(2:5), times=16)
acc   <- rbind( kcol, accuracy(tp.complete, fp.complete, tn.complete, fn.complete))
prec  <- rbind( kcol, precision(tp.complete, fp.complete))
sens  <- rbind( kcol, sensitivity(tp.complete, fn.complete))
spec  <- rbind( kcol, specificity(tn.complete, fp.complete))
fmeas <- fmeasure(prec, sens)

dres   <- "comparison/allcomp/"
fname  <- paste0(dres, "all_accuracy.txt")

for(i in 1:length(tissue.names)){
  cat(tissue.names[i],"\n", file = fname, append = TRUE)
  ba <- get.best.accuracy(i)
  cat(paste0(ba, collapse = " & "), "\n", file = fname, append = TRUE)
}


fname <- paste0(dres, "all_precision.txt")

for(i in 1:length(tissue.names)){
  cat(tissue.names[i],"\n", file = fname, append = TRUE)
  ba <- get.best.precision(i)
  cat(paste0(ba, collapse = " & "), "\n", file = fname, append = TRUE)
}


fname <- paste0(dres, "all_sensitivity.txt")

for(i in 1:length(tissue.names)){
  cat(tissue.names[i],"\n", file = fname, append = TRUE)
  ba <- get.best.sensitivity(i)
  cat(paste0(ba, collapse = " & "), "\n", file = fname, append = TRUE)
}


fname <- paste0(dres, "all_fmeasure.txt")

for(i in 1:length(tissue.names)){
  cat(tissue.names[i],"\n", file = fname, append = TRUE)
  ba <- get.best.fmeasure(i)
  cat(paste0(ba, collapse = " & "), "\n", file = fname, append = TRUE)
}

rm(tp.complete, tn.complete, fn.complete, fp.complete)
rm(kcol, acc, prec, sens, spec, fmeas)

########################
## Script and results
########################

tissue.nets <- lapply(tissue.names, function(tissue){
  
  preds <- tissue.common.predictions[[tissue]]
  cat("Compute tissue graph\n")
  uproteins <- unique.genes(preds)
  pcolors   <- rep(1, length(uproteins))
  pcolors[ endsWith(uproteins, suffix = "_P") ] <- 2
  names(pcolors) <- uproteins
  
  if(is.null(preds)){
    kgraph <- NULL
    cat("Prediction graph is empty, skip to next tissue...\n")
  }
  else{
    kgraph <- graph.edge(preds)
    write.motif.network(g = kgraph, vcolor = pcolors[names(V(kgraph))], file = paste0("motif/4common/", tissue, ".txt"))
  }
})

tissue.edges <- lapply(tissue.names, function(tissue){
  preds <- tissue.common.predictions[[tissue]]
  write.table(preds, paste0("motif/pera6common/", tissue, "_edges.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
})

# tcp <- data.frame(tissue=names(tissue.common.predictions))
# 
# for(i in 1:length(tissue.common.predictions)){
#   temp <- tissue.common.predictions[[i]]
#   tcp$n[i] <- nrow(temp)
#   tcp$proteins[i] <- paste(unique.genes(temp), collapse = ",")
#   tcp$predictions[i] <- paste0(temp[,1], "-", temp[,2], collapse = ", ")
# }
# write.table(tcp, file="comparison/all.common.predictions.tsv", quote = FALSE, row.names = FALSE, sep = "\t")

# for(i in 1:length(tissue.common.predictions)){
#   temp <- tissue.common.predictions[[i]]
#   cat(names(tissue.common.predictions[i]), "\n")
#   cat(paste0("  ", unique.genes(temp), collapse = "\n"), "\n\n")
# }

g.commons <- lapply(tissue.common.predictions, function(tcp){
  cat("Computing graph statistics...\n")
  g <- graph.edge(tcp)
  cmds <- mds(g)
  cci   <- collective.influence(g)
  list(common.graph=g, common.mds=cmds, common.ci=cci)
})

names(g.commons) <- tissue.names
# tname <- tissue.names[1]

lapply(tissue.names, function(tname){
  g       <- g.commons[[tname]]$common.graph
  vmds    <- rep("orange", length(V(g)))
  names(vmds) <- names(V(g))
  
  mds.ids <- g.commons[[tname]]$common.mds$mds
  vmds[mds.ids] <- "red"
  
  cat("Save graph for", tname, "\n")
  on.exit(dev.off())
  png(filename = paste0("comparison/allcomp/all.graph.plot/", tname, "_common.graph.png"), width = 8, height = 8, unit = "in", res = 600)
  plot.graph(g, g.commons[[tname]]$common.ci$ci, vcolor = vmds)
})

