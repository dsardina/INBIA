source("lib/method.functions.R")
source("conf.R")

require(igraph)
require(parallel)

#########################
## Compute True Positive
#########################

compute.tp.values <- function(methods.predictions){
  tp.values <- vector(mode = "list", length = length(methods.predictions))
  
  for(i in 1:length(methods.predictions)){
    predictions <- methods.predictions[[i]]
    
    tp.values[[i]] <- lapply(2:5, function(k){
      c(
        k = k,
        pearson.TP    = sum(predictions$pearson.ds$path.length <= k, na.rm = TRUE),
        spearman.TP   = sum(predictions$spearman.ds$path.length <= k, na.rm = TRUE),
        spc.TP        = sum(predictions$spc.ds$path.length <= k, na.rm = TRUE),
        genenet.TP    = sum(predictions$genenet.ds$path.length <= k, na.rm = TRUE),
        glasso.TP     = sum(predictions$glasso.ds$path.length <= k, na.rm = TRUE),
        pls.TP        = sum(predictions$pls.ds$path.length <= k, na.rm = TRUE),
        ridgenet.TP   = sum(predictions$ridge.ds$path.length <= k, na.rm = TRUE),
        lasso.TP      = sum(predictions$lasso.ds$path.length <= k, na.rm = TRUE),
        elasticnet.TP = sum(predictions$elasticnet.ds$path.length <= k, na.rm = TRUE),
        aracnea.TP    = sum(predictions$aracnea.ds$path.length <= k, na.rm = TRUE),
        aracnem.TP    = sum(predictions$aracnem.ds$path.length <= k, na.rm = TRUE),
        clr.TP        = sum(predictions$clr.ds$path.length <= k, na.rm = TRUE),
        mrnet.TP      = sum(predictions$mrnet.ds$path.length <= k, na.rm = TRUE),
        wgcna.TP      = sum(predictions$wgcna.ds$path.length <= k, na.rm = TRUE)
      )
    })
  }
  
  tp <- lapply(tp.values, function(vals){
    do.call(cbind, vals)
  })
  
  tp.complete_induced               <- as.data.frame(do.call(cbind, tp))
  rownames(tp.complete_induced)[-1] <- methods.names
  dir.create("results/classification/induced", showWarnings = FALSE, recursive = TRUE)
  save(tp.complete_induced, file="results/classification/induced/true.positive_induced.RData")
  # write.table(tp.complete_induced, "results/classification/induced/true_positive_induced.csv", quote = FALSE, sep = ";", col.names = FALSE)
  return(tp.complete_induced)
}


#############################
## Collect all predictions
#############################

compute.total.preds <- function(methods.intersections){
  total.predictions <- lapply(methods.intersections, function(mint){
    c(
      pearson.total.predictions    = nrow(mint$annotated.pearson),
      spearman.total.predictions   = nrow(mint$annotated.spearman),
      spc.total.predictions        = nrow(mint$annotated.spc),
      genenet.total.predictions    = nrow(mint$annotated.genenet),
      glasso.total.predictions     = nrow(mint$annotated.glasso),
      pls.total.predictions        = nrow(mint$annotated.pls),
      ridgenet.total.predictions   = nrow(mint$annotated.ridgenet),
      lasso.total.predictions      = nrow(mint$annotated.lasso),
      elasticnet.total.predictions = nrow(mint$annotated.elasticnet),
      aracnea.total.predictions    = nrow(mint$annotated.aracnea),
      aracnem.total.predictions    = nrow(mint$annotated.aracnem),
      clr.total.predictions        = nrow(mint$annotated.clr),
      mrnet.total.predictions      = nrow(mint$annotated.mrnet),
      wgcna.total.predictions      = nrow(mint$annotated.wgcna)
    )
  })
  
  total.prediction.ds           <- do.call(cbind, total.predictions)
  rownames(total.prediction.ds) <- methods.names
  
  return(total.prediction.ds)
}


##########################
## Compute False Positive
##########################

compute.fp.values <- function(methods.predictions, methods.intersections, tp.values){
  total.prediction.ds <- compute.total.preds(methods.intersections)
  
  fp.values <- vector(mode = "list", length = length(methods.predictions))
  tp.values <- tp.values[-1,]
  
  if(nrow(tp.values) != nrow(total.prediction.ds)) stop("Cannot compute false positives: different number of methods.\n")
  
  ids <- 1:4
  for(i in 1:length(methods.predictions)){
    fp.values[[i]] <- total.prediction.ds[,i] - tp.values[,ids]
    ids <- ids +4
  }
  
  fp.complete_induced              <- as.data.frame(do.call(cbind, fp.values))
  fp.complete_induced              <- rbind(rep(2:5, length(methods.predictions)), fp.complete_induced)
  rownames(fp.complete_induced)[1] <- "k"
  dir.create("results/classification/induced", showWarnings = FALSE, recursive = TRUE)
  save(fp.complete_induced, file="results/classification/induced/false.positive_induced.RData")
  # write.table(fp.complete_induced, "results/classification/induced/false.positive_induced.csv", quote = FALSE, sep = ";", col.names = FALSE)
  return(fp.complete_induced)
}


##########################################################
## Compute all shortest path for the induced GS graph
##########################################################
compute.gs.shortest.paths <- function(induced.human.gs){
  cat("### COMPUTING ALL SHORTEST PATHS WITHIN GS ###\n")
  cat("..Extracting genes\n")
  gs.genes_ind     <- unique.genes(induced.human.gs)
  induced.gs.graph <- graph.edge(edges = as.data.frame(induced.human.gs))
  cat("..Computing all interactions\n")
  all.gs.edges_ind <- compute.all.interactions(gs.genes_ind)
  cat("..Computing all shortest paths\n")
  gs.allsp_ind     <- compute.shortest.paths(all.gs.edges_ind, induced.gs.graph)

  all.paths.gs_ind <- as.data.frame(cbind(all.gs.edges_ind, gs.allsp_ind, get.path.length.distr(gs.allsp_ind)))
  colnames(all.paths.gs_ind) <- c("node1", "node2", "shortest.path", "path.length")
  return(all.paths.gs_ind)
}


###############################################################
## Compute True Negative
## ============================================================
## Collect all predictions for current method and tissue.
## Consider the prediction (A,B), it must not exist a path of
## length <= k in the GS, for each k=2,3,4,5.
###############################################################

compute.tn.values <- function(methods.predictions, induced.human.gs, all.paths.gs_ind){
  cat("### COMPUTING TRUE NEGATIVE VALUES ###\n")
  cat("..Creating cluster for parallel computation\n")
  cl <- makeCluster(detectCores()-1)
  on.exit(stopCluster(cl))
  
  ## Compute TN in all tissue for all methods
  tn.values <- vector(mode = "list", length = length(methods.predictions))
  gs.genes_ind <- unique.genes(induced.human.gs)
  
  cat("..Computing for each element and k path length\n")
  for(i in 1:length(methods.predictions)){
    cat("..", i)
    predictions <- methods.predictions[[i]]

    cat(" ... k =")
    tn.values[[i]] <- lapply(2:5, function(k){
      cat(k, ",")
      gs.paths <- all.paths.gs_ind[as.numeric(all.paths.gs_ind$path.length) <= k, ]
      
      c(
        k = k,
        pearson.TN    = get.tn(gs.paths, predictions$pearson.ds, gs.genes_ind, cl),
        spearman.TN   = get.tn(gs.paths, predictions$spearman.ds, gs.genes_ind, cl),
        spc.TN        = get.tn(gs.paths, predictions$spc.ds, gs.genes_ind, cl),
        genenet.TN    = get.tn(gs.paths, predictions$genenet.ds, gs.genes_ind, cl),
        glasso.TN     = get.tn(gs.paths, predictions$glasso.ds, gs.genes_ind, cl),
        pls.TN        = get.tn(gs.paths, predictions$pls.ds, gs.genes_ind, cl),
        ridgenet.TN   = get.tn(gs.paths, predictions$ridge.ds, gs.genes_ind, cl),
        lasso.TN      = get.tn(gs.paths, predictions$lasso.ds, gs.genes_ind, cl),
        elasticnet.TN = get.tn(gs.paths, predictions$elasticnet.ds, gs.genes_ind, cl),
        aracnea.TN    = get.tn(gs.paths, predictions$aracnea.ds, gs.genes_ind, cl),
        aracnem.TN    = get.tn(gs.paths, predictions$aracnem.ds, gs.genes_ind, cl),
        clr.TN        = get.tn(gs.paths, predictions$clr.ds, gs.genes_ind, cl),
        mrnet.TN      = get.tn(gs.paths, predictions$mrnet.ds, gs.genes_ind, cl),
        wgcna.TN      = get.tn(gs.paths, predictions$wgcna.ds, gs.genes_ind, cl)
      )
    })
    cat("\n")
  }
  
  tn <- lapply(tn.values, function(vals){
    do.call(cbind, vals)
  })
  
  tn.complete_induced               <- as.data.frame(do.call(cbind, tn))
  rownames(tn.complete_induced)[-1] <- methods.names
  dir.create("results/classification/induced", showWarnings = FALSE, recursive = TRUE)
  save(tn.complete_induced, file = "results/classification/induced/true.negative_induced.RData")
  # write.table(tn.complete_induced, "results/classification/induced/true.negative_induced.csv", quote = FALSE, sep = ";", col.names = FALSE)
  return(tn.complete_induced)
}


#####################################################################
## Compute False Negative - computational method
## ==================================================================
## Collect all unpredicted interactions and check how many of them
## are within a path of length <= k in GS, for each k=2,3,4,5.
#####################################################################

compute.fn.values <- function(methods.predictions, induced.human.gs, all.paths.gs_ind){
  cat("### COMPUTING FALSE NEGATIVE VALUES ###\n")
  cat("..Creating cluster for parallel computation\n")
  cl <- makeCluster(detectCores()-1)
  on.exit(stopCluster(cl))
  
  ## Compute FN in all tissue for all methods
  fn.values <- vector(mode = "list", length = length(methods.predictions))
  gs.genes_ind <- unique.genes(induced.human.gs)
  
  cat("..Computing for each element and k path length\n")
  for(i in 1:length(methods.predictions)){
    cat("..", i)
    predictions <- methods.predictions[[i]]
    cat(" ... k =")
    fn.values[[i]] <- lapply(2:5, function(k){
      cat(k, ",")
      gs.paths <- all.paths.gs_ind[as.numeric(all.paths.gs_ind$path.length) <= k, ]
      
      c(
        k = k,
        pearson.FN    = get.fn(gs.paths, predictions$pearson.ds, gs.genes_ind, cl),
        spearman.FN   = get.fn(gs.paths, predictions$spearman.ds, gs.genes_ind, cl),
        spc.FN        = get.fn(gs.paths, predictions$spc.ds, gs.genes_ind, cl),
        genenet.FN    = get.fn(gs.paths, predictions$genenet.ds, gs.genes_ind, cl),
        glasso.FN     = get.fn(gs.paths, predictions$glasso.ds, gs.genes_ind, cl),
        pls.FN        = get.fn(gs.paths, predictions$pls.ds, gs.genes_ind, cl),
        ridgenet.FN   = get.fn(gs.paths, predictions$ridge.ds, gs.genes_ind, cl),
        lasso.FN      = get.fn(gs.paths, predictions$lasso.ds, gs.genes_ind, cl),
        elasticnet.FN = get.fn(gs.paths, predictions$elasticnet.ds, gs.genes_ind, cl),
        aracnea.FN    = get.fn(gs.paths, predictions$aracnea.ds, gs.genes_ind, cl),
        aracnem.FN    = get.fn(gs.paths, predictions$aracnem.ds, gs.genes_ind, cl),
        clr.FN        = get.fn(gs.paths, predictions$clr.ds, gs.genes_ind, cl),
        mrnet.FN      = get.fn(gs.paths, predictions$mrnet.ds, gs.genes_ind, cl),
        wgcna.FN      = get.fn(gs.paths, predictions$wgcna.ds, gs.genes_ind, cl)
      )
    })
    cat("\n")
  }
  
  fn <- lapply(fn.values, function(vals){
    do.call(cbind, vals)
  })
  
  fn.complete_induced <- as.data.frame(do.call(cbind, fn))
  rownames(fn.complete_induced)[-1] <- methods.names
  dir.create("results/classification/induced", showWarnings = FALSE, recursive = TRUE)
  save(fn.complete_induced, file = "results/classification/induced/false.negative_induced.RData")
  # write.table(fn.complete, "results/classification/induced/false.negative_induced.csv", quote = FALSE, sep = ";", col.names = FALSE)
  return(fn.complete_induced)
}

########################
## Save all TP FP TN FN
########################
# save(tp.complete_induced, fp.complete_induced,
     # tn.complete_induced, fn.complete_induced, file="classification/wgcna-all/induced/classification.results_induced.RData")
