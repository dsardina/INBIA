setwd("D:/Dottorato/material/magi")
options(stringsAsFactors = FALSE)

require(ROCR)

###########################
## COMPUTE and LABEL PPIs
###########################

get.tissuenet.labels <- function(methods.predictions, cancer.types, best.methods, nname, tissuenet.dir, tissuenet.files, ctm){
  
  # Create map table for tissuenet conversion purposes
  conv.table    <- read.delim("tissuenet/tissuenet.map.txt")
  map           <- unique(conv.table[,2:3])
  rownames(map) <- map$gene_symbol
  rm(conv.table)
  
  # Compute all interactions and label them with
  # specific tissuenet normal counterpart (present/not present)
  cat("Compute all interactions...\n")
  all.genes <- as.character(map[,2])
  all.edges <- compute.all.interactions(all.genes)
  
  # Label all predictions for each tissue counterpart
  lapply(cancer.types, function(tissue){
    cat("Compute intersection for", tissue, "\n")
    # Bind together the predictions from paths
    tissue.predictions <- as.data.frame(all.edges)
    
    # Select current term for file selection
    term     <- ctm[tissue, 2]
    ds.files <- paste0(tissuenet.dir, tissuenet.files[grep(term, tissuenet.files, fixed = TRUE)])
    
    lapply(best.methods, function(mname){
      cat(" intersection with", mname, "\n")
      k <- 1
      best.preds <- get.best.predictions(methods.predictions, tissue, mname)
      preds.ids <- best.preds$path.length == (k+1)
      k.preds   <- best.preds[preds.ids, paste0(nname, ".shortest.path")]
      predictions <- do.call(rbind, lapply(k.preds, get.predictions))
      cbests <- apply(predictions, 2, convert, map)
      cbests <- as.data.frame(na.omit(cbests))
      
      # add column for current best method
      tissue.predictions[,mname] <<- as.numeric(edges.intersection(tissue.predictions, cbests))
    })
    
    # Compute score for current best methods
    cat(" compute ensembl score...\n")
    tissue.predictions$ensembl.score <- apply(tissue.predictions[, best.methods], 1, sum) / length(best.methods)
    
    # Compute predictions intersection with Tissue Net
    preds.tissuenet <- lapply(ds.files, function(ds){
      
      # Load current dataset from Tissue Net
      cmp.ds <- load.PPI(ds)
      
      # get only tcpa genes
      c1 <- cmp.ds[,1] %in% map$ensembl_gene
      c2 <- cmp.ds[,2] %in% map$ensembl_gene
      sub.ids <- c1 & c2
      cmp.ds <- cmp.ds[sub.ids, ]
      rm(c1, c2, sub.ids)
      
      curr.predictions <- tissue.predictions
      
      # Compute labels
      cat("  retrieving labels for current normal tissue...\n")
      curr.predictions$label <- as.numeric(edges.intersection(curr.predictions, cmp.ds))
      
      return(curr.predictions)
    })
    
    return(preds.tissuenet)
    
  })
}


###########################
## COMPUTE FOR INBIA
###########################

load("rdata/methods.predictions_wgcna_all.RData")
inbia.methods.predictions <- methods.predictions
names(inbia.methods.predictions) <- tissue.names
rm(methods.predictions)

inbia.bmr <- read.delim("bestanalysis/results.best.methods_induced.txt")
rownames(inbia.bmr) <- inbia.bmr$CancerType

tissuenet.ctm <- read.delim("comparison/cancer.tissuenet.map.txt")
rownames(tissuenet.ctm) <- tissuenet.ctm$cancer
tissuenet.dir   <- "tissuenet/hpa-protein/"
tissuenet.files <- list.files(tissuenet.dir)

inbia.best.methods <- unique(inbia.bmr$k.1)

## Compute path statistics
inbia.tnet <- get.tissuenet.labels(inbia.methods.predictions, tissue.names,
                                   inbia.best.methods, "iref", tissuenet.dir,
                                   tissuenet.files, tissuenet.ctm)
names(inbia.tnet) <- tissue.names


###########################
## COMPUTE FOR PERA
###########################

load("pera/methods.predictions.pera_all.RData")
pera.methods.predictions <- methods.predictions
names(pera.methods.predictions) <- tissue.names
rm(methods.predictions)

pera.bmr <- read.delim("pera/bestmethods/results.best.methods.txt")
rownames(pera.bmr) <- pera.bmr$CancerType

pera.best.methods <- unique(pera.bmr$k.1)

## Compute path statistics
pera.tnet <- get.tissuenet.labels(pera.methods.predictions, tissue.names,
                                  pera.best.methods, "pera", tissuenet.dir, tissuenet.files, tissuenet.ctm)
names(pera.tnet) <- tissue.names


##################
## ROCR ANALYSIS
##################

plot.rocr <- function(inbia.nets, pera.nets, tissue, type){
  
  lapply(1:length(inbia.nets[[tissue]]), function(i){
    inbia.ctissue <- inbia.nets[[tissue]][[i]]
    inbia.preds <- prediction(inbia.ctissue$ensembl.score, inbia.ctissue$label)
    
    inbia.roc.perf <- performance(inbia.preds, "tpr", "fpr")
    inbia.pr.perf <- performance(inbia.preds, "prec", "rec")
    
    pera.ctissue <- pera.nets[[tissue]][[i]]
    pera.preds <- prediction(pera.ctissue$ensembl.score, pera.ctissue$label)
    
    pera.roc.perf <- performance(pera.preds, "tpr", "fpr")
    pera.pr.perf <- performance(pera.preds, "prec", "rec")
    
    if(type == "PR"){
      plot(inbia.pr.perf, main=paste0("P-R - ", tissue), col="orange", xlim=c(0,1), ylim=c(0,1))
      plot(pera.pr.perf, add=TRUE, col="blue")
      legend("topright", legend = c("INBIA", "PERA"), lty = 1, col=c("orange", "blue"))
    } else if(type == "ROC"){
      plot(inbia.roc.perf, main=paste0("ROC - ", tissue), col="orange", xlim=c(0,1), ylim=c(0,1))
      plot(pera.roc.perf, add=TRUE, col="blue")
      abline(a=0, b= 1)
      legend("bottomright", legend = c("INBIA", "PERA"), lty = 1, col=c("orange", "blue"))
    }
  })
}

png("classification/pr.tissue.curve_all.png", width = 10, height = 12, units = "in", res = 600)
par(mfrow=c(5,4))
lapply(tissue.names, function(tissue){
  plot.rocr(inbia.tnet, pera.tnet, tissue, "PR")
})
dev.off()

png("classification/roc.tissue.curve_all.png", width = 10, height = 12, units = "in", res = 600)
par(mfrow=c(5,4))
lapply(tissue.names, function(tissue){
  plot.rocr(inbia.tnet, pera.tnet, tissue, "ROC")
})
dev.off()

