source("conf.R")

save.results <- function(ds, dir, file, measure){
  fname  <- paste0(dir, file)
  res <- as.data.frame(do.call(rbind, ds))
  colnames(res) <- c("method", measure, "tissue", "k")

  write.table(res, fname, quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
}

build.best.methods <- function(ds, dir, file, tissue.names){
  res <- as.data.frame(do.call(rbind, ds))
  rownames(res) <- 1:nrow(res)
  colnames(res) <- c("method", "measure", "tissue", "k")
  
  fname  <- paste0(dir, file)
  bm.res <- lapply(tissue.names, function(tissue){
    temp <- res[res$tissue == tissue,]
    rownames(temp) <- temp$k
    names(methods.dsnames) <- methods.names
    return(c(tissue, methods.dsnames[temp$method]))
  })
  bm <- as.data.frame(do.call(rbind, bm.res))
  colnames(bm) <- c("TissueType", "k=1", "k=2", "k=3", "k=4")
  write.table(bm, fname, quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
}

#########################
## Function definitions
#########################

accuracy    <- function(tp, fp, tn, fn) return((tp + tn)/(tp + fp + tn + fn))
error       <- function(accs) return(1- accs)
precision   <- function(tp, fp) return(tp/(tp + fp))
sensitivity <- function(tp, fn) return(tp / (tp + fn)) # also called recall
specificity <- function(tn, fp) return(tn/(tn + fp))
fmeasure    <- function(precision, recall) return(2*precision*recall/(precision+recall))

get.sample.measures <- function(id, mtype){
  tmp <- (id*4)
  ids <- (tmp-3):tmp
  return( t(mtype)[ids,-1] )
}

getbestaccuracy <- function(i, acc){
  ba <- get.sample.measures(i, acc)
  tmp.method   <- apply(ba, 1, which.max)
  tmp.k        <- apply(ba, 1, max, na.rm=TRUE)
  names(tmp.k) <- methods.names[tmp.method]
  return(tmp.k)
}

getbestprecision <- function(i, prec){
  bp <- get.sample.measures(i, prec)
  tmp.method   <- apply(bp, 1, which.max)
  tmp.k        <- apply(bp, 1, max, na.rm=TRUE)
  names(tmp.k) <- methods.names[tmp.method]
  return(tmp.k)
}

getbestsensitivity <- function(i, sens){
  bs <- get.sample.measures(i, sens)
  tmp.method   <- apply(bs, 1, which.max)
  tmp.k        <- apply(bs, 1, max, na.rm=TRUE)
  names(tmp.k) <- methods.names[tmp.method]
  return(tmp.k)
}

getbestfmeasure <- function(i, fmeas){
  bfm <- get.sample.measures(i, fmeas)
  tmp.method   <- apply(bfm, 1, which.max)
  tmp.k        <- apply(bfm, 1, max, na.rm=TRUE)
  names(tmp.k) <- methods.names[tmp.method]
  return(tmp.k)
}

getbesterror <- function(i, err){
  be <- get.sample.measures(i, err)
  tmp.method   <- apply(be, 1, which.min)
  tmp.k        <- apply(be, 1, min, na.rm=TRUE)
  names(tmp.k) <- methods.names[tmp.method]
  return(tmp.k)
}

## Script section

# IREF ========================

# load("classification/iref/classification.results.RData")
# 
# kcol  <- tp.complete[1,]
# acc   <- rbind( kcol, accuracy(tp.complete[-1, ], fp.complete[-1, ], tn.complete[-1, ], fn.complete[-1, ]))
# err   <- rbind( kcol, error(acc[-1, ]))
# prec  <- rbind( kcol, precision(tp.complete[-1, ], fp.complete[-1, ]))
# sens  <- rbind( kcol, sensitivity(tp.complete[-1, ], fn.complete[-1, ]))
# spec  <- rbind( kcol, specificity(tn.complete[-1, ], fp.complete[-1, ]))
# fmeas <- fmeasure(prec, sens)
# 
# 
# dres   <- "results/classification/GS/"
# fname  <- paste0(dres, "best_accuracy_corrected.txt")
# 
# for(i in 1:length(tissue.names)){
#   cat(tissue.names[i],"\n", file = fname, append = TRUE)
#   ba <- getbestaccuracy(i)
#   cat(names(ba), "\n", sep = "\t", file = fname, append = TRUE)
#   cat(ba,"\n", sep = " & ", file = fname, append = TRUE)
# }
# 
# 
# fname <- paste0(dres, "best_precision_corrected.txt")
# 
# for(i in 1:length(tissue.names)){
#   cat(tissue.names[i],"\n", file = fname, append = TRUE)
#   ba <- getbestprecision(i)
#   cat(names(ba), "\n", sep = "\t", file = fname, append = TRUE)
#   cat(ba,"\n", sep = " & ", file = fname, append = TRUE)
# }
# 
# 
# fname <- paste0(dres, "best_sensitivity_corrected.txt")
# 
# for(i in 1:length(tissue.names)){
#   cat(tissue.names[i],"\n", file = fname, append = TRUE)
#   ba <- getbestsensitivity(i)
#   cat(names(ba), "\n", sep = "\t", file = fname, append = TRUE)
#   cat(ba,"\n", sep = " & ", file = fname, append = TRUE)
# }
# 
# 
# fname <- paste0(dres, "best_fmeasure_corrected.txt")
# 
# for(i in 1:length(tissue.names)){
#   cat(tissue.names[i],"\n", file = fname, append = TRUE)
#   ba <- getbestfmeasure(i)
#   cat(names(ba), "\n", sep = "\t", file = fname, append = TRUE)
#   cat(ba,"\n", sep = " & ", file = fname, append = TRUE)
# }
# 
# rm(tp.complete, tn.complete, fn.complete, fp.complete)
# rm(kcol, acc, err, prec, sens, spec, fmeas)
# 

# ## INDUCED GS =========================
# 
# load("classification/induced/classification.results_induced.RData")

save.best.results <- function(tp, tn, fp, fn, tissue.names){
  kcol  <- tp[1,]
  acc   <- rbind( kcol, accuracy(tp[-1, ], fp[-1, ], tn[-1, ], fn[-1, ]))
  err   <- rbind( kcol, error(acc[-1, ]))
  prec  <- rbind( kcol, precision(tp[-1, ], fp[-1, ]))
  sens  <- rbind( kcol, sensitivity(tp[-1, ], fn[-1, ]))
  spec  <- rbind( kcol, specificity(tn[-1, ], fp[-1, ]))
  fmeas <- fmeasure(prec, sens)
  
  res.dir <- "results/classification/induced/"
  
  ds.fmeasure <- lapply(1:length(tissue.names), function(i){
    ba <- getbestfmeasure(i, fmeas)
    return(cbind(names(ba), ba, tissue.names[i], 1:4))
  })
  
  ds.precision <- lapply(1:length(tissue.names), function(i){
    ba <- getbestprecision(i, prec)
    return(cbind(names(ba), ba, tissue.names[i], 1:4))
  })
  
  ds.recall <- lapply(1:length(tissue.names), function(i){
    ba <- getbestsensitivity(i, sens)
    return(cbind(names(ba), ba, tissue.names[i], 1:4))
  })
  
  save.results(ds.fmeasure, res.dir, "fmeasure.txt", "fmeasure")
  save.results(ds.precision, res.dir, "precision.txt", "precision")
  save.results(ds.recall, res.dir, "sensitivity.txt", "sensitivity")
  build.best.methods(ds.fmeasure, res.dir, "best.methods_induced.txt", tissue.names)
}

# save.best.results <- function(tp, tn, fp, fn, tissue.names){
#   kcol  <- tp[1,]
#   acc   <- rbind( kcol, accuracy(tp[-1, ], fp[-1, ], tn[-1, ], fn[-1, ]))
#   err   <- rbind( kcol, error(acc[-1, ]))
#   prec  <- rbind( kcol, precision(tp[-1, ], fp[-1, ]))
#   sens  <- rbind( kcol, sensitivity(tp[-1, ], fn[-1, ]))
#   spec  <- rbind( kcol, specificity(tn[-1, ], fp[-1, ]))
#   fmeas <- fmeasure(prec, sens)
# 
#   dres   <- "results/classification/induced/"
#   fname  <- paste0(dres, "best_accuracy_induced.txt")
# 
#   for(i in 1:length(tissue.names)){
#     cat(tissue.names[i],"\n", file = fname, append = TRUE)
#     ba <- getbestaccuracy(i, acc)
#     cat(names(ba), "\n", sep = "\t", file = fname, append = TRUE)
#     cat(ba,"\n", sep = " & ", file = fname, append = TRUE)
#   }
# 
#   fname <- paste0(dres, "best_precision_induced.txt")
# 
#   for(i in 1:length(tissue.names)){
#     cat(tissue.names[i],"\n", file = fname, append = TRUE)
#     ba <- getbestprecision(i, prec)
#     cat(names(ba), "\n", sep = "\t", file = fname, append = TRUE)
#     cat(ba,"\n", sep = " & ", file = fname, append = TRUE)
#   }
# 
# 
#   fname <- paste0(dres, "best_sensitivity_induced.txt")
# 
#   for(i in 1:length(tissue.names)){
#     cat(tissue.names[i],"\n", file = fname, append = TRUE)
#     ba <- getbestsensitivity(i, sens)
#     cat(names(ba), "\n", sep = "\t", file = fname, append = TRUE)
#     cat(ba,"\n", sep = " & ", file = fname, append = TRUE)
#   }
# 
#   fname <- paste0(dres, "best_fmeasure_induced.txt")
# 
#   for(i in 1:length(tissue.names)){
#     cat(tissue.names[i],"\n", file = fname, append = TRUE)
#     ba <- getbestfmeasure(i, fmeas)
#     cat(names(ba), "\n", sep = "\t", file = fname, append = TRUE)
#     cat(ba,"\n", sep = " & ", file = fname, append = TRUE)
#   }
# }
