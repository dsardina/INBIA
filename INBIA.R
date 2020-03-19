options(stringsAsFactors = FALSE)

#################################
##
## LOAD functions and packages 
##
#################################

source("lib/method.functions.R")
source("lib/compute.classification.quality.R")
source("lib/manage.classification.results.R")
source("lib/extract.networks.best.methods.R")
# source("lib/network.statistics.R")
# source("lib/all.common.predictions.R")
source("conf.R")

require(WGCNA)
require(parcor)
require(igraph)
require(ggplot2)
require(parmigene)
require(GeneNet)
require(glasso)

cat(" ### INBIA STARTED ###\t\t\t", paste(Sys.time()),"\n")
cat("============================================================\n")
cat("INBIA current working directory:\n")
cat("\t", getwd(), "\n\n")

#################################
##
## Load GS (IRefIndex)
##
#################################

# Load GS PPI network
human.gs <- read.delim("iref/IRefNet_New_HGNC.txt", header = FALSE)
human.gs.graph <- graph.edge(edges = as.data.frame(human.gs))

#################################
##
## FUNCTION to run all methods
## and save the results.
##
#################################

execute.all <- function(dataset, ds.name){
  View(dataset)
  
  cat(" ### Executing methods for", ds.name, "###\n")
  dataset <- na.omit(dataset) # delete any NA value from dataset
  
  # === PEARSON =========================
  cat("\n..Pearson\n")
  pearson.matrix <- cor(dataset, method = "pearson")
  pearson.edges <- network.test.edges(pearson.matrix, verbose=FALSE, plot=FALSE) # compute significance edges
  
  # === SPEARMAN =========================
  cat("\n..Spearman\n")
  spearman.matrix <- cor(dataset, method = "spearman")
  spearman.edges <- network.test.edges(spearman.matrix, verbose=FALSE, plot=FALSE) # compute significance edges
  
  # === SIMPLEPARCOR =====================
  cat("\n..Simple partial correlation\n")
  spc.matrix <- cor2pcor(cor(dataset, method = "pearson"))
  spc.edges <- network.test.edges(spc.matrix, verbose=FALSE, plot=FALSE)
  
  # === GENENET ==========================
  cat("\n..GeneNet\n")
  genenet.matrix <- ggm.estimate.pcor(dataset)
  genenet.edges <- network.test.edges(genenet.matrix, verbose=FALSE, plot=FALSE)
  
  # === GLASSO ===========================
  cat("\n..Graphical lasso\n")
  sicm <- glasso(cov(dataset), rho=0.01)
  glasso.matrix <- invcov2parcor(sicm$wi)
  glasso.edges <- network.test.edges(glasso.matrix, verbose=FALSE, plot=FALSE)
  
  # === PLSNET ===========================
  cat("\n..Partial least square\n")
  pls.matrix <- pls.net(as.matrix(dataset))
  pls.edges <- network.test.edges(pls.matrix$pcor, verbose=FALSE, plot=FALSE)
    
  # === RIDGE NET ========================
  cat("\n..RidgeNet\n")
  ridge.matrix <- ridge.net(dataset)
  ridge.edges <- network.test.edges(ridge.matrix$pcor, verbose=FALSE, plot=FALSE) # compute significance edges
  
  # === LASSONET =========================
  cat("\n..Lasso\n")
  lasso.matrix <- adalasso.net(dataset, use.Gram = TRUE, both = FALSE)
  lasso.edges <- network.test.edges(lasso.matrix$pcor.lasso, verbose=FALSE, plot=FALSE)
  
  # === ELASTICNET =======================
  cat("\n..ElasticNet\n")
  alphaVec <- seq(0.01, 0.99, length.out=10)
  elasticnet.matrix <- elasticNetwork(dataset, nfold = 3, alpha = alphaVec, verbose=TRUE)
  elasticnet.edges <-  network.test.edges(elasticnet.matrix$pcor, verbose=FALSE, plot=FALSE)

  ## Compute mi matrix from dataset for mutual information methods
  mi <- knnmi.all(t(dataset))
  
  # === ARACNE ===========================
  cat("\n..Aracne\n")
  aracnea.matrix <- aracne.a(mi)
  aracnem.matrix <- aracne.m(mi)
  
  aa.edges <- getEdges(aracnea.matrix)
  am.edges <- getEdges(aracnem.matrix)
  
  # === CLR ==============================
  cat("\n..CLR\n")
  clr.matrix <- clr(mi)
  clr.edges <- getEdges(clr.matrix)
  
  # === MRNET ============================
  cat("\n..MRNet\n")
  mrnet.matrix <- mrnet(mi)
  mrnet.edges <- getEdges(mrnet.matrix)
  
  # === WGCNA ============================
  cat("\n..WGCNA\n")
  wgcna.matrix <- wgcna(dataset, ds.name)
  wgcna.edges  <- network.test.edges(wgcna.matrix, verbose=FALSE, plot=FALSE)
  
  cat(" ### Saving results ###\n\n")
  return(list(dataset=dataset,
              pearson.matrix=pearson.matrix,       pearson.edges=pearson.edges,
              spearman.matrix=spearman.matrix,     spearman.edges=spearman.edges,
              spc.edges=spc.edges,                 spc.matrix=spc.matrix,
              genenet.matrix=genenet.matrix,       genenet.edges=genenet.edges,
              glasso.matrix=glasso.matrix,         glasso.edges=glasso.edges,
              pls.matrix=pls.matrix,               pls.edges=pls.edges,
              ridge.matrix=ridge.matrix,           ridge.edges=ridge.edges,
              lasso.matrix=lasso.matrix,           lasso.edges=lasso.edges,
              elasticnet.matrix=elasticnet.matrix, elasticnet.edges=elasticnet.edges,
              aracnea.matrix=aracnea.matrix,       aa.edges=aa.edges,
              aracnem.matrix=aracnem.matrix,       am.edges=am.edges,
              clr.matrix=clr.matrix,               clr.edges=clr.edges,
              mrnet.matrix=mrnet.matrix,           mrnet.edges=mrnet.edges,
              wgcna.matrix=wgcna.matrix,           wgcna.edges=wgcna.edges))
}

##################################################
##
## Function for running all methods and collect
## the results.
##
##################################################

# ress <- lapply(multi.ds, execute.all)
# save(ress, file = "rdata/ress.RData")

##################################################
##
## Comparison between methods predictions and GS.
##
##################################################

get.methods.intersections <- function(ress, convTab, induced.human.gs.converted){
  lapply(ress, function(tissue){
    cat(" ### Computing common PPIs ###\n")
    
    tissue.pearson    <- tissue$pearson.edges
    tissue.spearman   <- tissue$spearman.edges
    tissue.spc        <- tissue$spc.edges
    tissue.genenet    <- tissue$genenet.edges
    tissue.glasso     <- tissue$glasso.edges
    tissue.pls        <- tissue$pls.edges
    tissue.ridgenet   <- tissue$ridge.edges
    tissue.lasso      <- tissue$lasso.edges
    tissue.elasticnet <- tissue$elasticnet.edges
    tissue.aracnea    <- tissue$aa.edges
    tissue.aracnem    <- tissue$am.edges
    tissue.clr        <- tissue$clr.edges
    tissue.mrnet      <- tissue$mrnet.edges
    tissue.wgcna      <- tissue$wgcna.edges
    
    cat("..Computing significant edges\n")
    
    # Get significant edges with symbols from prediction methods
    final.pearson    <- significant.edges(tissue.pearson, colnames(tissue$data))
    final.spearman   <- significant.edges(tissue.spearman, colnames(tissue$data))
    final.spc        <- significant.edges(tissue.spc, colnames(tissue$data))
    final.genenet    <- significant.edges(tissue.genenet, colnames(tissue$data))
    final.glasso     <- significant.edges(tissue.glasso, colnames(tissue$data))
    final.pls        <- significant.edges(tissue.pls, colnames(tissue$data))
    final.ridgenet   <- significant.edges(tissue.ridgenet, colnames(tissue$data))
    final.lasso      <- significant.edges(tissue.lasso, colnames(tissue$data))
    final.elasticnet <- significant.edges(tissue.elasticnet, colnames(tissue$data))
    final.aracnea    <- net.from.index(tissue.aracnea, colnames(tissue$data))
    final.aracnem    <- net.from.index(tissue.aracnem, colnames(tissue$data))
    final.clr        <- net.from.index(tissue.clr, colnames(tissue$data))
    final.mrnet      <- net.from.index(tissue.mrnet, colnames(tissue$data))
    final.wgcna      <- significant.edges(tissue.wgcna, colnames(tissue$data))
    
    cat("..Converting edges\n")
    
    # Annotate phosphoproteins
    annotated.pearson     <- get.phospho.converted.interactions(final.pearson[,-1], convTab)
    annotated.spearman    <- get.phospho.converted.interactions(final.spearman[,-1], convTab)
    annotated.spc         <- get.phospho.converted.interactions(final.spc[,-1], convTab)
    annotated.genenet     <- get.phospho.converted.interactions(final.genenet[,-1], convTab)
    annotated.glasso      <- get.phospho.converted.interactions(final.glasso[,-1], convTab)
    annotated.pls         <- get.phospho.converted.interactions(final.pls[,-1], convTab)
    annotated.ridgenet    <- get.phospho.converted.interactions(final.ridgenet[,-1], convTab)
    annotated.lasso       <- get.phospho.converted.interactions(final.lasso[,-1], convTab)
    annotated.elasticnet  <- get.phospho.converted.interactions(final.elasticnet[,-1], convTab)  
    annotated.aracnea     <- get.phospho.converted.interactions(final.aracnea, convTab)
    annotated.aracnem     <- get.phospho.converted.interactions(final.aracnem, convTab)
    annotated.clr         <- get.phospho.converted.interactions(final.clr, convTab)
    annotated.mrnet       <- get.phospho.converted.interactions(final.mrnet, convTab)
    annotated.wgcna       <- get.phospho.converted.interactions(final.wgcna[,-1], convTab)
    
    cat("..Computing intersection with the GS\n")
    
    # Create tissue graph from current method predictions
    pearson.intersection    <- edges.intersection(annotated.pearson, induced.human.gs.converted)
    spearman.intersection   <- edges.intersection(annotated.spearman, induced.human.gs.converted)
    spc.intersection        <- edges.intersection(annotated.spc, induced.human.gs.converted)
    genenet.intersection    <- edges.intersection(annotated.genenet, induced.human.gs.converted)
    glasso.intersection     <- edges.intersection(annotated.glasso, induced.human.gs.converted)
    pls.intersection        <- edges.intersection(annotated.pls, induced.human.gs.converted)
    ridgenet.intersection   <- edges.intersection(annotated.ridgenet, induced.human.gs.converted)
    lasso.intersection      <- edges.intersection(annotated.lasso, induced.human.gs.converted)
    elasticnet.intersection <- edges.intersection(annotated.elasticnet, induced.human.gs.converted)
    aracnea.intersection    <- edges.intersection(annotated.aracnea, induced.human.gs.converted)
    aracnem.intersection    <- edges.intersection(annotated.aracnem, induced.human.gs.converted)
    clr.intersection        <- edges.intersection(annotated.clr, induced.human.gs.converted)
    mrnet.intersection      <- edges.intersection(annotated.mrnet, induced.human.gs.converted)
    wgcna.intersection      <- edges.intersection(annotated.wgcna, induced.human.gs.converted)
    
    cat("..Saving results\n\n")
    return(list(final.pearson = final.pearson, annotated.pearson = annotated.pearson, pearson.intersection = pearson.intersection,
                final.spearman = final.spearman, annotated.spearman = annotated.spearman, spearman.intersection = spearman.intersection,
                final.spc = final.spc, annotated.spc = annotated.spc, spc.intersection = spc.intersection,
                final.genenet = final.genenet, annotated.genenet = annotated.genenet, genenet.intersection = genenet.intersection,
                final.glasso = final.glasso, annotated.glasso = annotated.glasso, glasso.intersection = glasso.intersection,
                final.pls = final.pls, annotated.pls = annotated.pls, pls.intersection = pls.intersection,
                final.ridgenet = final.ridgenet, annotated.ridgenet = annotated.ridgenet, ridgenet.intersection = ridgenet.intersection,
                final.lasso = final.lasso, annotated.lasso = annotated.lasso, lasso.intersection = lasso.intersection,
                final.elasticnet = final.elasticnet, annotated.elasticnet = annotated.elasticnet, elasticnet.intersection = elasticnet.intersection,
                final.aracnea = final.aracnea, annotated.aracnea = annotated.aracnea, aracnea.intersection = aracnea.intersection,
                final.aracnem = final.aracnem, annotated.aracnem = annotated.aracnem, aracnem.intersection = aracnem.intersection,
                final.clr = final.clr, annotated.clr = annotated.clr, clr.intersection = clr.intersection,
                final.mrnet = final.mrnet, annotated.mrnet = annotated.mrnet, mrnet.intersection = mrnet.intersection,
                final.wgcna = final.wgcna, annotated.wgcna = annotated.wgcna, wgcna.intersection = wgcna.intersection))
  })
}

# save(methods.intersections, file = "rdata/methods.intersection.results_all_new.RData")


##################################################
##
## Compute shortest path for methods' predictions
## and perform statistical analysis of the
## comparisons with GS.
##
##################################################

## Init a list for predictions results
get.methods.predictions <- function(methods.intersections,
                                    tissue.names=NULL,
                                    induced.human.gs.converted,
                                    save.plots=TRUE,
                                    save.stats=FALSE){
  methods.predictions <- vector("list", length(methods.intersections))
  
  if(length(methods.intersections) != length(tissue.names)){
    warning("Tissue names are not the same length of methods intersections. Use default names.\n")
    tissue.names <- paste0("tissue_", 1:length(methods.intersections))
  }
  
  ## Load GS genes and graph
  gs.genes <- unique.genes(induced.human.gs.converted)
  gs.graph <- graph.edge(edges = as.data.frame(induced.human.gs.converted))
  
  for(i in 1:length(methods.predictions)){
    tissue <- tissue.names[i]
    dir.create(paste0("results/", tissue), showWarnings = FALSE)
    
    cat(" ### Predictions results for", tissue, "###\n")
    
    # Load predictions
    elem <- methods.intersections[[i]]
    
    pearson.ds    <- elem$annotated.pearson
    spearman.ds   <- elem$annotated.spearman
    spc.ds        <- elem$annotated.spc
    genenet.ds    <- elem$annotated.genenet
    glasso.ds     <- elem$annotated.glasso
    pls.ds        <- elem$annotated.pls
    ridge.ds      <- elem$annotated.ridgenet
    lasso.ds      <- elem$annotated.lasso
    elasticnet.ds <- elem$annotated.elasticnet
    aracnea.ds    <- elem$annotated.aracnea
    aracnem.ds    <- elem$annotated.aracnem
    clr.ds        <- elem$annotated.clr
    mrnet.ds      <- elem$annotated.mrnet
    wgcna.ds      <- elem$annotated.wgcna
    
    cat("..Comparing methods predictions\n")
    
    # Set if they are present in GS
    pearson.ds$present.gs    <- elem$pearson.intersection
    spearman.ds$present.gs   <- elem$spearman.intersection
    spc.ds$present.gs        <- elem$spc.intersection
    genenet.ds$present.gs    <- elem$genenet.intersection
    glasso.ds$present.gs     <- elem$glasso.intersection
    pls.ds$present.gs        <- elem$pls.intersection
    ridge.ds$present.gs      <- elem$ridgenet.intersection
    lasso.ds$present.gs      <- elem$lasso.intersection
    elasticnet.ds$present.gs <- elem$elasticnet.intersection
    aracnea.ds$present.gs    <- elem$aracnea.intersection
    aracnem.ds$present.gs    <- elem$aracnem.intersection
    clr.ds$present.gs        <- elem$clr.intersection
    mrnet.ds$present.gs      <- elem$mrnet.intersection
    wgcna.ds$present.gs      <- elem$wgcna.intersection
    
    cat("..Filtering loop and predictions not in GS\n")
    
    # Filter all the predictions if proteins are not in GS or the interaction is a loop
    pearson.ds    <- pearson.ds[filter.predictions(pearson.ds, gs.genes), ]
    spearman.ds   <- spearman.ds[filter.predictions(spearman.ds, gs.genes), ]
    spc.ds        <- spc.ds[filter.predictions(spc.ds, gs.genes), ]
    genenet.ds    <- genenet.ds[filter.predictions(genenet.ds, gs.genes), ]
    glasso.ds     <- glasso.ds[filter.predictions(glasso.ds, gs.genes), ]
    pls.ds        <- pls.ds[filter.predictions(pls.ds, gs.genes), ]
    ridge.ds      <- ridge.ds[filter.predictions(ridge.ds, gs.genes), ]
    lasso.ds      <- lasso.ds[filter.predictions(lasso.ds, gs.genes), ]
    elasticnet.ds <- elasticnet.ds[filter.predictions(elasticnet.ds, gs.genes), ]
    aracnea.ds    <- aracnea.ds[filter.predictions(aracnea.ds, gs.genes), ]
    aracnem.ds    <- aracnem.ds[filter.predictions(aracnem.ds, gs.genes), ]
    clr.ds        <- clr.ds[filter.predictions(clr.ds, gs.genes), ]
    mrnet.ds      <- mrnet.ds[filter.predictions(mrnet.ds, gs.genes), ]
    wgcna.ds      <- wgcna.ds[filter.predictions(wgcna.ds, gs.genes), ]
    
    cat("..Computing shortest paths")
    
    # Compute shortest paths with induced GS
    pearson.ds$gs.shortest.path    <- compute.shortest.paths(pearson.ds, gs.graph)
    cat(" ...", round(50/14*1), "%")
    spearman.ds$gs.shortest.path   <- compute.shortest.paths(spearman.ds, gs.graph)
    cat(" ...", round(50/14*2), "%")
    spc.ds$gs.shortest.path        <- compute.shortest.paths(spc.ds, gs.graph)
    cat(" ...", round(50/14*3), "%")
    genenet.ds$gs.shortest.path    <- compute.shortest.paths(genenet.ds, gs.graph)
    cat(" ...", round(50/14*4), "%")
    glasso.ds$gs.shortest.path     <- compute.shortest.paths(glasso.ds, gs.graph)
    cat(" ...", round(50/14*5), "%")
    pls.ds$gs.shortest.path        <- compute.shortest.paths(pls.ds, gs.graph)
    cat(" ...", round(50/14*6), "%")
    ridge.ds$gs.shortest.path      <- compute.shortest.paths(ridge.ds, gs.graph)
    cat(" ...", round(50/14*7), "%")
    lasso.ds$gs.shortest.path      <- compute.shortest.paths(lasso.ds, gs.graph)
    cat(" ...", round(50/14*8), "%")
    elasticnet.ds$gs.shortest.path <- compute.shortest.paths(elasticnet.ds, gs.graph)
    cat(" ...", round(50/14*9), "%")
    aracnea.ds$gs.shortest.path    <- compute.shortest.paths(aracnea.ds, gs.graph)
    cat(" ...", round(50/14*10), "%")
    aracnem.ds$gs.shortest.path    <- compute.shortest.paths(aracnem.ds, gs.graph)
    cat(" ...", round(50/14*11), "%")
    clr.ds$gs.shortest.path        <- compute.shortest.paths(clr.ds, gs.graph)
    cat(" ...", round(50/14*12), "%")
    mrnet.ds$gs.shortest.path      <- compute.shortest.paths(mrnet.ds, gs.graph)
    cat(" ...", round(50/14*13), "%")
    wgcna.ds$gs.shortest.path      <- compute.shortest.paths(wgcna.ds, gs.graph)
    cat(" ...", 50, "%")
    
    # Compute shortest paths with human GS 
    pearson.ds$human.gs.shortest.path    <- compute.shortest.paths(pearson.ds, human.gs.graph)
    cat(" ...", round(50+50/14*1), "%")
    spearman.ds$human.gs.shortest.path   <- compute.shortest.paths(spearman.ds, human.gs.graph)
    cat(" ...", round(50+50/14*2), "%")
    spc.ds$human.gs.shortest.path        <- compute.shortest.paths(spc.ds, human.gs.graph)
    cat(" ...", round(50+50/14*3), "%")
    genenet.ds$human.gs.shortest.path    <- compute.shortest.paths(genenet.ds, human.gs.graph)
    cat(" ...", round(50+50/14*4), "%")
    glasso.ds$human.gs.shortest.path     <- compute.shortest.paths(glasso.ds, human.gs.graph)
    cat(" ...", round(50+50/14*5), "%")
    pls.ds$human.gs.shortest.path        <- compute.shortest.paths(pls.ds, human.gs.graph)
    cat(" ...", round(50+50/14*6), "%")
    ridge.ds$human.gs.shortest.path      <- compute.shortest.paths(ridge.ds, human.gs.graph)
    cat(" ...", round(50+50/14*7), "%")
    lasso.ds$human.gs.shortest.path      <- compute.shortest.paths(lasso.ds, human.gs.graph)
    cat(" ...", round(50+50/14*8), "%")
    elasticnet.ds$human.gs.shortest.path <- compute.shortest.paths(elasticnet.ds, human.gs.graph)
    cat(" ...", round(50+50/14*9), "%")
    aracnea.ds$human.gs.shortest.path    <- compute.shortest.paths(aracnea.ds, human.gs.graph)
    cat(" ...", round(50+50/14*10), "%")
    aracnem.ds$human.gs.shortest.path    <- compute.shortest.paths(aracnem.ds, human.gs.graph)
    cat(" ...", round(50+50/14*11), "%")
    clr.ds$human.gs.shortest.path        <- compute.shortest.paths(clr.ds, human.gs.graph)
    cat(" ...", round(50+50/14*12), "%")
    mrnet.ds$human.gs.shortest.path      <- compute.shortest.paths(mrnet.ds, human.gs.graph)
    cat(" ...", round(50+50/14*13), "%")
    wgcna.ds$human.gs.shortest.path      <- compute.shortest.paths(wgcna.ds, human.gs.graph)
    cat(" ...", 100, "%. End.\n")
    
    cat("..Computing path lengths\n")
    
    # Plot path length distribution with induced GS
    pearson.length.distr    <- get.path.length.distr(pearson.ds$gs.shortest.path)
    spearman.length.distr   <- get.path.length.distr(spearman.ds$gs.shortest.path)
    spc.length.distr        <- get.path.length.distr(spc.ds$gs.shortest.path)
    genenet.length.distr    <- get.path.length.distr(genenet.ds$gs.shortest.path)
    glasso.length.distr     <- get.path.length.distr(glasso.ds$gs.shortest.path)
    pls.length.distr        <- get.path.length.distr(pls.ds$gs.shortest.path)
    ridge.length.distr      <- get.path.length.distr(ridge.ds$gs.shortest.path)
    lasso.length.distr      <- get.path.length.distr(lasso.ds$gs.shortest.path)
    elasticnet.length.distr <- get.path.length.distr(elasticnet.ds$gs.shortest.path)
    aracnea.length.distr    <- get.path.length.distr(aracnea.ds$gs.shortest.path)
    aracnem.length.distr    <- get.path.length.distr(aracnem.ds$gs.shortest.path)
    clr.length.distr        <- get.path.length.distr(clr.ds$gs.shortest.path)
    mrnet.length.distr      <- get.path.length.distr(mrnet.ds$gs.shortest.path)
    wgcna.length.distr      <- get.path.length.distr(wgcna.ds$gs.shortest.path)
    
    pearson.ds$path.length    <- pearson.length.distr
    spearman.ds$path.length   <- spearman.length.distr
    spc.ds$path.length        <- spc.length.distr
    genenet.ds$path.length    <- genenet.length.distr
    glasso.ds$path.length     <- glasso.length.distr
    pls.ds$path.length        <- pls.length.distr
    ridge.ds$path.length      <- ridge.length.distr
    lasso.ds$path.length      <- lasso.length.distr
    elasticnet.ds$path.length <- elasticnet.length.distr
    aracnea.ds$path.length    <- aracnea.length.distr
    aracnem.ds$path.length    <- aracnem.length.distr
    clr.ds$path.length        <- clr.length.distr
    mrnet.ds$path.length      <- mrnet.length.distr
    wgcna.ds$path.length      <- wgcna.length.distr
    
    
    # Plot path length distribution with induced GS
    pearson.length.distr.human    <- get.path.length.distr(pearson.ds$human.gs.shortest.path)
    spearman.length.distr.human   <- get.path.length.distr(spearman.ds$human.gs.shortest.path)
    spc.length.distr.human        <- get.path.length.distr(spc.ds$human.gs.shortest.path)
    genenet.length.distr.human    <- get.path.length.distr(genenet.ds$human.gs.shortest.path)
    glasso.length.distr.human     <- get.path.length.distr(glasso.ds$human.gs.shortest.path)
    pls.length.distr.human        <- get.path.length.distr(pls.ds$human.gs.shortest.path)
    ridge.length.distr.human      <- get.path.length.distr(ridge.ds$human.gs.shortest.path)
    lasso.length.distr.human      <- get.path.length.distr(lasso.ds$human.gs.shortest.path)
    elasticnet.length.distr.human <- get.path.length.distr(elasticnet.ds$human.gs.shortest.path)
    aracnea.length.distr.human    <- get.path.length.distr(aracnea.ds$human.gs.shortest.path)
    aracnem.length.distr.human    <- get.path.length.distr(aracnem.ds$human.gs.shortest.path)
    clr.length.distr.human        <- get.path.length.distr(clr.ds$human.gs.shortest.path)
    mrnet.length.distr.human      <- get.path.length.distr(mrnet.ds$human.gs.shortest.path)
    wgcna.length.distr.human      <- get.path.length.distr(wgcna.ds$human.gs.shortest.path)
    
    pearson.ds$path.length.human    <- pearson.length.distr.human
    spearman.ds$path.length.human   <- spearman.length.distr.human
    spc.ds$path.length.human        <- spc.length.distr.human
    genenet.ds$path.length.human    <- genenet.length.distr.human
    glasso.ds$path.length.human     <- glasso.length.distr.human
    pls.ds$path.length.human        <- pls.length.distr.human
    ridge.ds$path.length.human      <- ridge.length.distr.human
    lasso.ds$path.length.human      <- lasso.length.distr.human
    elasticnet.ds$path.length.human <- elasticnet.length.distr.human
    aracnea.ds$path.length.human    <- aracnea.length.distr.human
    aracnem.ds$path.length.human    <- aracnem.length.distr.human
    clr.ds$path.length.human        <- clr.length.distr.human
    mrnet.ds$path.length.human      <- mrnet.length.distr.human
    wgcna.ds$path.length.human      <- wgcna.length.distr.human
    
    cat("..Saving predictions results\n")
    
    methods.predictions[[i]] <- list(pearson.ds=pearson.ds,
                                     spearman.ds=spearman.ds,
                                     spc.ds=spc.ds,
                                     genenet.ds=genenet.ds,
                                     glasso.ds=glasso.ds,
                                     pls.ds=pls.ds,
                                     ridge.ds=ridge.ds,
                                     lasso.ds=lasso.ds,
                                     elasticnet.ds=elasticnet.ds,
                                     aracnea.ds=aracnea.ds,
                                     aracnem.ds=aracnem.ds,
                                     clr.ds=clr.ds,
                                     mrnet.ds=mrnet.ds,
                                     wgcna.ds=wgcna.ds )
    
    if(save.plots){
      complete.ds <- do.call(rbind, methods.predictions[[i]])
      
      labels <- lapply(1:length(methods.predictions[[i]]), function(mid){
        rep(methods.names[mid], nrow(methods.predictions[[i]][[mid]]))
      })
      
      length.ds       <- data.frame(length=complete.ds$path.length, method=unlist(labels))
      length.ds.human <- data.frame(length=complete.ds$path.length.human, method=unlist(labels))
      
      cat("..Create and save plots\n")
      
      ggsave(paste0("results/", tissue, "/methods.path.length_stacked.jpg"), methods.sp.plot(length.ds))
      ggsave(paste0("results/", tissue, "/methods.path.length.jpg"), methods.sp.plot(length.ds, stacked = FALSE))
      
      ggsave(paste0("results/", tissue, "/methods.path.length.human_stacked.jpg"), methods.sp.plot(length.ds.human))
      ggsave(paste0("results/", tissue, "/methods.path.length.human.jpg"), methods.sp.plot(length.ds.human, stacked = FALSE))
    }
    
    if(save.stats){
      cat("..Computing statistics\n")
      
      # statistics about protein interactions and phosphorilation
      both <- c(sum(both.phospho(pearson.ds)),
                sum(both.phospho(spearman.ds)),
                sum(both.phospho(spc.ds)),
                sum(both.phospho(genenet.ds)),
                sum(both.phospho(glasso.ds)),
                sum(both.phospho(pls.ds)),
                sum(both.phospho(ridge.ds)),
                sum(both.phospho(lasso.ds)),
                sum(both.phospho(elasticnet.ds)),
                sum(both.phospho(aracnea.ds)),
                sum(both.phospho(aracnem.ds)),
                sum(both.phospho(clr.ds)),
                sum(both.phospho(mrnet.ds)),
                sum(both.phospho(wgcna.ds)) )
      
      atleastone <- c(sum(atleastone.phospho(pearson.ds)),
                      sum(atleastone.phospho(spearman.ds)),
                      sum(atleastone.phospho(spc.ds)),
                      sum(atleastone.phospho(genenet.ds)),
                      sum(atleastone.phospho(glasso.ds)),
                      sum(atleastone.phospho(pls.ds)),
                      sum(atleastone.phospho(ridge.ds)),
                      sum(atleastone.phospho(lasso.ds)),
                      sum(atleastone.phospho(elasticnet.ds)),
                      sum(atleastone.phospho(aracnea.ds)),
                      sum(atleastone.phospho(aracnem.ds)),
                      sum(atleastone.phospho(clr.ds)),
                      sum(atleastone.phospho(mrnet.ds)),
                      sum(atleastone.phospho(wgcna.ds)) )
      
      none <- c(sum(none.phospho(pearson.ds)),
                sum(none.phospho(spearman.ds)),
                sum(none.phospho(spc.ds)),
                sum(none.phospho(genenet.ds)),
                sum(none.phospho(glasso.ds)),
                sum(none.phospho(pls.ds)),
                sum(none.phospho(ridge.ds)),
                sum(none.phospho(lasso.ds)),
                sum(none.phospho(elasticnet.ds)),
                sum(none.phospho(aracnea.ds)),
                sum(none.phospho(aracnem.ds)),
                sum(none.phospho(clr.ds)),
                sum(none.phospho(mrnet.ds)),
                sum(none.phospho(wgcna.ds)) )
      
      na.amount <- c(sum(is.na(pearson.ds$gs.shortest.path)),
                     sum(is.na(spearman.ds$gs.shortest.path)),
                     sum(is.na(spc.ds$gs.shortest.path)),
                     sum(is.na(genenet.ds$gs.shortest.path)),
                     sum(is.na(glasso.ds$gs.shortest.path)),
                     sum(is.na(pls.ds$gs.shortest.path)),
                     sum(is.na(ridge.ds$gs.shortest.path)),
                     sum(is.na(lasso.ds$gs.shortest.path)),
                     sum(is.na(elasticnet.ds$gs.shortest.path)),
                     sum(is.na(aracnea.ds$gs.shortest.path)),
                     sum(is.na(aracnem.ds$gs.shortest.path)),
                     sum(is.na(clr.ds$gs.shortest.path)),
                     sum(is.na(mrnet.ds$gs.shortest.path)),
                     sum(is.na(wgcna.ds$gs.shortest.path)) )
      
      stat <- as.data.frame(rbind(both, atleastone, none, na.amount))
      colnames(stat) <- methods.names
      
      cat("..Saving stats\n\n")
      
      write.table(stat, paste0("results/", tissue, "/phosho.interactions.statistics.txt"), quote = FALSE, sep = "\t")
    }
  }
  return(methods.predictions)
}

# save(methods.predictions, file = "rdata/methods.predictions_all.RData")
