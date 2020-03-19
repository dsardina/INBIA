setwd("D:/Dottorato/material/magi")
options(stringsAsFactors = FALSE)

library(igraph)
source("scripts/method.functions.R")
source("scripts/tissue.methods.names.R")

load("bestanalysis/tissuenet_k1.RData")
iref.tissue.nets <- tissue.nets
rm(tissue.nets)

load("pera/bestmethods/tissuenet_k1.RData")
pera.tissue.nets <- tissue.nets
rm(tissue.nets)

iref.net.measures <- lapply(iref.tissue.nets, function(tnets){
  g <- tnets$tissue.kgraph
  
  edensity <- edge_density(g)
  transitiv <- transitivity(g, type="global")
  diam <- diameter(g, directed = FALSE)
  ddistr <- degree.distribution(g, cumulative = TRUE, mode = "all")
  cdegree <- centr_degree(g, mode="all", normalized=T)
  closen <- centr_clo(g, mode="all", normalized=T)
  betwewn <- centr_betw(g, directed=F, normalized=T)
  hub.score <- hub_score(g, weights=NA)$vector
  ascore <- authority_score(g, weights=NA)$vector
  mdist <- mean_distance(g, directed=F)
  all.cliques <- cliques(g)
  lcliques <- largest_cliques(g)
  return(list(g=g, edensity=edensity, transitiv=transitiv, diam=diam, ddistr=ddistr,
              cdegree=cdegree, closen=closen, betwewn=betwewn, hub.score=hub.score,
              ascore=ascore, mdist=mdist, all.cliques=all.cliques, lcliques=lcliques))
})
names(iref.net.measures) <- tissue.names

pera.net.measures <- lapply(pera.tissue.nets, function(tnets){
  g <- tnets$tissue.kgraph
  edensity <- edge_density(g)
  transitiv <- transitivity(g, type="global")
  diam <- diameter(g, directed = FALSE)
  ddistr <- degree.distribution(g, cumulative = TRUE, mode = "all")
  cdegree <- centr_degree(g, mode="all", normalized=T)
  closen <- centr_clo(g, mode="all", normalized=T)
  betwewn <- centr_betw(g, directed=F, normalized=T)
  hub.score <- hub_score(g, weights=NA)$vector
  ascore <- authority_score(g, weights=NA)$vector
  mdist <- mean_distance(g, directed=F)
  all.cliques <- cliques(g)
  lcliques <- largest_cliques(g)
  return(list(g=g, edensity=edensity, transitiv=transitiv, diam=diam, ddistr=ddistr,
              cdegree=cdegree, closen=closen, betwewn=betwewn, hub.score=hub.score,
              ascore=ascore, mdist=mdist, all.cliques=all.cliques, lcliques=lcliques))
})
names(pera.net.measures) <- tissue.names

# save(iref.net.measures, pera.net.measures, file = "rdata/iref.pera.netmeasuers.RData")


print.statistics <- function(ds, type){
  
  png(paste0("mdsanalysis/", type, "_tissue.degree.distribution.png"), width = 8, height = 8, unit = "in", res = 600)
  par(mfrow=c(4,4))
  for(i in 1:length(tissue.names)){
    tissue <- tissue.names[i]
    cat(tissue, "\n")
    cat(tissue, "\n", file = "mdsanalysis/bnet.stats.txt", append = TRUE)
    tstat <- ds[[tissue]]
    cat(length(V(tstat$g)), "\n")
    cat(length(E(tstat$g)), "\n")
    cat("Edge density: ", tstat$edensity, "\n", file = "mdsanalysis/bnet.stats.txt", append = TRUE)
    cat("Transitivity:", tstat$transitiv, "\n", file = "mdsanalysis/bnet.stats.txt", append = TRUE)
    cat("Diameter:", tstat$diam, "\n", file = "mdsanalysis/bnet.stats.txt", append = TRUE)
    
    plot(tstat$ddistr, type = "l", main = paste(tissue, "degree distribution"), xlab = "nodes", ylab = "cumulative degree")
    
    # cat("Degree centrality:", tstat$cdegree$centralization, "\n", file = "mdsanalysis/bnet.stats.txt", append = TRUE)
    # cat("Closeness centrality:", tstat$closen$centralization, "\n", file = "mdsanalysis/bnet.stats.txt", append = TRUE)
    # cat("Betweenness centrality:", tstat$betwewn$centralization, "\n", file = "mdsanalysis/bnet.stats.txt", append = TRUE)
    # tstat$hub.score
    # tstat$ascore
    cat("Mean distance:", tstat$mdist, "\n", file = "mdsanalysis/bnet.stats.txt", append = TRUE)
    
    # png()
    # plot(tstat$all.cliques, type = "l", main = paste0(tissue, ""), xlab = "nodes", ylab = "clique size distribution")
    # dev.off()
    
    # tstat$lcliques
  }
  
  dev.off()
}

print.statistics(iref.net.measures, "iref")
print.statistics(pera.net.measures, "pera")

png(paste0("mdsanalysis/tissue.degree.distribution.png"), width = 8, height = 8, unit = "in", res = 600)
par(mfrow=c(4,4))
for(i in 1:length(tissue.names)){
  tissue <- tissue.names[i]
  iref.tstat <- iref.net.measures[[tissue]]
  pera.tstat <- pera.net.measures[[tissue]]
  
  plot(iref.tstat$ddistr, type = "l", col = "blue", main = paste(tissue, "degree distribution"), xlab = "nodes", ylab = "cumulative degree")
  lines(pera.tstat$ddistr, col = "orange")
}
dev.off()
