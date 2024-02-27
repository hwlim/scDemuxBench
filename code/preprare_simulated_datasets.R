library(Seurat)
library(dplyr)
source("functions.R")

####################################################################
# This script generates simulated datasets using a modified version of simulateTags function 
# from deMULTIplex2 package (Zhu et. at., 2024)
# Modified from https://github.com/Gartner-Lab/deMULTIplex2/blob/92130a626249194ef7f6c0b1aecad6ac210e258b/R/simulation.R to 
# 1) allow for different number of cells per tag
# 2) guarantee that one barcode's ambient contamination level is eual to `max.ambient.log`, 
# and other barcodes' ambient contamination levels are uniformly sampled from range [min.ambient.log , min.ambient.log]
####################################################################

# folder to save downloaded raw data
results.dir = "../results//simulated/"

if(!dir.exists(results.dir)) {
  dir.create(results.dir, recursive = T)
}

n.bc = 5
min.n.cell = 1000
max.n.cell = 1000
doublet.rate= 0.1 
min.size.log = log(10)
max.size.log = log(1000)
min.ambient.log = 0

max.ambients = seq(100,2000,200)
cell_bound = c(-4,-3,-2,-1)

seeds = c(1,2,3)

for(max.ambient in max.ambients){
  max.ambient.log = log(max.ambient)
  for(b0 in cell_bound){
    for(seed in seeds){
      results = simulateTags2(min.n.cell = min.n.cell , max.n.cell = max.n.cell,n.bc = n.bc,min.size.log = min.size.log,
                              max.size.log = max.size.log,min.ambient.log = min.ambient.log, 
                              max.ambient.log = max.ambient.log,
                              doublet.rate = doublet.rate, b0 = b0,seed = seed)
      
      tag_mtx = results$counts
      log(summary(rowSums(tag_mtx)))
      seurat_object = CreateSeuratObject(counts = tag_mtx %>% t(),assay = "HTO")
      seurat_object$ground_truth = as.character(seurat_object$orig.ident)
      
      seurat_object$ground_truth[seurat_object$ground_truth == "doublet"] = "Doublet"
      
      seurat_object@misc[["simulation_params"]] = list(min.n.cell = min.n.cell,
                                                       max.n.cell =max.n.cell,
                                                       num_tags = n.bc,
                                                       doublet.rate = doublet.rate,
                                                       min.size.log = min.size.log,
                                                       max.size.log = max.size.log,
                                                       min.ambient.log = min.ambient.log,
                                                       max.ambient.log = max.ambient.log,
                                                       b0=b0,
                                                       seed = seed)
      
      saveRDS(seurat_object, file = sprintf("%s/%s_%g.%s_%g.%s_%g.rds",results.dir,"max.ambient",max.ambient,"b0", -b0,"rep",seed))   
    }
  }
}
