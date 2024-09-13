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
results.dir = "~/projects/demux-benchmarking/manuscript/results/revision//run_time_simulated/seurat_objects/"

if(!dir.exists(results.dir)) {
  dir.create(results.dir, recursive = T)
}

n.bcs = c(5,10,20,40)
max.n.cells = c(100, 200,400,800)

doublet.rate= 0.1 
min.size.log = log(10)
max.size.log = log(1000)
min.ambient.log = 0

max.ambient = 100
b0 = -4

seed = 1

for(n.bc in n.bcs){
  for(max.n.cell in max.n.cells){
      max.ambient.log = log(max.ambient)
      min.n.cell = max.n.cell
      
      results = simulateTags2(min.n.cell = min.n.cell , max.n.cell = max.n.cell,n.bc = n.bc,min.size.log = min.size.log,
                              max.size.log = max.size.log,min.ambient.log = min.ambient.log, 
                              max.ambient.log = max.ambient.log,
                              doublet.rate = doublet.rate, b0 = b0,seed = seed)
      
      tag_mtx = results$counts
      counts = tag_mtx %>% t()
      seurat_object = CreateSeuratObject(counts = counts,assay = "HTO")
      
      labels = stringr::str_split(string = colnames(counts), pattern = "_")
      labels = lapply(labels,FUN = function(x) x[1]) %>% unlist()
      seurat_object$ground_truth = as.character(labels)
      
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
      
      saveRDS(seurat_object, file = sprintf("%s/%s_%g.%s_%g.rds",results.dir,"ncells",max.n.cell,"n.bc", n.bc))   
  }
}

