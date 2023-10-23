#!/usr/bin/env Rscript

library(deMULTIplex2)
library(doParallel)
library(cellhashR)
library(demuxmix)
library(Seurat)
library(dplyr)

source("./sample_demultiplexing.R")
source("./benchmarking.R")

#source("./functions.R")

set.seed(1)

###### load Seurat objects
args = commandArgs(trailingOnly = TRUE)
data.dir <- args[1]
samples <- list.files(data.dir, pattern = "*.rds", full.names=F)
results_file = args[2]
gmm_results_dir = args[3]

results_dir = dirname(results_file)
doublet_calling_output = paste0(results_dir,"/doublet_calling.txt")

###### column name in Seurat object with ground truth labels
ground_truth_column = "orig.ident"

# report performance of all tags or just the average
report_all_SBOs = FALSE

###### cluster_based demultiplexing parameters
cluster_based_lfc = 0.1
cluster_based_resolution = 1
mode = "auto"                     # "auto" determines logfc parameters automatically 
expected_doublet_rate = NULL  # if NULL, expected rate is estimated based on number of cells 
knn = 20

###### tools to be compared
norm.methods = c( "CLR_1", "CLR_2", "two_step")

demux_tools = c( "Seurat", "MULTIseqDemux","cluster_based")

combinations = expand.grid(norm.methods,demux_tools)
normalization_based = paste0(combinations$Var1,"_",combinations$Var2)

cellhashR_methods =  c("dropletutils", "bff_raw", "bff_cluster")

workflows = list("normalization_based" = normalization_based,"cellhashR_methods" = cellhashR_methods,
                 "demuxmix" = "demuxmix", "deMULTIplex2" = "deMULTIplex2", "GMM_Demux" = "GMM_Demux")

# 
McGinnis_2019_tag_mapping = list( "Bar1" = "HEK", "Bar2" = "HEK", "Bar3" = "MEF","Bar4" ="MEF",
                         "Bar5" = "Jurkats","Bar6" = "Jurkats","Bar7" = "Jurkats",
                         "Bar8" = "Jurkats","Bar9" = "Jurkats","Bar10" = "Jurkats",
                         "Bar11" = "Jurkats","Bar12" = "Jurkats", "Negative" = "Negative",
                         "Doublet"="Doublet")

for (sample in samples) {
  seurat_object = readRDS(paste0(data.dir,"/", sample))
  sample_name = stringr::str_remove_all(sample,pattern = ".rds")
  print(sample_name)
  
  ### A) Run demultiplexing tools that expect normalized counts as an input (HTOdemux, MULTIseqdemux and cluster_based)
  if("normalization_based" %in% names(workflows)){
    tmp = compare_normalization_methods(seurat_object, sbo_assay = "HTO",norm.methods = norm.methods, 
                                  demux_tools = demux_tools,mode = mode,expected_doublet_rate = expected_doublet_rate,
                                  lfc = cluster_based_lfc, resol = cluster_based_resolution, knn = knn)
    seurat_object = tmp
  }

  # ### B) Run demultiplexing tools that don’t expect normalized counts as an input
  # # 1) deMULTIplex2
  if("deMULTIplex2" %in% names(workflows)){
    tag_mtx <- seurat_object@assays$HTO@counts %>% t()
    
    tryCatch(expr = {
      demultiplex2_res <- demultiplexTags(tag_mtx, # Required, the tag count matrix from your experiment, can be either dense or sparse
                                          plot.path = "./", # Where to output a summary plot
                                          plot.name = sample_name, # text append to the name of the summary plot file
                                          plot.diagnostics =F,plot.umap = "none") # Whether to output diagnostics plots for each tag
      
      seurat_object$deMULTIplex2 = demultiplex2_res$final_assign
      seurat_object$deMULTIplex2[seurat_object$deMULTIplex2  == "negative"] = "Negative"
      seurat_object$deMULTIplex2[seurat_object$deMULTIplex2  == "multiplet"] = "Doublet"
      
    },error = function(e){ 
      
      seurat_object$deMULTIplex2 = "NA"
      
      })
  }
  # 
  # ## 2) demuxmix
  if("demuxmix" %in% names(workflows)){
    label = rep(NA,ncol(seurat_object))
    
    tryCatch(expr = {
      dmm = demuxmix(seurat_object@assays$HTO@counts %>% as.matrix(),model = "naive")
      classes <- dmmClassify(dmm)
      label = classes$HTO
      names(label) = rownames(classes)
      
      label[classes$Type %in% c("uncertain","negative")] = "Negative"
      label[classes$Type == "multiplet"] = "Doublet"
      seurat_object$demuxmix = label
    },error = function(e){ message("demuxmix failed") }, 
    finally =  { 
      seurat_object$demuxmix = label
      }
 )
  }
  ## 3) run tools supported by cellhashR
  
  if("cellhashR_methods" %in%  names(workflows)){
    
    df <- GenerateCellHashingCalls(barcodeMatrix = seurat_object@assays$HTO@counts,
                                   methods =cellhashR_methods,
                                   doTSNE = F,doHeatmap=F)
    
    rownames(df) = df$cellbarcode
    df =df[colnames(seurat_object),]
   
   for (tool in cellhashR_methods) {
     seurat_object@meta.data[,tool] = as.character(df[,tool])
   }
  }
  
  if("GMM_Demux" %in% names(workflows)){
    seurat_object$GMM_Demux = process_GMM_Demux_output(seurat_object, results_dir = paste0(gmm_results_dir, sample_name) )
  }
  
  all_workflows = workflows %>% unlist()
  
  if(sample_name == "McGinnis_2019"){
    for (w in all_workflows) {
      calls <- seurat_object@meta.data[,w]
      calls = sapply(calls, FUN = function(x) {McGinnis_2019_tag_mapping[[x]] })  
      seurat_object@meta.data[,w] = calls
      
    }
  }
  
  ### calculate classification accuracy metric (Fsore, recall and precision) and write results to a file
  write_results(seurat_object ,methods = all_workflows ,ground_truth_column = ground_truth_column,
                doublet_calling_output = doublet_calling_output,
                report_all_SBOs = report_all_SBOs,
                results_file = results_file ,sample = sample_name)
}






