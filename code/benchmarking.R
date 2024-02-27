library(Seurat)
library(hashDemux)
library(doParallel)
library(cellhashR)
library(demuxmix)
library(deMULTIplex2)
library(dplyr)

source("./functions.R")

set.seed(1)
###### load Seurat objects
data.dir <- "../results/simulated/"
samples <- list.files(data.dir, pattern = "*.rds", full.names=F)
results_dir = "../results/"
gmm_demux_results_dir = "../results/gmm/"

###### column name in Seurat object with ground truth labels
ground_truth_column = "ground_truth"

# report performance of all tags or just the average 
per_tag = FALSE

###### tools to be compared
norm.methods = c("LogNormalize", "CLR_exp_bias", "CLR_seq_depth", "RC")
demux_tools = c("Seurat", "MULTISeq","Clustering")

combinations = expand.grid(norm.methods,demux_tools)
normalization_based = paste0(combinations$Var2,"(",combinations$Var1,")")

cellhashR_methods =  c("dropletutils", "bff_raw", "bff_cluster")

#workflows = list("normalization_based" = normalization_based,"cellhashR_methods" = cellhashR_methods,
#                "demuxmix" = "demuxmix", "deMULTIplex2" = "deMULTIplex2", "GMM_Demux" = "GMM_Demux")

workflows = list("demuxmix" = "demuxmix")

# mapping tags in counts matrices to ground truth ID
McGinnis_2019_tag_mapping = list( "Bar1" = "HEK", "Bar2" = "HEK", "Bar3" = "MEF","Bar4" ="MEF",
                         "Bar5" = "Jurkats","Bar6" = "Jurkats","Bar7" = "Jurkats",
                         "Bar8" = "Jurkats","Bar9" = "Jurkats","Bar10" = "Jurkats",
                         "Bar11" = "Jurkats","Bar12" = "Jurkats", "Negative" = "Negative",
                         "Doublet"="Doublet")

for (sample in samples) {
  seurat_object = readRDS(paste0(data.dir,"/", sample))
  sample_name = stringr::str_remove_all(sample,pattern = ".rds")
  print(sample_name)
  
  ### A) Run demultiplexing tools that expect normalized counts as an input (HTOdemux, MULTIseqdemux and clustering_based)
  if("normalization_based" %in% names(workflows)){
    tmp = run_normalization_based_workflows(seurat_object, assay = "HTO",norm.methods = norm.methods, 
                                  demux_tools = demux_tools)
    seurat_object = tmp
  }

  # ### B) Run demultiplexing tools that don’t expect normalized counts as an input
  # # 1) deMULTIplex2
  if("deMULTIplex2" %in% names(workflows)){
    tag_mtx <- GetAssayData(seurat_object, slot = "counts", assay = "HTO") %>% as.matrix() %>% t()
    
    label = rep(NA,ncol(seurat_object))
    
    tryCatch(expr = {
      demultiplex2_res <- demultiplexTags(tag_mtx, # Required, the tag count matrix from your experiment, can be either dense or sparse
                                          plot.path = "./", # Where to output a summary plot
                                          plot.name = sample_name, # text append to the name of the summary plot file
                                          plot.diagnostics =F,plot.umap = "none") # Whether to output diagnostics plots for each tag
      
      label = demultiplex2_res$final_assign
      label[label == "negative"] = "Negative"
      label[label == "multiplet"] = "Doublet"

    },error = function(e){ message("deMULTIplex2 failed!") }, 
    finally =  { 
      seurat_object$deMULTIplex2 = as.character(label)
    } )
  }
  # 
  # ## 2) demuxmix
  if("demuxmix" %in% names(workflows)){
    label = rep(NA,ncol(seurat_object))
    
    tryCatch(expr = {
      dmm = demuxmix(GetAssayData(seurat_object, slot = "counts", assay = "HTO") %>% as.matrix(),model = "naive")
      classes <- dmmClassify(dmm)
      label = classes$HTO
      names(label) = rownames(classes)
      
      label[classes$Type %in% c("uncertain","negative")] = "Negative"
      label[classes$Type == "multiplet"] = "Doublet"
      seurat_object$demuxmix = label
    },error = function(e){ message("demuxmix failed!") }, 
    finally =  { 
      seurat_object$demuxmix = as.character(label)
      } )
  }
  ## 3) run tools supported by cellhashR
  
  if("cellhashR_methods" %in%  names(workflows)){
    
    df <- GenerateCellHashingCalls(barcodeMatrix = GetAssayData(seurat_object, slot = "counts", assay = "HTO"),
                                   methods =cellhashR_methods,
                                   doTSNE = F,doHeatmap=F)
    
    rownames(df) = df$cellbarcode
    df =df[colnames(seurat_object),]
   
   for (tool in cellhashR_methods) {
     if (tool == "dropletutils") {
       seurat_object$hashedDrops = as.character(df[,tool])
       
     }else{
       seurat_object@meta.data[,tool] = as.character(df[,tool])
     }
   }
  }
  
  if("GMM_Demux" %in% names(workflows)){
    seurat_object$GMM_Demux = process_GMM_Demux_output(seurat_object, results_dir = paste0(gmm_demux_results_dir, sample_name) )
  }
  
  all_workflows = workflows %>% unlist()
  all_workflows[all_workflows == "dropletutils"] = "hashedDrops"
  
  # for some datasets do tag mapping
  if(sample_name == "McGinnis_2019"){
    for (w in all_workflows) {
      calls <- seurat_object@meta.data[,w]
      calls = sapply(calls, FUN = function(x) {McGinnis_2019_tag_mapping[[x]] })  
      seurat_object@meta.data[,w] = calls
      
    }
  }
  
  ### calculate classification accuracy metric (Fsore, recall and precision) and doublet classification results and write results to a file
  write_results(seurat_object ,methods = all_workflows ,ground_truth_column = ground_truth_column,
                per_tag = per_tag, results_dir = results_dir ,sample = sample_name)
}






