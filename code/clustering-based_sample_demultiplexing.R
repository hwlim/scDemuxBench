library(Seurat)
library(dplyr)


# clustering-based sample demux
clustering_based_demux <- function(seurat_object,assay = "HTO",mode ="auto",
                                   expected_doublet_rate = NULL,resol = NULL,
                                   logfc.threshold = 0.5, adj_pVal =  0.05, knn = 20 )
{
  
  # build nearest-neighbor graph
  seurat_object <- Seurat::FindNeighbors(seurat_object, features = rownames(seurat_object), dims = NULL,
                                         assay = assay,k.param = knn )
  
  if(is.null(resol)){
    resol = ceiling(nrow(seurat_object) / 20)
  }
  
  # cluster cells
  seurat_object <- Seurat::FindClusters(seurat_object,
                                        graph.name = paste0(assay,"_snn") ,
                                        resolution = resol )
  # find cluster markers
  if(mode == "auto"){
    logfc.thresholds <- seq(0.05,0.5,0.05)
    
    doParallel::registerDoParallel(parallel::detectCores() -2)
    
    results <- foreach(logfc.threshold = logfc.thresholds ) %dopar% {
      out <- label_clusters(seurat_object,assay = assay,
                            logfc.threshold =logfc.threshold)
      prop.table(table(out$class))
    }
    dfs <- lapply(results, data.frame)
    df = dplyr::bind_rows(dfs)
    colnames(df) = c("classification" , "proportion")
    df$logfc.threshold = rep(logfc.thresholds, each = 3)
    
    if(is.null(expected_doublet_rate)){
      expected_doublet_rate = 8e-6 * ncol(seurat_object) * (nrow(seurat_object)-1)  / nrow(seurat_object)
    }
    
    df2 = df[df$classification == "Doublet",]
    index = which.min(abs(df2$proportion - expected_doublet_rate))
    threshold = df2$logfc.threshold[index]
    
  }else{
    threshold = logfc.threshold
  }
  labels = label_clusters(seurat_object,assay = assay,
                          logfc.threshold = threshold)
  return(labels)
  
}

## Assign a tag to each cluster based on marker tags
##
## output: a data frame with two columns,
##        1. sampleBC: assign a tag to each Singlet, tag1_tag2 in case of Doublets and Negative otherwise
##        2. class: classify cells to either Singlet, Doublet or Negative
label_clusters <- function(seurat_object,assay = "HTO",
                           logfc.threshold = 0.5, adj_pVal =  0.05 )
{
  # find cluster markers
  islet.markers <- FindAllMarkers(seurat_object, only.pos = TRUE,pseudocount.use = 0.1, # compatible with Seurat5
                                  assay = assay, verbose = F,
                                  logfc.threshold = logfc.threshold  )
  
  # marker expression plots
  if(nrow(islet.markers) == 0){
    df = data.frame(cell_barcode = colnames(seurat_object), sampleBC = "Negative", class = "Negative")
    rownames(df) <- df$cell_barcode
    return(df)
  }
  
  top_markers <- islet.markers %>% dplyr::filter(pct.1 == 1, p_val_adj < adj_pVal)  %>%
    dplyr::group_by(cluster) %>%
    dplyr::top_n(2, wt = avg_log2FC)
  
  ## assign cells to sample of origin
  cluster_barcode_dict <- top_markers %>% dplyr::select(c("cluster", "gene"))
  
  cell_cluster_dict <- seurat_object[["seurat_clusters"]]
  colnames(cell_cluster_dict) <- "cluster"
  cell_cluster_dict$cell_barcode <- rownames(cell_cluster_dict)
  
  cell_barcode_dict <- dplyr::left_join(cell_cluster_dict , cluster_barcode_dict , by = c("cluster") )
  
  cell_classification <- cell_barcode_dict %>% dplyr::group_by(cell_barcode) %>%
    dplyr::summarise(sampleBC = paste0(sort(gene), collapse = "_"))
  
  cell_classification$sampleBC[cell_classification$sampleBC == ""] = "Negative"
  
  cell_classification.global <- cell_classification %>%
    dplyr::mutate(class= if_else(sampleBC == "Negative" , "Negative",
                                 if_else(grepl("_",sampleBC) ,"Doublet" ,"Singlet")))
  
  cell_classification.global$class = factor(cell_classification.global$class,
                                            levels = c( "Doublet","Negative", "Singlet" ))
  rownames(cell_classification.global) <- cell_classification.global$cell_barcode
  
  cell_classification.global= cell_classification.global[colnames(seurat_object),]
  return( cell_classification.global)
  
}

