library(Seurat)
library(ComplexHeatmap)
library(dplyr)
library(stringr)

## generate Heatmap from columns in a dataframe
dataframeToHeatmap = function(df, rows.var, columns.var,value.var = "value",cluster_rows=F,cluster_columns=T,
                              legend.name="F-score",min_value_heatmap=0,max_value_heatmap = 1,
                              text.fontsize = 12, row_names.fontsize =14,column_names.fontsize =14,
                              hcl_palette = "BlueRed 2"    # choose one from colorspace::hcl_palettes()
                              ){
  
  mtrx = df %>% reshape2::dcast(formula = paste0(rows.var, "~", columns.var),
                                value.var = value.var)
  rownames(mtrx) = mtrx[,rows.var]
  mtrx = mtrx[,c(2:ncol(mtrx))]
  mtrx = mtrx %>% as.matrix()
  
  
  col_fun = circlize::colorRamp2(breaks = seq(min_value_heatmap,max_value_heatmap, 0.1 * (max_value_heatmap- min_value_heatmap)),
                                 hcl_palette = hcl_palette, reverse = F)
  ComplexHeatmap::Heatmap(mtrx %>% as.matrix(), col =col_fun, cluster_rows = cluster_rows, cluster_columns = cluster_columns,
                          row_names_gp = gpar(fontsize = row_names.fontsize) ,
                          column_names_gp = gpar(fontsize = column_names.fontsize) ,
                          name = legend.name,
                          cell_fun = function(j, i, x, y, width, height, fill) 
                          {grid.text(sprintf("%0.2f", mtrx[i, j]), x, y, gp = gpar(fontsize = text.fontsize,
                                                                                   col = "black"))},
                          row_names_side = "right", column_names_side = "bottom" )
}


# Modified from https://github.com/Gartner-Lab/deMULTIplex2/blob/1b2942afeae24dd48fb747e81aa4fb65d2bf8ad2/R/simulation.R to
# 1) allow for different number of cells per tag
# 2) guarantee that one barcode's ambient contamination level is equal to `max.ambient.log`, 
# and other barcodes' ambient contamination levels are uniformly sampled from range [min.ambient.log , max.ambient.log]
simulateTags2 <- function(min.n.cell = 10,
                          max.n.cell = 1000,
                          n.bc = 10,
                          seed = 1,
                          nb.theta = 10,
                          min.size.log = 1,
                          max.size.log = 7,
                          min.ambient.log = 0,
                          max.ambient.log = 5,
                          b0 = -3,
                          doublet.rate = .1,
                          separate.sim = TRUE,
                          cell.contam = T,
                          ambient.contam = T,
                          return.all = FALSE) {
  set.seed(seed)
  bcs = paste0("bc", 1:n.bc)
  
  # modified to guarantee that one barcode's ambient contamination level is max.ambient.log, 
  # other barcodes' ambient contamination levels are uniformly sampled from range [min.ambient.log , min.ambient.log]
  bc.contam.count = round(exp(runif((n.bc-1), min = min.ambient.log, max = max.ambient.log)))
  bc.contam.count = c(bc.contam.count, round(exp(max.ambient.log)))
  names(bc.contam.count) <- bcs
  
  # 1. Simulation of true initial staining
  cell.true.umi <- list()
  for(bc in bcs) {
    n.cell = runif(1, min = min.n.cell, max = max.n.cell)
    log.bc.umi = runif(n.cell, min = min.size.log, max = max.size.log)
    df = data.frame(ceiling(exp(log.bc.umi)))
    rownames(df) <- paste0(bc, "_", 1:nrow(df))
    colnames(df) = bc
    cell.true.umi[[bc]] <- df
  }
  cell.true.umi.mtx = do.call(plyr::rbind.fill, cell.true.umi)
  rownames(cell.true.umi.mtx) <- as.character(unlist(lapply(cell.true.umi, function(x) rownames(x))))
  cell.true.umi.mtx[is.na(cell.true.umi.mtx)] = 0
  
  if(separate.sim ) {
    # 2. Simulation of floating barcodes contaminating cell surface
    cell.surface.size= log(rowSums(cell.true.umi.mtx)) # Note here direct sum because 0 on other entries
    
    if(cell.contam) {
      mu = exp(1*cell.surface.size + b0) # Assume same contamination level for all barcodes
      cell.contam.umi.mtx = sapply(bcs, function(bc) {
        MASS::rnegbin(length(cell.surface.size), mu = mu, theta = nb.theta)
      })
    } else {
      cell.contam.umi.mtx = cell.true.umi.mtx *0
    }
    
    # 3. Simulation of floating barcodes contaminating ambients
    if(ambient.contam) {
      ambient.contam.umi.mtx = sapply(bcs, function(bc) {
        MASS::rnegbin(n = length(cell.surface.size), mu = bc.contam.count[bc], theta = nb.theta)
      })
    } else {
      ambient.contam.umi.mtx = cell.true.umi.mtx *0
    }
    
    final.umi.mtx = cell.true.umi.mtx + cell.contam.umi.mtx + ambient.contam.umi.mtx
  } else {
    ## Combine 2 and 3 into a single NB? This can lead to a simpler model
    cell.surface.size= log(rowSums(cell.true.umi.mtx))
    mu_c = round(exp(1*cell.surface.size + b0))
    mu_list <- sapply(bcs, function(bc) {
      mu_c + bc.contam.count[bc]
    })
    cbn.contam.umi.mtx = sapply(bcs, function(bc) {
      rnegbin(length(cell.surface.size), mu = mu_list[,bc], theta = nb.theta)
    })
    final.umi.mtx = cell.true.umi.mtx + cbn.contam.umi.mtx
  }
  
  # Simulate doublet, inside cell.true.umi.mtx? TO BE DETERMINED
  if(doublet.rate > 0) {
    ndoublet = nrow(final.umi.mtx) * doublet.rate
    doub.idx1 = sample(1:nrow(final.umi.mtx), ndoublet)
    doub.idx2 = sample(1:nrow(final.umi.mtx), ndoublet)
    # doub.umi.mtx = final.umi.mtx[doub.idx1, ] + final.umi.mtx[doub.idx2, ]
    ambient.contam.idx = sample(1:nrow(ambient.contam.umi.mtx), ndoublet)
    doub.umi.mtx = cell.true.umi.mtx[doub.idx1, ] + cell.contam.umi.mtx[doub.idx1, ] +
      cell.true.umi.mtx[doub.idx2, ] + cell.contam.umi.mtx[doub.idx2, ] + ambient.contam.umi.mtx[ambient.contam.idx, ]
    rownames(doub.umi.mtx) <- paste0('doublet_', 1:nrow(doub.umi.mtx), "|", rownames(final.umi.mtx)[doub.idx1], "_", rownames(final.umi.mtx)[doub.idx2])
    final.umi.mtx <- rbind(final.umi.mtx, doub.umi.mtx)
  }
  
  if(return.all) {
    if(separate.sim) {
      res <- list(
        cell.true.umi.mtx = cell.true.umi.mtx,
        cell.contam.umi.mtx = cell.contam.umi.mtx,
        ambient.contam.umi.mtx = ambient.contam.umi.mtx,
        final.umi.mtx = final.umi.mtx
      )
    } else {
      res <- list(
        cell.true.umi.mtx = cell.true.umi.mtx,
        cbn.contam.umi.mtx = cbn.contam.umi.mtx,
        final.umi.mtx = final.umi.mtx
      )
    }
  } else {
    res <- final.umi.mtx
  }
  return(list(counts = res, bc.contam.count = bc.contam.count))
}


# centered-log ratio transformation
# from https://github.com/satijalab/seurat/blob/master/R/preprocessing.R
clr_function <- function(x) {
  return(log1p(x = x / (exp(x = sum(log1p(x = x[x > 0]), na.rm = TRUE) / length(x = x)))))
}


# two step normalization: first, relative count normalization and then per-tag CLR transformation
two_step = function(mtrx){
  mtrx2 = apply(mtrx, MARGIN = 2, FUN = function(x) (x/ sum(x)) * 10**4)
  mtrx2 = t(apply(mtrx2 , MARGIN = 1, FUN = clr_function))
  return(mtrx2)
}


## compute the proportion of true doublets predicted as Doublets, Singlets and Negatives
doublet_calling <- function(seurat_object ,calls_column, ground_truth_column){
  calls <- seurat_object@meta.data[,calls_column]
  calls[calls == "negative"] = "Negative"
  calls[calls %in% c("doublet", "multiplet") ] = "Doublet"
  calls[!calls %in% c("Doublet", "Negative") ] = "Singlet"
  calls = factor(calls, levels = c("Doublet","Singlet" ,"Negative") )
  ground_truth <- seurat_object@meta.data[,ground_truth_column]
  doublets = calls[ground_truth == "Doublet"]
  prop.table(table(doublets))
}


## calculates Recall, Precision, and Fscore metrics for classification performance
# modified from "calculate_HTO_fscore" function @ https://github.com/Oshlack/hashtag-demux-paper/blob/main/analysis/BAL_analysis.Rmd
evaluate_predictions <- function(seurat_object ,calls_column, ground_truth_column) {
  calls <- seurat_object@meta.data[,calls_column]
  ground_truth <- seurat_object@meta.data[,ground_truth_column]
  f <- NULL
  recall <- NULL
  precision <- NULL
  labels = unique(as.character(ground_truth))
  labels = labels[! stringr::str_to_lower(labels) %in% c("doublet", "negative", "multiplet") ]
  for (tag in labels) {
    tp <- sum(calls == tag & ground_truth == tag) #True positive rate
    fp <- sum(calls == tag & ground_truth != tag) #False positive rate
    fn <- sum(calls != tag & ground_truth == tag) #False negative rate
    
    f <- c(f, tp / (tp + 0.5 * (fp + fn)))
    recall<- c(recall, tp / (tp + fn) )
    precision<- c(precision, tp / (tp + fp) )
  }
  
  f <- c(f,mean(f)) #Add mean F score
  recall <- c(recall, mean(recall)) #Add mean F score
  precision <- c(precision,mean(precision)) #Add mean F score
  
  names(f) <- c(labels, "Average")
  names(recall) <- c( labels,"Average")
  names(precision) <- c( labels,"Average")
  return(list(Fscore = f, recall = recall ,precision = precision))
}

## run different normalization methods followed by either Seurat's HTOdemux, Seurat's MULTIseqDemux or clustering-based
compare_normalization_methods = function(seurat_object,sbo_assay = "HTO",
                                   norm.methods = c("LogNormalize", "CLR_1", "CLR_2", "RC", "two_step","SCTransform"),
                                   demux_tools = c( "Seurat", "MULTIseqDemux","cluster_based"),
                                  
                                    # cluster_based demultiplexing parameters
                                   expected_doublet_rate = NULL, knn = 20,
                                   lfc = 0.1, resol = 1,  
                                   mode = "auto")                      
{
  
  for (demux_tool in demux_tools) {
    for (norm.method in norm.methods) {
      
      if(norm.method == "CLR_1"){
        seurat_object = seurat_object %>% NormalizeData(normalization.method = "CLR", 
                                                        margin = 1 ,assay = sbo_assay)
      }else if(norm.method == "CLR_2"){
        seurat_object = seurat_object %>% NormalizeData(normalization.method = "CLR", 
                                                        margin = 2 ,assay = sbo_assay)
      }else if(norm.method == "SCTransform"){
        seurat_object = seurat_object %>% SCTransform(assay = sbo_assay)
        
      }else if(norm.method == "two_step"){
        mtrx = seurat_object[[sbo_assay]]@counts %>% as.matrix()
        seurat_object[[sbo_assay]]@data = two_step(mtrx)
        
      }else{
        seurat_object = seurat_object %>% NormalizeData(normalization.method = norm.method,assay = sbo_assay)
      }

      # demux 
      if(demux_tool == "Seurat"){
        if(norm.method == "SCTransform"){
          tmp <- HTODemux(seurat_object, assay = "SCT" ,positive.quantile = 0.99)
        }else{
          tmp <- HTODemux(seurat_object, assay = sbo_assay ,positive.quantile = 0.99)
        }
        
        seurat_object$hash.ID <- tmp$hash.ID
        
      }else if(demux_tool == "cluster_based"){
        
        if(norm.method == "SCTransform"){
          
          tmp = clustering_based_demux(seurat_object, assay = "SCT", resol=resol,knn = knn,
                                                         expected_doublet_rate= expected_doublet_rate,
                                                         logfc.threshold = lfc, mode = mode)        
        }else{
          tmp = clustering_based_demux(seurat_object, assay = sbo_assay,resol=resol,knn = knn,
                                                         expected_doublet_rate=expected_doublet_rate,
                                                         logfc.threshold = lfc, mode = mode) 
        }
        
        seurat_object$hash.ID = tmp$sampleBC
        seurat_object$hash.ID[tmp$class == "Doublet" ] = "Doublet"
        
      }else if(demux_tool == "MULTIseqDemux"){
        
        if(norm.method == "SCTransform"){
          seurat_object = MULTIseqDemux(seurat_object,assay = "SCT",autoThresh = T)
          seurat_object$hash.ID = seurat_object$MULTI_ID
          
        }else{
          seurat_object = MULTIseqDemux(seurat_object,autoThresh = T)
          seurat_object$hash.ID = seurat_object$MULTI_ID
        }
        
      }
      
      #save predictions
      seurat_object@meta.data[, paste0(norm.method,"_",demux_tool)] = as.character(seurat_object$hash.ID)
    }
  }
  return(seurat_object)
}


# write classification accuracy metrics (i.e recall, precision and Fscore) and doublet calling results to a text file
write_results = function(seurat_object, methods,ground_truth_column , results_file = "./results.txt", 
                         doublet_calling_output = "./doublet_calling.txt.txt",
                         sample = "sample",report_all_SBOs= FALSE){
  for(method in methods){
    doublet_calls = doublet_calling(seurat_object,calls_column =  method,ground_truth_column)
    
    for (i in 1:length(names(doublet_calls))) {
      cat(sprintf("%s\t%s\t%s\t%f\n",sample,method,names(doublet_calls)[i],doublet_calls[i]), 
          file = doublet_calling_output, append = TRUE)
    }
    
    Fscore = evaluate_predictions(seurat_object,calls_column =  method,ground_truth_column)

        #write to file
    if(report_all_SBOs){
      for(sbo in names(Fscore[[1]])){
        cat(sprintf("%s\t%s\t%s\t%s\t%f\n",sample,method,sbo ,"Fscore" ,Fscore[[1]][sbo]), 
            file = results_file, append = TRUE)
        cat(sprintf("%s\t%s\t%s\t%s\t%f\n",sample,method, sbo,"Recall",Fscore[[2]][sbo]), 
            file = results_file, append = TRUE)
        cat(sprintf("%s\t%s\t%s\t%s\t%f\n",sample,method,sbo, "Precision",Fscore[[3]][sbo]), 
            file = results_file, append = TRUE)
      }
    }else{
      sbo = "Average"
      cat(sprintf("%s\t%s\t%s\t%s\t%f\n",sample,method,sbo ,"Fscore" ,Fscore[[1]][sbo]), 
          file = results_file, append = TRUE)
      cat(sprintf("%s\t%s\t%s\t%s\t%f\n",sample,method, sbo,"Recall",Fscore[[2]][sbo]), 
          file = results_file, append = TRUE)
      cat(sprintf("%s\t%s\t%s\t%s\t%f\n",sample,method,sbo, "Precision",Fscore[[3]][sbo]), 
          file = results_file, append = TRUE)
    }
    
  }
  
}

# postprocess GMM-Demux output
process_GMM_Demux_output = function(seurat_object, results_dir ){
  gmm = read.csv(paste0(results_dir,"/GMM_full.csv"), row.names = 1)
  gmm$cell_barcode = rownames(gmm)
  
  config = read.csv(paste0(results_dir,"/GMM_full.config"), header = F,strip.white = T)
  colnames(config)  = c("Cluster_id", "class")
  
  df = dplyr::inner_join(gmm,config)
  
  df$calls = df$class
  df$calls[df$class == "negative"] = "Negative"
  df$calls[!df$class %in%  c("Negative", rownames(seurat_object))] = "Doublet"
  
  rownames(df)= df$cell_barcode
  df = df[colnames(seurat_object),]
  
  return(as.character(df$calls))
}




