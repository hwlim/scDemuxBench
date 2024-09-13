library(Seurat)
library(PRROC)
library(doParallel)
library(ComplexHeatmap)
library(dplyr)
library(stringr)
library(hashDemux)
#source("./demultiplexing.R")

# generate UMAP plot 
runUMAP = function(seurat_object,assay = "HTO",normalization_margin = 2){
  seurat_object = seurat_object %>% NormalizeData(assay = assay, normalization.method = "CLR",margin = normalization_margin) %>% 
    ScaleData() %>% RunPCA(features = rownames(seurat_object@assays[[assay]]) ,
                           npcs = nrow(seurat_object@assays[[assay]]) - 1) %>% 
    RunUMAP(dims = 1:(nrow(seurat_object@assays[[assay]]) - 1), seed.use = 1) 
  return(seurat_object)
}


# compute area under the Precision-recall curve (AUCPR) for each tag given ground truth (tag-positive and tag-negative cells)
separability_analysis <- function(seurat_object, tag, assay = "HTO", positive_cells,
                                               norm.methods = c("LogNormalize", "CLR_1", "CLR_2", "RC")){
  p = list()
  AUCPR = list()
  
  # using raw counts
  df = data.frame(class = seurat_object@meta.data[,positive_cells], expression = log1p(seurat_object[[assay]]@counts[tag,] ))
  # Compute roc
  AUCPR[["raw"]]  <- PRROC::pr.curve(scores.class0 = df$expression[df$class == paste0(tag,"+ cells")], scores.class1 = df$expression[df$class == paste0(tag,"- cells")])$auc.davis.goadrich
  
  #plot distribution
  p[["raw"]] = ggplot(df, aes(x = expression, fill = class)) +geom_density( alpha = 0.7,aes(y= after_stat(density) ))+
    xlab("log (raw counts + 1)")+ labs(fill = "classification") +
    theme(axis.title = element_text(size=13), legend.text =  element_text(size=13))
  
  p[["raw2"]] = ggplot(df, aes(x = expression, fill = "gray")) +geom_density( alpha = 0.7,aes(y= after_stat(density) ))+
    xlab("log (raw counts + 1)")+
    theme(axis.title = element_text(size=13))
  
  for (norm.method in norm.methods) {
    
    if(norm.method == "CLR_1"){
      seurat_object = seurat_object %>% NormalizeData(normalization.method = "CLR",
                                                      margin = 1 ,assay = assay)
    }else if(norm.method == "CLR_2"){
      seurat_object = seurat_object %>% NormalizeData(normalization.method = "CLR",
                                                      margin = 2 ,assay = assay)
    }else{
      seurat_object = seurat_object %>% NormalizeData(normalization.method = norm.method,assay = assay)
    }
    
    df = data.frame(class = seurat_object@meta.data[,positive_cells], expression = seurat_object[[assay]]@data[tag,] )

    p[[norm.method]] = ggplot(df, aes(x = expression, fill = class)) + geom_density(alpha = 0.7, aes(y= after_stat(density) )) +
      xlab("normalized counts")+ labs(fill = "classification")+
      theme(axis.title = element_text(size=13), legend.text =  element_text(size=13))
    
    # Compute AUC
    AUCPR[[norm.method]] <- PRROC::pr.curve(scores.class0 = df$expression[df$class == paste0(tag,"+ cells")], scores.class1 = df$expression[df$class == paste0(tag,"- cells")])$auc.davis.goadrich
    
  }
  return(list(histograms = p, AUCPR = AUCPR))
  
}

### generate simulated tag counts matrix
# Modified from https://github.com/Gartner-Lab/deMULTIplex2/blob/92130a626249194ef7f6c0b1aecad6ac210e258b/R/simulation.R to 
# 1) allow for different number of cells per tag
# 2) guarantee that one barcode's ambient contamination level is eual to `max.ambient.log`, 
# and other barcodes' ambient contamination levels are uniformly sampled from range [min.ambient.log , min.ambient.log]
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
evaluate_predictions <- function(seurat_object ,calls_column, ground_truth_column = "ground_truth") {
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
  recall <- c(recall, mean(recall)) #Add mean recall
  precision <- c(precision,mean(precision)) #Add mean precision
  
  names(f) <- c(labels, "Average")
  names(recall) <- c( labels,"Average")
  names(precision) <- c( labels,"Average")
  return(list(Fscore = f, recall = recall ,precision = precision))
}

## run different normalization methods followed by either Seurat's HTOdemux, Seurat's MULTIseqDemux or clustering-based
run_normalization_based_workflows = function(seurat_object,assay = "HTO",
                                   norm.methods = c("LogNormalize", "CLR_exp_bias", "CLR_seq_depth", "RC"),
                                   demux_tools = c( "Seurat", "MULTISeq","Clustering"))                      
{
  
  for (demux_tool in demux_tools) {
    for (norm.method in norm.methods) {
      
      label = rep(NA,ncol(seurat_object))
      
      tryCatch(expr = {    
        
        if(norm.method == "CLR_exp_bias"){
          seurat_object = seurat_object %>% NormalizeData(normalization.method = "CLR", 
                                                          margin = 1 ,assay = assay)
        }else if(norm.method == "CLR_seq_depth"){
          seurat_object = seurat_object %>% NormalizeData(normalization.method = "CLR", 
                                                          margin = 2 ,assay = assay)
        }else{
          seurat_object = seurat_object %>% NormalizeData(normalization.method = norm.method,assay = assay)
        }
  
        # demux 
        if(demux_tool == "Seurat"){
          
          tmp <- HTODemux(seurat_object, assay = assay)
          seurat_object$hash.ID <- tmp$hash.ID
          
        }else if(demux_tool == "Clustering"){
          
          tmp = clustering_based_demux(seurat_object, assay = assay ) 
          seurat_object$hash.ID = tmp$sampleBC
          seurat_object$hash.ID[tmp$classification == "Doublet" ] = "Doublet"
          
        }else if(demux_tool == "MULTISeq"){
           tmp = MULTIseqDemux(seurat_object,autoThresh = T)
            seurat_object$hash.ID = tmp$MULTI_ID
            }
        
        label = seurat_object$hash.ID
      
    },error = function(e){ message( paste0(demux_tool,"(",norm.method,") failed!")) }, 
    finally =  { 
      seurat_object@meta.data[, paste0(demux_tool,"(",norm.method,")")] = as.character(label)
    } )
      
    }
  }
  return(seurat_object)
}


# write classification accuracy metrics (i.e recall, precision and Fscore) and doublet calling results to a text file
write_results = function(seurat_object, methods,ground_truth_column , results_dir = "./", 
                         sample = "sample",per_tag= FALSE){
  
  for(method in methods){
    doublet_calls = doublet_calling(seurat_object,calls_column =  method,ground_truth_column)
    
    for (i in 1:length(names(doublet_calls))) {
      cat(sprintf("%s\t%s\t%s\t%f\n",sample,method,names(doublet_calls)[i],doublet_calls[i]), 
          file = paste0(results_dir,"doublet_calling.txt"), append = TRUE)
    }
    
    Fscore = evaluate_predictions(seurat_object,calls_column =  method,ground_truth_column)

        #write to file
    if(per_tag){
      for(tag in names(Fscore[[1]])){
        cat(sprintf("%s\t%s\t%s\t%s\t%f\n",sample,method,tag ,"Fscore" ,Fscore[[1]][tag]), 
            file = paste0(results_dir,"fscore_precision_recall.txt"), append = TRUE)
        cat(sprintf("%s\t%s\t%s\t%s\t%f\n",sample,method, tag,"Recall",Fscore[[2]][tag]), 
            file = paste0(results_dir,"fscore_precision_recall.txt"), append = TRUE)
        cat(sprintf("%s\t%s\t%s\t%s\t%f\n",sample,method,tag, "Precision",Fscore[[3]][tag]), 
            file = paste0(results_dir,"fscore_precision_recall.txt"), append = TRUE)
      }
    }else{
      tag = "Average"
      cat(sprintf("%s\t%s\t%s\t%s\t%f\n",sample,method,tag ,"Fscore" ,Fscore[[1]][tag]), 
          file = paste0(results_dir,"fscore_precision_recall.txt"), append = TRUE)
      cat(sprintf("%s\t%s\t%s\t%s\t%f\n",sample,method, tag,"Recall",Fscore[[2]][tag]), 
          file = paste0(results_dir,"fscore_precision_recall.txt"), append = TRUE)
      cat(sprintf("%s\t%s\t%s\t%s\t%f\n",sample,method,tag, "Precision",Fscore[[3]][tag]), 
          file = paste0(results_dir,"fscore_precision_recall.txt"), append = TRUE)
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


### investigate the effect of sequencing depth normalization on demultiplexing performance
seq_depth_effect <- function(seurat_object, assay = "HTO", tag = NULL,nbins=100, norm.method = "CLR_2",
                             ground_truth_column = NULL, predictions = NULL, pointSize = 0.5  )
{
  
  plots= list()
  
  mtrx = GetAssayData(object = seurat_object, assay = assay, slot = "counts") %>% as.matrix()
  df = mtrx %>% reshape2::melt()
  colnames(df) = c("Tag", "cell_id","counts")
  
  if(norm.method == "CLR_1"){
    seurat_object = seurat_object %>% NormalizeData(normalization.method = "CLR",
                                                    margin = 1 ,assay = assay)
  }else if(norm.method == "CLR_2"){
    seurat_object = seurat_object %>% NormalizeData(normalization.method = "CLR",
                                                    margin = 2 ,assay = assay)
  }else{
    seurat_object = seurat_object %>% NormalizeData(normalization.method = norm.method,assay = assay)
  }
  
  mtrx_norm = GetAssayData(object = seurat_object, assay = assay, slot = "data") %>% as.matrix()
  
  df2 = mtrx_norm %>% reshape2::melt()
  colnames(df2) = c("Tag", "cell_id","normalized_counts")
  
  df$normalized_counts = df2$normalized_counts
  
  # get read counts of a specific tag across all cells
  df3 = df[df$Tag == tag, ]
  rownames(df3) = df3$cell_id
  df3 = df3[colnames(seurat_object),]
  
  # add cells to buckets according to their raw/normalized counts of this tag
  df3 = df3 %>% mutate(normalized_percentile = ntile(normalized_counts, nbins),
                       raw_percentile = ntile(counts,nbins) )
  
  df3$diff = df3$normalized_percentile - df3$raw_percentile
  df3$abs_diff = abs(df3$normalized_percentile - df3$raw_percentile)
  
  plots[[1]] = ggplot(df3, aes( x=raw_percentile, y = normalized_percentile)) +
    geom_point(size = pointSize, color = "#56B4E9") + geom_abline(linetype ="dashed", color = "black", size = 1)+
    ylab("Percentile of normalized counts") +
    xlab("Percentile of raw counts")+
    theme(text = element_text(size = 12)) + theme_bw()
  
  plots[[2]] = ggplot(df3, aes( x=raw_percentile, y = normalized_percentile)) +
    geom_density_2d_filled() + geom_abline(linetype ="dashed", color = "red", size = 1)+
    ylab("Percentile of normalized counts") +
    xlab("Percentile of raw counts")+theme_classic()+
    theme(text = element_text(size = 12),legend.position = "none")
  
  if(!is.null(ground_truth_column) ){
    df3[,"ground_truth"] = seurat_object@meta.data[,ground_truth_column]
    
    df3$ground_truth2 =  df3$ground_truth
    df3$ground_truth2[df3$ground_truth2 != tag ] = "Others"
    df3$ground_truth2=as.factor(df3$ground_truth2)
    df3$ground_truth2 = factor(df3$ground_truth2, levels = c(tag, "Others"))
    
    plots[[5]] = ggplot(df3, aes( x=raw_percentile, y = normalized_percentile, color=ground_truth2 )) +
      geom_point(size = pointSize) + geom_abline(linetype ="dashed", color = "black", size = 1)+
      ylab("Percentile of normalized counts") +
      xlab("Percentile of raw counts")+
      theme(text = element_text(size = 12), legend.text=element_text(size=12) )+
      labs(color='Ground truth')  + theme_bw()
    
    if(!is.null(predictions) ){
      df3[,"predictions"] = seurat_object@meta.data[,predictions]
      
      df3$discordant = ifelse(df3$ground_truth != df3$predictions, "discordant", "concordant")
      df3$discordant2 = ifelse(df3$discordant == "discordant" , paste0(df3$ground_truth," -> ",df3$predictions), "concordant")
      
      sign_categories = c(paste0(tag," -> ", "Negative") , paste0("Negative"," -> ",tag) )
      
      df3$Concordnance = ifelse(df3$discordant2 %in% c(sign_categories,"concordant") ,df3$discordant2 , "other_discordant")
      df3$Concordnance = factor(df3$Concordnance , levels = c("concordant", "other_discordant",
                                                              paste0(tag," -> ", "Negative"),  paste0("Negative"," -> ",tag) ))
      plots[[3]] = ggplot(df3, aes( x=raw_percentile, y = normalized_percentile, color = Concordnance)) +
        geom_point(size = pointSize) +geom_abline(linetype ="dashed")+
        ylab("Percentile of normalized counts") +
        xlab("Percentile of raw counts")+
        theme(text = element_text(size = 12), legend.text=element_text(size=11) )+
        scale_color_manual(breaks =  c("concordant", "other_discordant",
                                       paste0(tag," -> ", "Negative"),  paste0("Negative"," -> ",tag) ),
                           values=c("#999999", "#E69F00", "#56B4E9", "red")) + theme_bw()
      
      plots[[4]] = ggplot(df3,  aes(y=diff, x = Concordnance,color =Concordnance) )+
        geom_boxplot(width = 0.4)+geom_hline(yintercept = 0,linetype = "dashed")+
        xlab("Discrepancy between ground truth and predictions") +
        ylab("Percentile of normalized counts - Percentile of raw counts")+theme_minimal()+
        theme(text = element_text(size = 12), axis.text.x =  element_text(angle=45, hjust=1, size = 11))
    }
  }
  return(list(plots = plots, data = df3, med_diff = median(df3$abs_diff)))
}

