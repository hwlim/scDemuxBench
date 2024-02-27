library(Seurat)
library(dplyr)
source("functions.R")

data.dir = "~/projects/demux-benchmarking/results/all_datasets/seurat/"
results_dir = "../results/"

AUCPR = list()
samples <- list.files(data.dir, pattern = "*.rds", full.names=F)
samples=samples[samples!= "McGinnis_2019.rds"]

norm.methods = c("CLR_1", "CLR_2","RC", "LogNormalize")

for (sample in samples) {
  seurat_object = readRDS(paste0(data.dir,"/", sample))
  sample_name = stringr::str_remove_all(sample,pattern = ".rds")
  tmp = subset(seurat_object, subset = ground_truth != "Doublet")
  
  for (tag in rownames(seurat_object)) {
    tmp$pos_neg= ifelse(grepl(pattern = tag,x = tmp$ground_truth ),paste0(tag,"+ cells"),paste0(tag,"- cells"))
    
    res = separability_analysis(tmp, tag = tag,positive_cells = "pos_neg", 
                                             norm.methods = norm.methods)
    for(m in c("raw",norm.methods)) {
      AUCPR[[sample_name]][[tag]][[m]] = res$AUCPR[[m]] 
    }
    
  }
}

df = data.frame(AUCPR = AUCPR %>% unlist())
df$Dataset = lapply(str_split(rownames(df), pattern = "\\."),FUN = function(x) x[1] ) %>% as.character()
df$tag = lapply(str_split(rownames(df), pattern = "\\."),FUN = function(x) x[2] ) %>% as.character()
df$Normalization = lapply(str_split(rownames(df), pattern = "\\."),FUN = function(x) x[3] ) %>% as.character()

write.csv(df, file = paste0(results_dir, "separability_analysis.csv"),quote = F,row.names = F)
