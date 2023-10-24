library(Seurat)
library(openxlsx)
library(dplyr)

# increase download timeout
options(timeout=1000)

# This script provides the instructions to download, filter and save 12 publicly-available datasets
# each tag counts matrix along with ground truth are saved as Seurat object


# folder to save downloaded raw data
data.dir = "../datasets/raw_data//"

# folder with created Seurat objects
seurat.dir = "../datasets/seurat_objects/"

########## Mylka (2022)##########################################################
#citation: Mylka et al. "Comparative analysis of antibody-and lipid-based multiplexing methods for single-cell RNA-seq." Genome Biology 23.1 (2022): 1-21.

## download data
samples = c("CMO_nuclei","LMO_custom_cells","LMO_MULTISeq_cells","TotalSeqA_cells", "TotalSeqA_cells_rep2",
            "TotalSeqA_nuclei", "TotalSeqC_cells")

for (sample_name in samples) {
  sample.dir = paste0(data.dir,sample_name,"/")
  dir.create(sample.dir)
  
  #count matrix file
  download.file(paste0("https://www.ebi.ac.uk/biostudies/files/E-MTAB-9964/matrix_",sample_name,".mtx.gz"),
                destfile = paste0(sample.dir,"/matrix.mtx.gz"))
  #barcodes file
  download.file(paste0("https://www.ebi.ac.uk/biostudies/files/E-MTAB-9964/barcodes_",sample_name,".tsv.gz"),
                destfile = paste0(sample.dir,"/barcodes.tsv.gz"))
  #features file
  download.file(paste0("https://www.ebi.ac.uk/biostudies/files/E-MTAB-9964/features_",sample_name,".tsv.gz"),
                destfile = paste0(sample.dir,"/features.tsv.gz"))
  
  # ground truth file
  download.file(paste0("https://ftp.ebi.ac.uk/biostudies/nfs/E-MTAB-/964/E-MTAB-9964/Files/",sample_name
                       ,"_freemuxlet_MULTI_HTO_GMM_annotations.csv"), destfile = paste0(sample.dir,"/ground_truth.csv") )
  
  # preprocess ground truth
  ground.truth = read.table( paste0(sample.dir,"/ground_truth.csv"), sep = "\t", header = T)
  ground.truth$freemuxlet_final = ground.truth$freemuxlet_DROPLET.TYPE
  ground.truth$freemuxlet_final[ground.truth$freemuxlet_DROPLET.TYPE == "Singlet"] = ground.truth$freemuxlet_BEST.GUESS[ground.truth$freemuxlet_DROPLET.TYPE == "Singlet"]
  ground.truth$freemuxlet_final[ground.truth$freemuxlet_final == "AMB"] = "Negative"
  table(ground.truth$freemuxlet_final)
  
  # read tag counts matrix
  raw_mtrx = Read10X(sample.dir ,  strip.suffix = T)$`Antibody Capture` 
  rownames(raw_mtrx)
  
  # update tag names to be consistent with ground truth
  xx = table(ground.truth$MULTI_ID[!ground.truth$MULTI_ID %in% c("Doublet","Negative")], ground.truth$MULTI_classification[!ground.truth$MULTI_ID %in% c("Doublet","Negative")])
  xx
  rownames(raw_mtrx) = apply(xx,MARGIN = 2,FUN = function(x) names(x)[x == max(x)] )
  
  raw_mtrx=raw_mtrx[,ground.truth$cell]
  
  valid_tags = rownames(raw_mtrx)[rowSums(raw_mtrx > 0) > 30]
  raw_mtrx=raw_mtrx[valid_tags,]
  
  # create Seurat object
  seurat_object = CreateSeuratObject(counts = raw_mtrx,assay = "HTO",)
  
  # filter cells based on total tag counts
  seurat_object = seurat_object %>% subset(subset = nCount_HTO > 0)
  
  # add ground truth column to Seurat object meta.data
  ground.truth=ground.truth[ground.truth$cell %in% colnames(seurat_object),]
  seurat_object$ground_truth = ground.truth$freemuxlet_final
  
  saveRDS(seurat_object, file = paste0(seurat.dir,sample_name,".rds"))
}

########## Howitt (2022)#########################################################
#citation: Howitt et al. "Benchmarking single-cell hashtag oligo demultiplexing methods." NAR Genomics and Bioinformatics 5.4 (2023)
# data from 
data.source = "https://raw.githubusercontent.com/Oshlack/hashtag-demux-paper/f91952b926159d245a954796a16fa398f58e56c3/data/"

# donor-tag mapping
donor_list_batch1 <- list("donor_A" = "Human-HTO-3", "donor_B" = "Human-HTO-1", "donor_C" = "Human-HTO-4", "donor_D" = "Human-HTO-2",
                          "donor_E" = "Human-HTO-8", "donor_F" = "Human-HTO-6", "donor_G" = "Human-HTO-5", "donor_H" = "Human-HTO-7", 
                          "doublet" = "Doublet", "unassigned" = "Negative")
donor_list_batch2 <- list("donor_A" = "Human-HTO-7", "donor_B" = "Human-HTO-13", "donor_C" = "Human-HTO-15", 
                          "donor_D" = "Human-HTO-12", "donor_E" = "Human-HTO-6", "donor_F" = "Human-HTO-9", 
                          "donor_G" = "Human-HTO-14", "donor_H" = "Human-HTO-10","doublet" = "Doublet",
                          "unassigned" = "Negative")
donor_list_batch3 <- list("donor_A" = "Human-HTO-6", "donor_B" = "Human-HTO-10", "donor_C" = "Human-HTO-14", "donor_D" = "Human-HTO-13",
                          "donor_E" = "Human-HTO-15", "donor_F" = "Human-HTO-7", "donor_G" = "Human-HTO-12", "donor_H" = "Human-HTO-9",
                          "doublet" = "Doublet", "unassigned" = "Negative")


donor_LMO_list <- list("donor_A" = "BC1", "donor_B" = "BC2", "donor_C" = "BC3", 
                       "Doublet" = "Doublet", "Negative" = "Negative")

batches = paste0("batch",c(1:3))
captures = c("_c1", "_c2")
samples = paste0(rep(batches, each =2),rep(captures, 3))

samples = c(samples, "lmo")
list.obs = list()

for (sample_name in samples) {
  sample.dir = paste0(data.dir,sample_name,"/")
  dir.create(sample.dir)
  
  #count matrix file
  if(sample_name == "lmo"){
  download.file(paste0(data.source,sample_name,"_counts.csv"),
                destfile = paste0(sample.dir,"/hto_counts.csv"))
  }else{
    download.file(paste0(data.source,sample_name,"_hto_counts.csv"),
                  destfile = paste0(sample.dir,"/hto_counts.csv"))
  }
  
  # ground truth file
  download.file(paste0(data.source,sample_name,"_donors.csv"),
                destfile = paste0(sample.dir,"/donors.csv"))
  
  raw_mtrx = read.csv(paste0(sample.dir,"/hto_counts.csv"),
                        check.names = FALSE, row.names = 1)

  ground.truth = read.table(paste0(sample.dir,"/donors.csv"), header = T, sep = ",")
  colnames(ground.truth) = c("cell_ID", "label")
  rownames(ground.truth) = ground.truth$cell_ID
  
  # create Seurat object
  seurat_object = CreateSeuratObject(counts = raw_mtrx,assay = "HTO",)
  
  # filter cells based on total tag counts
  seurat_object = seurat_object %>% subset(subset = nCount_HTO > 0)
  
  # add ground truth column to Seurat object meta.data
  ground.truth=ground.truth[ground.truth$cell_ID %in% colnames(seurat_object),]
  seurat_object$genetic_donor = ground.truth$label
  
  if(substr(sample_name,start = 6,6) == "3" ){            # sample_name[6] is batch number
    donor_list  = donor_list_batch3
  }else if (substr(sample_name,start = 6,6) == "2" ){
    donor_list  = donor_list_batch2
  }else if (substr(sample_name,start = 6,6) == "1" ){
    donor_list  = donor_list_batch1
  }else{
    donor_list  = donor_LMO_list
  }
  
  seurat_object$ground_truth = sapply(seurat_object$genetic_donor, 
                                      FUN = function(x) {donor_list[[x]] })
  list.obs[[sample_name]] = seurat_object
}

saveRDS( merge(list.obs[[1]], list.obs[[2]]), file = paste0(seurat.dir,"BAL1",".rds"))
saveRDS( merge(list.obs[[3]], list.obs[[4]]), file = paste0(seurat.dir,"BAL2",".rds"))
saveRDS( merge(list.obs[[5]], list.obs[[6]]), file = paste0(seurat.dir,"BAL3",".rds"))
saveRDS( list.obs[[7]], file = paste0(seurat.dir,"lung_cell_line",".rds"))


########## McGinnis (2019)##########################################################
# citation: McGinnis et al. "MULTI-seq: sample multiplexing for single-cell RNA sequencing using lipid-tagged indices." Nature methods 16.7 (2019): 619-626.
sample_name = "McGinnis_2019"
sample.dir = paste0(data.dir,sample_name,"/")
dir.create(sample.dir)

# download ground truth
download.file("https://static-content.springer.com/esm/art%3A10.1038%2Fs41592-019-0433-8/MediaObjects/41592_2019_433_MOESM3_ESM.xlsx",
              destfile = paste0(sample.dir,"/ground_truth.xlsx"))

human_cells = openxlsx::read.xlsx(paste0(sample.dir,"/ground_truth.xlsx"), sheet = 2, startRow = 2)
mouse_cells = openxlsx::read.xlsx(paste0(sample.dir,"/ground_truth.xlsx"), sheet = 3, startRow = 2)

all_cells = dplyr::bind_rows(human_cells, mouse_cells)
all_cells = all_cells[all_cells$nUMI_Bar > 0, ]

raw_mtrx = all_cells[!duplicated( all_cells$CellID) , ]
rownames(raw_mtrx) = raw_mtrx$CellID
raw_mtrx = dplyr::select(raw_mtrx, starts_with("Bar")) %>% as.matrix() %>% t()

ground.truth = all_cells[,c("CellID","CellType")]
ground.truth = ground.truth %>% group_by(CellID) %>% summarise(CellType=paste0(CellType, collapse = ","))
rownames(ground.truth) = ground.truth$CellID

seurat_object = CreateSeuratObject(counts = raw_mtrx,assay = "HTO",)
seurat_object = seurat_object %>% subset(subset = nCount_HTO > 0)

ground.truth = ground.truth[colnames(seurat_object),]

seurat_object$ground_truth = ground.truth$CellType
seurat_object$ground_truth[grepl(pattern = ",", x = seurat_object$ground_truth)] = "Doublet"
table(seurat_object$ground_truth)
saveRDS( seurat_object, file = paste0(seurat.dir,sample_name,".rds"))





