#data:https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE272922
#=================================
#Load required Packages
#=================================
library(Seurat)
library(ggplot2)
library(tidyverse)
library(gridExtra)
library(patchwork)
library(SeuratObject)
install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(paste("Installing package:", pkg))
    install.packages(pkg, dependencies = TRUE)
  }
  library(pkg, character.only = TRUE)
}
#BiocManager::install("Seurat")

#install.packages("vctrs")
library(celldex)
library(SingleR)
#update.packages(ask = FALSE, checkConfluence = TRUE)
#remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')
#========================================================
#set paths and Sample info
#========================================================
getwd()
base_path = "D:/BI_prj/scrna_proj/neutrophils_BC/GSE272922_RAW"
setwd("D:/BI_prj/scrna_proj/neutrophils_BC")
sample_files = list.files(path = base_path, pattern = "matrix\\.mtx(\\.gz)?$"
, full.names = FALSE)
sample_files = sub("_matrix\\.mtx(\\.gz)?$", "", sample_files)

#========================================================
#Read Each sample and create seurat object
#========================================================

get_file_path = function(base,sample,suffix){
  gz_file = file.path(base, paste0(sample,"_",suffix,".gz"))
  unc_file = file.path(base, paste0(sample,"_",suffix))
  
  if (file.exists(gz_file)){
    return(gz_file)
  } else if (file.exists(unc_file)){
    return(unc_file)
  } else {
    return(NA)
  }
}


for (sample in sample_files){
  counts = ReadMtx(
    mtx = get_file_path(base_path,sample,"matrix.mtx"),
    features = get_file_path(base_path,sample,"features.tsv"),
    cells = get_file_path(base_path,sample,"barcodes.tsv")
  )
  seurat_obj = CreateSeuratObject(counts = counts, project = "MBC", min.cells = 3, min.features = 200)
  sample<-gsub("GSM[0-9]+_","", sample)
  assign(sample,seurat_obj)
  
}
#========================================================
#Merging datasets
#========================================================
ls()

data_merged = merge(Control,
             y = list(Phenformin,PolyIC,Combination),
             add.cell.ids = c("Control", "Phenformin", "PolyIC", "Combination"),
             project = "MBC"
)
View(data_merged@meta.data)
rm("Control", "Phenformin", "PolyIC", "Combination")
gc()
#Merging sample data
data_merged$sample <-rownames(data_merged@meta.data)
data_merged@meta.data<- separate(data_merged@meta.data, col = 'sample', into = c('Treatment','Barcode'),
         sep = '_')

#low quality bases - high mit contamination
data_merged[["percent.mt"]] <- PercentageFeatureSet(data_merged, pattern = "^mt-") #mouse is ^mt-

sum(grepl("^mt-", rownames(data_merged)))
head(rownames(data_merged), 20)

#========================================================
#QC and filtering https://www.sc-best-practices.org/preprocessing_visualization/quality_control.html
#========================================================
#summary table
qc_log <- data.frame(
  Step = c("Before", "Manual", "SoupX", "Doublet Finder"),
  Removed = 0,
  Remaining = 0
)

#Count cell before filtering
n_before <- ncol(data_merged)
qc_log[qc_log$Step == "Before", "Remaining"] <- n_before

QCp1<-VlnPlot(data_merged, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave("results/QC_vlnplot.png", bg = "white", plot = QCp1, width = 8, height = 6, dpi = 300)

QCp2<-FeatureScatter(data_merged, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")+
  geom_smooth(method = 'lm')
ggsave("results/QC_scatterplot.png", bg = "white", plot = QCp2, width = 8, height = 6, dpi = 300)

quantile(data_merged$nFeature_RNA, probs = c(0.01, 0.99))
quantile(data_merged$nCount_RNA, probs = c(0.01, 0.99))
quantile(data_merged$percent.mt, probs = c(0.01, 0.99))

# Filtering data using thresholds
data_filtered <- subset(
  data_merged,
  subset = nFeature_RNA > 200 & percent.mt < 10)

#Count cell after manual filtering
n_man <- ncol(data_filtered)
rem_man <- n_before - n_man
qc_log[qc_log$Step == "Manual", "Remaining"] <- n_man
qc_log[qc_log$Step == "Manual", "Removed"] <- rem_man


View(data_filtered@meta.data)
#========================================================
#Treatment-wise - Ambient RNA and Doublet filtering
#========================================================
library(DoubletFinder)
library(SoupX)

sample_list <- SplitObject(data_filtered, split.by = "Treatment")

samples <- names(sample_list)

# Count cells after ambiant RNA removal
n_amb <- 0
# Count cells after doublet removal
n_doublets <- 0

# Filtered per sample
doublet_removed <- list()

# Sample-wise preprocessing and doublet detection
for (sample in samples) {
  raw_name <- paste0(sample,"_object")
  norm_name <- paste0(sample,"_normalised")
  
  sample_obj <- sample_list[[sample]]
  og_obj <- sample_obj
  #========================================================
  #Ambient RNA removal using SoupX
  #========================================================
  #For raw counts
  n_amb_rna <- 0
  prefix <- list.files(path = base_path, pattern = paste0(sample,"_matrix"), full.names = FALSE)
  prefix <- sub("_matrix\\.mtx(\\.gz)?$", "", prefix)
  
  raw_counts = ReadMtx(
    mtx = get_file_path(base_path,prefix,"matrix.mtx"),
    features = get_file_path(base_path,prefix,"features.tsv"),
    cells = get_file_path(base_path,prefix,"barcodes.tsv")
  )
  raw_counts<- raw_counts[rownames(sample_obj),]
  
  # Intialising soup Channel
  sc = SoupChannel(raw_counts,sample_obj[["RNA"]]$counts)
  
  tmp_obj <- NormalizeData(sample_obj) #Default = logTransformation
  tmp_obj <- FindVariableFeatures(tmp_obj) #exclude housekeeping genes
  tmp_obj <- ScaleData(tmp_obj)
  tmp_obj <- RunPCA(tmp_obj)
  
    
  #Find significant PCs
  stdv <- tmp_obj[["pca"]]@stdev
  percent_var <- (stdv^2/sum(stdv^2)) * 100
  cumulative_var <- cumsum(percent_var)
  co1 <- which(cumulative_var > 90)[1]
  co2 <- which(diff(percent_var) < 0.1)[1] + 1
  min_pc <- min(co1, co2)
  
  tmp_obj <- RunUMAP(tmp_obj, dims = 1:min_pc)
  tmp_obj <- FindNeighbors(object = tmp_obj, dims = 1:min_pc)              
  tmp_obj <- FindClusters(object = tmp_obj, resolution = 0.1)
  
  sc = setClusters(sc, setNames(tmp_obj$seurat_clusters, colnames(tmp_obj)))
  sc = setContaminationFraction(sc, 0.05) # Automatically estimates contamination rate
  out = adjustCounts(sc) # Produces "cleaned" counts
  
  # Create the actual working object with the cleaned counts
  sample_obj <- CreateSeuratObject(counts = out, project = sample)
  n_amb <- n_amb + (ncol(sample_obj))
  # Transfer mitochondria
  sample_obj[["percent.mt"]] <- og_obj[["percent.mt"]]
  
  rm(og_obj,raw_counts, tmp_obj)
  gc()
  #========================================================
  #Doublet detection using doublet finder
  #========================================================
  
  sample_obj = NormalizeData(sample_obj) #Default = logTransformation
  sample_obj = FindVariableFeatures(sample_obj) #exclude housekeeping genes
  sample_obj = ScaleData(sample_obj)
  sample_obj = RunPCA(sample_obj)
  
  #Find significant PCs
  stdv <- sample_obj[["pca"]]@stdev
  percent_var <- (stdv^2/sum(stdv^2)) * 100
  cumulative_var <- cumsum(percent_var)
  co1 <- which(cumulative_var > 90)[1]
  co2 <- which(diff(percent_var) < 0.1)[1] + 1
  min_pc <- min(co1, co2)
  
  sample_obj <- RunUMAP(sample_obj, dims = 1:min_pc)
  sample_obj <- FindNeighbors(object = sample_obj, dims = 1:min_pc)              
  sample_obj <- FindClusters(object = sample_obj, resolution = 0.1)
  
  #Find the optimal pK value using paramsweep
  sweep_list <- paramSweep(sample_obj, PCs = 1:min_pc, sct = FALSE)   
  sweep_stats <- summarizeSweep(sweep_list)
  bcmvn <- find.pK(sweep_stats) # computes a metric to find the optimal pK value (max mean variance normalised by modality coefficient)
  # Optimal pK is the max of the bimodality coefficient (BCmvn) distribution
  optimal.pk <- bcmvn %>% 
    dplyr::filter(BCmetric == max(BCmetric)) %>%
    dplyr::select(pK)
  optimal.pk <- as.numeric(as.character(optimal.pk[[1]]))
  
  # Homotypic Doublet proportion estimate - for nexp
  annotations <- sample_obj@meta.data$seurat_clusters 
  homotypic.prop <- modelHomotypic(annotations) # get proportions of homotypic doublets
  
  ## Recovered cells ~10000
  multiplet_rate <- 0.076
  nExp.poi <- round(multiplet_rate * nrow(sample_obj@meta.data)) # multiply by number of cells to get the number of expected multiplets
  nExp.poi.adj <- round(nExp.poi * (1 - homotypic.prop)) # expected number of doublets
  
  # Running doublet finder
  sample_obj <- doubletFinder_v3(seu = sample_obj, 
                          PCs = 1:min_pc, 
                          pK = optimal.pk, # the neighborhood size used to compute the number of artificial nearest neighbours
                          nExp = nExp.poi.adj) # number of expected real doublets
  
  # DoubletFinder adds a metadata column starting with "DF.classifications"
  df_column <- grep("DF.classifications", colnames(sample_obj@meta.data), value = TRUE)
  sample_obj <- subset(sample_obj, subset = get(df_column) == "Singlet")
  
  doublet_removed[[sample]] <- sample_obj
  gc()
  }

#Merging doublet removed datasets
data_merged_DR <- merge(doublet_removed[[1]],
                    y = doublet_removed[[2:length(doublet_removed)]],
                    add.cell.ids = names(doublet_removed),
                    project = "MBC_DR"
                    )

rem_amb <- n_man - n_amb
qc_log[qc_log$Step == "SoupX", "Remaining"] <- n_amb
qc_log[qc_log$Step == "SoupX", "Removed"] <- rem_amb

n_droplet <- ncol(data_merged_DR)
rem_droplet <- n_amb - n_droplet
qc_log[qc_log$Step == "Droplet Finder", "Remaining"] <- n_droplet
qc_log[qc_log$Step == "Droplet Finder", "Removed"] <- rem_droplet


#========================================================
#Normalisation
#========================================================
data_filtered <- data_merged_DR
data_filtered = JoinLayers(data_filtered) #multiple treatment
data_filtered = NormalizeData(data_filtered, normalization.method = 'LogNormalize', scale.factor = 10000) #Default = logTransformation
data_filtered = FindVariableFeatures(data_filtered, selection.method = 'vst', nfeatures = 2000) #exclude housekeeping genes
data_filtered = ScaleData(data_filtered, vars.to.regress = c('nCount_RNA','percent.mt')) #equal weight; so that highly expressed genes doesnt dominate

#========================================================
#Clustering
#========================================================

data_filtered = RunPCA(data_filtered) #Default = logTransformation
ElbowPlot(data_filtered)

#without integration
data_filtered = FindNeighbors(data_filtered, dims = 1:30, reduction = "pca", k.param = 20) #include clusters with 30 dims
data_filtered = FindClusters(data_filtered, resolution = 0.5, cluster.name = "unintegrated clusters")

#UMAP - w/o integration
data_filtered = RunUMAP(data_filtered, dims = 1:30, reduction = "pca", reduction.name = "unintegrated.UMAP")
DimPlot(data_filtered, reduction = "unintegrated.UMAP", group.by = "Treatment")
ggsave("results/UMAP_unintegrated_by_treatment_filter.png")
View(data_filtered@meta.data)

DimPlot(data_filtered, reduction = "unintegrated.UMAP")
ggsave("results/UMAP_filter.png")

#========================================================
#Cell annotation using SingleR
#========================================================
install_if_missing("celldex")
ref <- celldex::ImmGenData() #for mice data
ref$label.main

library(SingleR)
mbc_counts <-data_filtered[["RNA"]]$data
res<-SingleR(test = mbc_counts, ref = ref, labels = ref$label.main)

table(res$labels)

data_filtered$labels <- res[match(rownames(data_filtered@meta.data),rownames(res)),'labels']
View(data_filtered@meta.data)

DimPlot(data_filtered, reduction = "unintegrated.UMAP", group.by = 'labels', label = TRUE)

res

# Manual Annotation confirmation
FeaturePlot(data_filtered, 
            features = c("S100a8", "Ly6g", "Cd3e", "Adgre1"), 
            reduction = "unintegrated.UMAP")




#========================================================
#DEGS
#========================================================
data_deg = data_filtered 
DefaultAssay(data_deg) <- "RNA"

Idents(data_deg) <- "treatment"

deg_phenf = FindMarkers(data_deg, group.by = "treatment",
                        ident.1 = "Phenformin", ident.2 = "Control")
deg_polyic = FindMarkers(data_deg, group.by = "treatment",
                         ident.1 = "PolyIC", ident.2 = "Control")
deg_comb = FindMarkers(data_deg, group.by = "treatment",
                       ident.1 = "Combination", ident.2 = "Control")
# Save DEG results to CSV files
write.csv(deg_phenf, file = "deg_phenformin_vs_control.csv")
write.csv(deg_polyic, file = "deg_polyIC_vs_control.csv")
write.csv(deg_comb, file = "deg_combination_vs_control.csv")




