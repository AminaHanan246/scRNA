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
#=================================
#set paths and Sample info
#==================================
base_path = "D:/bioinfo_prj/GSE272922_RAW"
sample_files = list.files(path = base_path, pattern = "matrix\\.mtx(\\.gz)?$"
, full.names = FALSE)
sample_files = sub("matrix\\.mtx(\\.gz)?$", "", sample_files)

#==================================
#Read Each sample and create seurat object
#==================================

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
  seurat_obj = CreateSeuratObject(counts = counts, project = sample)
  assign(sample,seurat_obj)
}
dim(GSM8414917_Combination)
View(GSM8414918_Control@meta.data)

#=======================================================
#Merging datasets
#======================================================
data_merged = merge(GSM8414918_Control,
             y = list(GSM8414919_Phenformin,GSM8414920_PolyIC,GSM8414917_Combination),
             add.cell.ids = c("Control", "Phenformin", "PolyIC", "Combination"),
             projec = "MergedData"
)
View(data_merged@meta.data)

#==========================================
#Splitting geo id and treat
#===========================================

split_id = do.call(rbind,strsplit(data_merged$orig.ident,"_"))

data_merged$geo_id = split_id[,1]
data_merged$treatment = split_id[,2]

#===========================================
#QC and filtering https://www.sc-best-practices.org/preprocessing_visualization/quality_control.html
#===========================================

#low quality bases - high mit contamination
data_merged[["percent.mt"]] <- PercentageFeatureSet(data_merged, pattern = "^mt-") #mouse is ^mt-

sum(grepl("^mt-", rownames(data_merged)))
head(rownames(data_merged), 20)


VlnPlot(data_merged, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(data_merged, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

quantile(data_merged$nFeature_RNA, probs = c(0.01, 0.99))
quantile(data_merged$nCount_RNA, probs = c(0.01, 0.99))
quantile(data_merged$percent.mt, probs = c(0.01, 0.99))


data_filtered <- subset(
  data_merged,
  subset = nFeature_RNA > 500 & nFeature_RNA < 7500 & percent.mt < 10)

# Visualize QC metrics as a violin plot
VlnPlot(data_filtered, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

nrow(data_filtered)           # total number of genes/features
sum(grepl("^mt-", rownames(data_filtered)))  # how many are mitochondrial

#=================================
#Normalising before integration
#=================================

#without integration 
data_filtered = NormalizeData(data_filtered) #Defaul = logTransformation
data_filtered = FindVariableFeatures(data_filtered) #exclude housekeeping genes

data_filtered = ScaleData(data_filtered) #equal weight; so that highly expressed genes doesnt dominate
data_filtered = RunPCA(data_filtered)  #PCA score

#Clustering
data_filtered = FindNeighbors(data_filtered, dims = 1:30, reduction = "pca") #include clusters with 30 dims
data_filtered = FindClusters(data_filtered, resolution = 0.8, cluster.name = "unintegrated clusters")

#UMAP
data_filtered = RunUMAP(data_filtered, dims = 1:30, reduction = "pca", reduction.name = "unintegrated.UMAP")
DimPlot(data_filtered, reduction = "unintegrated.UMAP", group.by = "seurat_clusters")

#Integration and normalisation

data_filtered = IntegrateLayers(object = data_filtered, method = CCAIntegration, orig.reduction = "pca",new.reduction = "integrated.cca", verbose = FALSE)

data_filtered = FindNeighbors(data_filtered, dims = 1:30, reduction = "integrated.cca")
data_filtered = FindClusters(data_filtered, resolution = 0.8, cluster.name= "integrated Clusters")

data_filtered = RunUMAP(data_filtered, dims = 1:30, reduction = "integrated.cca", reduction.name = "integrated.UMAP")
DimPlot(data_filtered, reduction = "integrated.UMAP", group.by = "seurat_clusters")

DimPlot(data_filtered, reduction = "integrated.UMAP", split.by = "treatment", label = TRUE)
