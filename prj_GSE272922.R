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
#========================================================
#set paths and Sample info
#========================================================
base_path = "D:/BI_prj/scrna_proj/neutrophils_BC/GSE272922_RAW"
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
  assign(sample,seurat_obj)
}
dim(GSM8414917_Combination)
View(GSM8414918_Control@meta.data)
#========================================================
#Merging datasets
#========================================================
ls()
data_merged = merge(GSM8414918_Control,
             y = list(GSM8414919_Phenformin,GSM8414920_PolyIC,GSM8414917_Combination),
             add.cell.ids = c("Control", "Phenformin", "PolyIC", "Combination"),
             project = "MBC"
)
View(data_merged@meta.data)
data_merged$sample <-rownames(data_merged@meta.data)
data_merged@meta.data<- separate(data_merged@meta.data, col = 'sample', into = c('Treatment','Barcode'),
         sep = '_')

#========================================================
#QC and filtering https://www.sc-best-practices.org/preprocessing_visualization/quality_control.html
#========================================================

#low quality bases - high mit contamination
data_merged[["percent.mt"]] <- PercentageFeatureSet(data_merged, pattern = "^mt-") #mouse is ^mt-

sum(grepl("^mt-", rownames(data_merged)))
head(rownames(data_merged), 20)


QCp1<-VlnPlot(data_merged, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave("results/QC_vlnplot.png", bg = "white", plot = QCp1, width = 8, height = 6, dpi = 300)

QCp2<-FeatureScatter(data_merged, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  geom_smooth(method = 'lm')
ggsave("results/QC_scatterplot.png", bg = "white", plot = QCp2, width = 8, height = 6, dpi = 300)

quantile(data_merged$nFeature_RNA, probs = c(0.01, 0.99))
quantile(data_merged$nCount_RNA, probs = c(0.01, 0.99))
quantile(data_merged$percent.mt, probs = c(0.01, 0.99))


data_filtered <- subset(
  data_merged,
  subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Visualize QC metrics as a violin plot
VlnPlot(data_filtered, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

nrow(data_filtered)           # total number of genes/features
sum(grepl("^mt-", rownames(data_filtered)))  # how many are mitochondrial
Assays(data_filtered)

#========================================================
#Normalisation
#========================================================
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
ggsave("results/UMAP_unintegrated_by_treatment.png")
View(data_filtered@meta.data)

DimPlot(data_filtered, reduction = "unintegrated.UMAP")
ggsave("results/UMAP.png")

#========================================================
#Cell annotation
#========================================================




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




