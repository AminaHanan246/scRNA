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

View(GSM8414918_Control@meta.data)
