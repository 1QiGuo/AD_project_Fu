#load packages
dyn.load("/apps/hdf5-serial/gnu/9.1/1.12.0/lib/libhdf5_hl.so.200")
library(hdf5r)
library(Seurat)
library(qs)
library(ggplot2)
setwd("/fs/ess/PCON0022/guoqi/NC-snrna/codes")

#load data
setwd("/fs/ess/PCON0022/Yuzhou/Shane_ATAC/Primary_data_processing")
all_file_name <- list.files()
object <- list()
for(i in 1:8){
  a <-
    Read10X_h5(paste0(all_file_name[i], "/filtered_feature_bc_matrix.h5"))
  x <-
    CreateSeuratObject(counts = a[[1]],project =strsplit(all_file_name[i],"_")[[1]][1] ,min.cells = 3, min.features = 200)
  x[["percent.mt"]] <- PercentageFeatureSet(x, pattern = "^MT-")
  object[i]<-x
}
names(object)<-all_file_name
object_qc<-list()
object_qc[[1]] <- subset(object[[1]], subset = nFeature_RNA > 200 & nFeature_RNA < 10000 & percent.mt < 10)
object_qc[[2]] <- subset(object[[2]], subset = nFeature_RNA > 200 & nFeature_RNA < 7000 & percent.mt < 10)
object_qc[[3]] <- subset(object[[3]], subset = nFeature_RNA > 200 & nFeature_RNA < 10000 & percent.mt < 10)
object_qc[[4]] <- subset(object[[4]], subset = nFeature_RNA > 200 & nFeature_RNA < 11000 & percent.mt < 10)
object_qc[[5]] <- subset(object[[5]], subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 10)
object_qc[[6]] <- subset(object[[6]], subset = nFeature_RNA > 200 & nFeature_RNA < 10000 & percent.mt < 10)
object_qc[[7]] <- subset(object[[7]], subset = nFeature_RNA > 200 & nFeature_RNA < 9000 & percent.mt < 10)
object_qc[[8]] <- subset(object[[8]], subset = nFeature_RNA > 200 & nFeature_RNA < 10000 & percent.mt < 10)
object_pre_mer<-merge(object_qc[[1]],list(object_qc[[2]],object_qc[[3]],object_qc[[4]],object_qc[[5]],object_qc[[6]],object_qc[[7]],object_qc[[8]]),
                      add.cell.ids=c("1_1","1_7","18_64","2_10","2_3","2_5","2_8","T4857"))
object_list<-SplitObject(object_pre_mer,split.by = "orig.ident")
object_list <- lapply(X = object_list, FUN = SCTransform)
features <- SelectIntegrationFeatures(object.list = object_list, nfeatures = 3000)
object_list <- PrepSCTIntegration(object.list = object_list, anchor.features = features)
anchors <- FindIntegrationAnchors(object.list = object_list, normalization.method = "SCT",
                                  anchor.features = features)
object.int <- IntegrateData(anchorset = anchors, normalization.method = "SCT")
object.int <- RunPCA(object.int, verbose = FALSE)
object.int <- RunUMAP(object.int, reduction = "pca", dims = 1:30)
object.int <- FindNeighbors(object.int, reduction = "pca", dims = 1:30)
object.int <- FindClusters(object.int, resolution = 0.01)
p1 <- DimPlot(object.int, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(object.int, reduction = "umap", ,group.by = "seurat_clusters",label = TRUE, repel = TRUE)
p1+p2
qsave(object.int,"/fs/ess/PCON0022/guoqi/NC-snrna/atac_output/new_rna_integration.qs")
