#Integrate six RNA samples
library(qs)
library("Seurat")
set.seed(1234)
dyn.load("/apps/hdf5-serial/gnu/9.1/1.12.0/lib/libhdf5_hl.so.200")
library(hdf5r)
#load data from OSC
setwd("/fs/ess/PCON0022/Yuzhou/Shane_ATAC/Primary_data_processing")
all_file_name <- list.files()
all_file_name<-all_file_name[-c(2,4)]
object <- list()
for(i in 1:6){
  a <-
    Read10X_h5(paste0(all_file_name[i], "/filtered_feature_bc_matrix.h5"))
  x <-
    CreateSeuratObject(counts = a[[1]],project =strsplit(all_file_name[i],"_",)[[1]][1] ,min.cells = 3, min.features = 200)
  x[["percent.mt"]] <- PercentageFeatureSet(x, pattern = "^MT-")
  object[i]<-x
}
names(object)<-all_file_name
setwd("/fs/ess/PCON0022/guoqi/NC/Revision/Output")
qsave(object,"rna6objects_filter_preqc.qs")
#annotation
object[[1]]$orig.ident<-"1-1"
object[[2]]$orig.ident<-"18-64"
object[[3]]$orig.ident<-"2-3"
object[[4]]$orig.ident<-"2-5"
object[[5]]$orig.ident<-"2-8"
object[[6]]$orig.ident<-"T4857"

object[[1]]$stage<-"Control"
object[[2]]$stage<-"Control"
object[[3]]$stage<-"Mid-AD"
object[[4]]$stage<-"Control"
object[[5]]$stage<-"Mid-AD"
object[[6]]$stage<-"Mid-AD"

object[[1]]$condition<-"Control"
object[[2]]$condition<-"Control"
object[[3]]$condition<-"AD"
object[[4]]$condition<-"Control"
object[[5]]$condition<-"AD"
object[[6]]$condition<-"AD"
#QC
object_qc<-list()
VlnPlot(object[[1]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
object_qc[[1]] <- subset(object[[1]], subset = nFeature_RNA > 200 & nFeature_RNA < 10000 & percent.mt < 15)
VlnPlot(object[[2]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
object_qc[[2]] <- subset(object[[2]], subset = nFeature_RNA > 200 & nFeature_RNA < 12000 & percent.mt < 15)
VlnPlot(object[[3]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
object_qc[[3]] <- subset(object[[3]], subset = nFeature_RNA > 200 & nFeature_RNA < 10000 & percent.mt < 15)
VlnPlot(object[[4]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
object_qc[[4]] <- subset(object[[4]], subset = nFeature_RNA > 200 & nFeature_RNA < 10000 & percent.mt < 15)
VlnPlot(object[[5]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
object_qc[[5]] <- subset(object[[5]], subset = nFeature_RNA > 200 & nFeature_RNA < 10000 & percent.mt < 15)
VlnPlot(object[[6]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
object_qc[[6]] <- subset(object[[6]], subset = nFeature_RNA > 200 & nFeature_RNA < 10000 & percent.mt < 15)
#Normalization
object_qc<-lapply(object_qc,function(x){RenameCells(x,add.cell.id = x$orig.ident[1])})
object_sct<-lapply(object_qc,SCTransform)
qsave(object_sct,"rna6_sct.qs")
#Integration
features <- SelectIntegrationFeatures(object.list = object_sct, nfeatures = 3000)
object_sct <- PrepSCTIntegration(object.list = object_sct, anchor.features = features)
anchors <- FindIntegrationAnchors(object.list = object_sct, normalization.method = "SCT",
                                  anchor.features = features)
object.int <- IntegrateData(anchorset = anchors, normalization.method = "SCT")
object.int <- RunPCA(object.int, verbose = FALSE)
object.int <- RunUMAP(object.int, reduction = "pca", dims = 1:30)
DimPlot(object.int, reduction = "umap", group.by = "orig.ident")
qsave(object.int,"rna6_integration.qs")

