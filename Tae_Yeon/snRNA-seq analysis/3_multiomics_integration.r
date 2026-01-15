library(ggplot2)
atac_int<-qread("/home/qiguo/AD/output/new_qc/atac_8_integration.qs")
#atac_int<-FindClusters(object = atac_int,resolution = 0.2, algorithm = 3)
integrate<-qread("/home/qiguo/AD/output/new_qc/rna8_integration.qs")
name<-colnames(integrate)
name<-gsub("-","_",name)
integrate <- RenameCells(
  integrate,
  new.names = name)
name<-colnames(atac_int)
name<-gsub("-","_",name)
atac_int <- RenameCells(
  atac_int,
  new.names = name)
#integrate<-FindClusters(object = integrate,resolution = 0.2, algorithm = 3)
DefaultAssay(integrate)<-"integrated"
DefaultAssay(atac_int)<-"peaks"
integrate<-subset(integrate,cells=intersect(colnames(integrate),colnames(atac_int)))
atac_int<-subset(atac_int,cells=intersect(colnames(integrate),colnames(atac_int)))
integrate[["peaks"]]<-atac_int@assays$peaks
DefaultAssay(integrate) <- "peaks"
integrate@reductions$lsi<-atac_int@reductions$integrated_lsi
integrate@reductions$umap.atac<-atac_int@reductions$umap
integrate <- FindMultiModalNeighbors(integrate, reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:50))
integrate <- RunUMAP(integrate, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
integrate<-FindClusters(integrate,resolution = 0.5, graph.name = "wsnn", algorithm = 3)
p3 <- DimPlot(integrate, reduction = "wnn.umap", group.by = "wsnn_res.0.5", label = TRUE, label.size = 6, repel = TRUE) + ggtitle("WNN")
qsave(integrate,"wnn.qs")
ggsave(
  plot = p3,
  filename = "./heatmap_rna.tiff",
  device = "tiff",
  dpi = 150,
  width = 9,
  height = 12,
  units = "in"
)
