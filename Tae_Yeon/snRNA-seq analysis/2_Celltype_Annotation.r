#Integrate six RNA samples
library(qs)
library("Seurat")
set.seed(1234)
library(hdf5r)

object.int<-qread("rna6_integration.qs")
#annotation
#function
annotation<-function(marker,object,assay){
  library(pheatmap)
  marker<-as.data.frame(marker)
  marker_f<-melt(marker,measure.vars = colnames(marker),variable_name = "variable",value.name="marker")
  marker_f<-na.omit(marker_f)
  avg_data<-data.frame(rep(0,length(intersect(rownames(object),marker_f$marker))))
  print(length(intersect(rownames(object),marker_f$marker)))
  for(i in sort(unique(Idents(object)))){
    object_temp<-subset(object,idents=i)
    df<-AverageExpression(object_temp, assays = assay,features = marker_f$marker)
    df<-as.data.frame(df[[1]])
    avg_data<-cbind(avg_data,df)
  }
  avg_data<-avg_data[,-1]
  colnames(avg_data)<-sort(unique(Idents(object)))
  if(length(intersect(rownames(object),marker_f$marker))!=length(marker_f$marker)){
    dep_gene<-setdiff(marker_f$marker,rownames(object))
    for(i in dep_gene){
      marker_f<-marker_f[-which(marker_f$marker==i),]
    }
  }
  marker_f$variable<-factor(marker_f$variable,levels = unique(marker_f$variable))
  sample = data.frame(sample = marker_f$variable)
  color = sample
  levels(color) <- Polychrome::dark.colors(length(unique(marker_f$variable)))
  color <- list(sample = levels(color))
  names(color$sample)<- levels(sample$sample)
  separation_sequence <- cumsum(table(marker_f$variable))
  gaps_row = separation_sequence
  return_list<-list(avg_data,sample,color,separation_sequence)
  names(return_list)<-c("avg_data","sample","color","separation_sequence")
  return(return_list)
}

#cluster
#enhance cluster so that cluster 4 can be divided
object.int <- FindNeighbors(object.int, dims = 1:10)
temp<-FindClusters(object.int,resolution = 0.8)
#FeaturePlot(deepmaps, features =c("VWF","CLDN5","MBP","MOBP") )
#FeaturePlot(temp, features =c("VWF","CLDN5","MBP","MOBP") )
p1 <- DimPlot(temp, group.by = "seurat_clusters",label=T,label.size = 7, repel = TRUE) 
ggsave(
  plot = p1,
  filename = "UMAP_cluster0.8_withlabel.tiff",
  device = "tiff",
  dpi = 150,
  width = 11,
  height = 10,
  units = "in"
)
p2 <- DimPlot(temp, group.by = "seurat_clusters", repel = TRUE)
ggsave(
  plot = p2,
  filename = "UMAP_cluster0.8_nolabel.tiff",
  device = "tiff",
  dpi = 150,
  width = 11,
  height = 10,
  units = "in"
)
#heatmap
DefaultAssay(temp)<-"SCT"
Idents(temp)<-temp$seurat_clusters
#load markers
setwd("/fs/ess/PCON0022/guoqi/Data_raw/Brain_Fu_lab")
marker <- read_excel("New gene list.xlsx", sheet = 1)
marker<-as.data.frame(marker)
marker_nop<-marker[,-8]
test<-annotation(marker_nop,temp,"SCT")
p0<-pheatmap(test[[1]][,c("1","2","11","16","19","21","5","10","17","18","3","4","0",
                          "6","7","8","9","12","14","15","13","20")],
             color = colorRampPalette(c("blue","white","red"))(100),
             cluster_rows = F,
             annotation_row = test[[2]],
             annotation_colors = test[[3]],
             cluster_cols = F,
             scale = "row",border_color = "NA",
             gaps_row = test[[4]],fontsize = 15
)

ggsave(
  plot = p0,
  filename = "Heatmap_6samples_defaultparameter.tiff",
  device = "tiff",
  dpi = 150,
  width = 9,
  height = 12,
  units = "in"
)
setwd("/fs/ess/PCON0022/guoqi/NC/Revision/Output")
qsave(temp,"rna6_integration.qs")
temp<-qread("rna6_integration.qs")
#annotation
temp$celltype<-"unkown"
temp$celltype[which(temp$seurat_clusters==1)]<-"EX"
temp$celltype[which(temp$seurat_clusters==2)]<-"EX"
temp$celltype[which(temp$seurat_clusters==11)]<-"EX"
temp$celltype[which(temp$seurat_clusters==16)]<-"EX"
temp$celltype[which(temp$seurat_clusters==19)]<-"EX"
temp$celltype[which(temp$seurat_clusters==21)]<-"EX"

temp$celltype[which(temp$seurat_clusters==5)]<-"IN"
temp$celltype[which(temp$seurat_clusters==10)]<-"IN"
temp$celltype[which(temp$seurat_clusters==17)]<-"IN"
temp$celltype[which(temp$seurat_clusters==18)]<-"IN"
temp$celltype[which(temp$seurat_clusters==3)]<-"AG"
temp$celltype[which(temp$seurat_clusters==4)]<-"MG"
temp$celltype[which(temp$seurat_clusters==0)]<-"Oligo"
temp$celltype[which(temp$seurat_clusters==6)]<-"Oligo"
temp$celltype[which(temp$seurat_clusters==7)]<-"Oligo"
temp$celltype[which(temp$seurat_clusters==8)]<-"Oligo"
temp$celltype[which(temp$seurat_clusters==9)]<-"Oligo"
temp$celltype[which(temp$seurat_clusters==12)]<-"Oligo"
temp$celltype[which(temp$seurat_clusters==14)]<-"Oligo"
temp$celltype[which(temp$seurat_clusters==15)]<-"Oligo"
temp$celltype[which(temp$seurat_clusters==13)]<-"OPC"
temp$celltype[which(temp$seurat_clusters==20)]<-"Endo"
p3 <- DimPlot(temp, group.by = "celltype",label=T,label.size = 7, repel = TRUE) 
ggsave(
  plot = p3,
  filename = "UMAP_celltype_withlabel.tiff",
  device = "tiff",
  dpi = 150,
  width = 11,
  height = 10,
  units = "in"
)
p4 <- DimPlot(temp, group.by = "orig.ident",label=T,label.size = 7, repel = TRUE) 
ggsave(
  plot = p4,
  filename = "UMAP_sample_withlabel.tiff",
  device = "tiff",
  dpi = 150,
  width = 11,
  height = 10,
  units = "in"
)

umap<-temp@reductions$umap@cell.embeddings
umap<-as.data.frame(umap)
umap<-cbind(umap,celltype=temp$celltype)
umap$`temp$celltype`<-temp$seurat_clusters
colnames(umap)[3]<-"cluster"
write.csv(umap,"seurat_6integrate_coordinate.csv")
umap<-read.csv("seurat_6integrate_coordinate.csv")

#ari
install.packages("MLVSBM")
library("MLVSBM")
dim(temp)
setwd("/fs/ess/PCON0022/guoqi/NC-snrna/atac_output/New_results_qc")
int<-qread("wnn.qs")
dim(int)
Idents(int)<-int$orig.ident
int_6<-subset(int,cells=colnames(temp))
new_cell<-gsub("-","_",colnames(temp))
temp<-RenameCells(temp,new.names = new_cell)
cell_int<-intersect(colnames(temp),colnames(int))
temp_eva<-subset(temp,cells=cell_int)
int_6_eva<-subset(int,cells=cell_int)
ARI(temp_eva$celltype,int_6_eva$celltype)
0.9300079
annotation<-data.frame(barcode=colnames(temp),celltype=temp$celltype)
rownames(annotation)<-NULL
write.csv(annotation,"/fs/ess/PCON0022/guoqi/NC/Revision/Output/scRNA_6_annotation.csv",row.names=FALSE,quote=FALSE)
read.csv("")
#output h5 for deconvolution
write10xCounts(
  '/fs/ess/PCON0022/guoqi/NC/Revision/Output/scRNA_6.h5',
  temp@assays$SCT@counts,
  barcodes = colnames(temp),
  gene.id = rownames(temp),
  gene.symbol = rownames(temp),
  gene.type = "Gene Expression",
  overwrite = FALSE,
  type = c("HDF5"),
  genome = "unknown",
  version = c("2"),
  chemistry = "Single Cell 3' v3",
  original.gem.groups = 1L,
  library.ids = "custom"
)
#annotation
rna6<-temp
# temp$celltype[which(temp$celltype=="AG")]<-"Astrocyte"
# temp$celltype[which(temp$celltype=="Endo")]<-"Endothelial cells&Pericyte"
# temp$celltype[which(temp$celltype=="EX")]<-"Excitatory neurons"
# temp$celltype[which(temp$celltype=="IN")]<-"Inhibitory neurons"
# temp$celltype[which(temp$celltype=="MG")]<-"Microglia"
# temp$celltype[which(temp$celltype=="Oligo")]<-"Oligodendrocytes"
# ARI(temp$celltype,int_6)

#test for deconvolution
test <- Read10X_h5("./Output/scRNA_6.h5")
test <- CreateSeuratObject(counts = test)

write.csv(table(temp$celltype,temp$condition),"snrna_6_ctproportion.csv")
