#------from OSC
###########new qc_0.2
int<-FindClusters(int,resolution = 0.2, graph.name = "wsnn", algorithm = 3)
p3 <- DimPlot(int, reduction = "wnn.umap", group.by = "wsnn_res.0.2", label = TRUE, label.size = 6, repel = TRUE) + ggtitle("WNN")
ggsave(
  plot = p3,
  filename = "/fs/ess/PCON0022/guoqi/NC-snrna/atac_output/New_results_qc/UMAP_wnn_0.2_cluster.tiff",
  device = "tiff",
  dpi = 150,
  width = 11,
  height = 10,
  units = "in"
)
p4 <- DimPlot(int, reduction = "wnn.umap", group.by = "wsnn_res.0.2", repel = TRUE) + ggtitle("WNN")
ggsave(
  plot = p4,
  filename = "/fs/ess/PCON0022/guoqi/NC-snrna/atac_output/New_results_qc/UMAP_wnn_0.2_cluster2.tiff",
  device = "tiff",
  dpi = 150,
  width = 11,
  height = 10,
  units = "in"
)
p2 <- DimPlot(int, reduction = "wnn.umap", group.by = "stage",  repel = TRUE) + ggtitle("WNN")
ggsave(
  plot = p2,
  filename = "/fs/ess/PCON0022/guoqi/NC-snrna/atac_output/New_results_qc/UMAP_wnn_0.2_stage.tiff",
  device = "tiff",
  dpi = 150,
  width = 11,
  height = 10,
  units = "in"
)
p1 <- DimPlot(int, reduction = "wnn.umap", group.by = "orig.ident") + ggtitle("WNN")
ggsave(
  plot = p1,
  filename = "/fs/ess/PCON0022/guoqi/NC-snrna/atac_output/New_results_qc/UMAP_wnn_0.2_id.tiff",
  device = "tiff",
  dpi = 150,
  width = 11,
  height = 10,
  units = "in"
)
p0 <- DimPlot(int, reduction = "wnn.umap", group.by = "condition") + ggtitle("WNN")
ggsave(
  plot = p0,
  filename = "/fs/ess/PCON0022/guoqi/NC-snrna/atac_output/New_results_qc/UMAP_wnn_0.2_condition.tiff",
  device = "tiff",
  dpi = 150,
  width = 11,
  height = 10,
  units = "in"
)

##########heatmap
#newqc_new marker list
setwd("/fs/ess/PCON0022/guoqi/Data_raw/Brain_Fu_lab")
marker <- read_excel("New gene list.xlsx", sheet = 1)
marker<-as.data.frame(marker)
marker_f<-melt(marker,measure.vars = colnames(marker),variable_name = "Celltype",value.name="marker")
marker_f<-na.omit(marker_f)
marker_f<-marker_f[-which(marker_f$Celltype=="Pericyte"),]

avg_data<-data.frame(rep(0,23))
for(i in sort(unique(int$wsnn_res.0.2))){
  object<-subset(int,idents=i)
  df<-AverageExpression(object, assays = "SCT",features = marker_f$value)
  df<-as.data.frame(df$SCT)
  avg_data<-cbind(avg_data,df)
}
avg_data<-avg_data[,-1]
colnames(avg_data)<-sort(unique(int$wsnn_res.0.2))
# marker_f<-marker_f[-which(marker_f$value==c("AMBP")),]
# marker_f<-marker_f[-which(marker_f$value==c("COX4I2")),]
#create pheatmap data
marker_f$Celltype<-factor(marker_f$Celltype,levels = unique(marker_f$Celltype))

sample = data.frame(sample = marker_f$Celltype)
color = sample
#install.packages("Polychrome")
levels(color) <- Polychrome::dark.colors(7)
color <- list(sample = levels(color))
names(color$sample)<- levels(sample$sample)

separation_sequence <- cumsum(table(marker_f$Celltype))
gaps_row = separation_sequence
p <- pheatmap(avg_data[,c("1","8","11","12","14","16","17","7","9","10","13","15","3","6","0","2","4","5","18")],
              color = colorRampPalette(c("blue","white","red"))(100),
              cluster_rows = F,
              annotation_row = sample,
              annotation_colors = color,
              cluster_cols = F,
              scale = "row",border_color = "NA",
              gaps_row = separation_sequence,fontsize = 15
)
ggsave(
  plot = p,
  filename = "/fs/ess/PCON0022/guoqi/NC-snrna/atac_output/New_results_qc/Heatmap_0.2_newmarker.tiff",
  device = "tiff",
  dpi = 150,
  width = 9,
  height = 12,
  units = "in"
)
qsave(int,"/fs/ess/PCON0022/guoqi/NC-snrna/atac_output/New_results_qc/wnn.qs")

###############annotation
int$celltype<-"unkown"
int$celltype[which(int$wsnn_res.0.2==3)]<-"Astrocytes"
int$celltype[which(int$wsnn_res.0.2==18)]<-"Endothelial cells&Pericyte"
int$celltype[which(int$wsnn_res.0.2==16)]<-"Excitatory neurons"
int$celltype[which(int$wsnn_res.0.2==17)]<-"Excitatory neurons"
int$celltype[which(int$wsnn_res.0.2==1)]<-"Excitatory neurons"
int$celltype[which(int$wsnn_res.0.2==8)]<-"Excitatory neurons"
int$celltype[which(int$wsnn_res.0.2==11)]<-"Excitatory neurons"
int$celltype[which(int$wsnn_res.0.2==12)]<-"Excitatory neurons"
int$celltype[which(int$wsnn_res.0.2==14)]<-"Excitatory neurons"
int$celltype[which(int$wsnn_res.0.2==7)]<-"Inhibitory neurons"
int$celltype[which(int$wsnn_res.0.2==9)]<-"Inhibitory neurons"
int$celltype[which(int$wsnn_res.0.2==10)]<-"Inhibitory neurons"
int$celltype[which(int$wsnn_res.0.2==13)]<-"Inhibitory neurons"
int$celltype[which(int$wsnn_res.0.2==15)]<-"Inhibitory neurons"
#int$celltype[which(int$wsnn_res.0.2==18)]<-"Inhibitory neurons"
int$celltype[which(int$wsnn_res.0.2==6)]<-"Microglia"
int$celltype[which(int$wsnn_res.0.2==0)]<-"Oligodendrocytes"
int$celltype[which(int$wsnn_res.0.2==2)]<-"Oligodendrocytes"
int$celltype[which(int$wsnn_res.0.2==4)]<-"Oligodendrocytes"
int$celltype[which(int$wsnn_res.0.2==5)]<-"OPC"
