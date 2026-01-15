library("Seurat")
library(SeuratData)
#library(SeuratDisk)
#remotes::install_github("mojaveazure/seurat-disk") not work
#devtools::install_github("cellgeni/sceasy")
setwd("/fs/ess/PCON0022/guoqi/AD/Spatial_6samples(original name:NC)/output")
load(
  "/fs/ess/PCON0022/guoqi/AD/Spatial_6samples(original name:NC)/output/sample6_clusters.RData"
)
sample6.combined_withnonoise
table(sample6.combined_withnonoise$patientID)
table(sample6.combined_withnonoise$stage)

#load pathology information
setwd(
  "/fs/ess/PCON0022/guoqi/AD/Spatial_6samples(original name:NC)/output/pathology(withoutnoise)"
)
spot_2_3_nonoise <- read.csv("./2-3 Aβ annotation.csv")
spot_2_3_nonoise[, 1] <- sub("1", "2_3", spot_2_3_nonoise[, 1])
colnames(spot_2_3_nonoise) <- c("Barcode", "AB")
table(sample6.combined_withnonoise$patientID)
length(intersect(
  spot_2_3_nonoise$Barcode,
  colnames(sample6.combined_withnonoise)
))#4289 correct
#2_8
spot_2_8 <- read.csv("./2-8 Aβ annotation.csv")
spot_2_8[, 1] <- sub("1", "2_8", spot_2_8[, 1])
colnames(spot_2_8) <- colnames(spot_2_3_nonoise)
#T4857
spot_T4857 <- read.csv("./T4857 Aβ annotation.csv")
spot_T4857[, 1] <- sub("1", "T4857", spot_T4857[, 1])
colnames(spot_T4857) <- colnames(spot_2_3_nonoise)
spot_annotation <- rbind(spot_2_3_nonoise, spot_2_8, spot_T4857)
#change blank into level4
spot_annotation[, 2][which(spot_annotation[, 2] == "")] <-
  "non_Aβ_level4"
spot_annotation_new <-
  spot_annotation[-which(spot_annotation[, 2] == "non_Aβ_level4"),]


#function
violin <- function(object, metadata) {
  object$pathology <- "control"
  ct_object <- subset(object, subset = stage == "control")
  #midad
  midad_object <- subset(object, cells = metadata$Barcode)
  midad_object$pathology[match(metadata[, 1], colnames(midad_object))] <-metadata[, 2]
  category <- sort(unique(midad_object$pathology))
  category <- c(category, "control")
  object<-merge(ct_object,midad_object)
  object$pathology <-
    factor(object$pathology, levels = category)
  library(ggpubr)
  library(ggplot2)
  # ggsave(
  #   plot = p,
  #   filename = paste0(
  #     "/fs/ess/PCON0022/guoqi/AD/Wenjie/results/",
  #     outputname,
  #     "_violin.tiff"
  #   ),
  #   device = "tiff",
  #   dpi = 150,
  #   width = 15,
  #   height = 10,
  #   units = "in"
  # )
  return(object)
}
sample6.combined_withnonoise_nolevel4<-violin(sample6.combined_withnonoise,spot_annotation_new)
Idents(sample6.combined_withnonoise_nolevel4)<-sample6.combined_withnonoise_nolevel4$pathology
sample6.combined_withnonoise_nolevel4_nocontrol<-subset(sample6.combined_withnonoise_nolevel4,idents = "control", invert = TRUE)

DefaultAssay(sample6.combined_withnonoise_nolevel4_nocontrol)<-"Spatial"
sample6.combined_withnonoise_nolevel4_nocontrol_scaled<-ScaleData(sample6.combined_withnonoise_nolevel4_nocontrol)
# VlnPlot(sample6.combined_withnonoise_nolevel4,assay = "Spatial",slot = "counts",features = c("PLCG2"))
#VlnPlot(sample6.combined_withnonoise_nolevel4,assay = "Spatial",slot="data",features = c("PLCG2"))
VlnPlot(sample6.combined_withnonoise_nolevel4,assay = "SCT",slot = 'data',features = c("PLCG2"),pt.size = 0)
# VlnPlot(sample6.combined_withnonoise_nolevel4,assay = "SCT",features = c("PLCG2"), slot = 'scale.data',pt.size = 0)

violin_data<-data.frame(category=sample6.combined_withnonoise_nolevel4_nocontrol$pathology,
                        expression_PLCG2=sample6.combined_withnonoise_nolevel4_nocontrol@assays$SCT@data["PLCG2",])
#DefaultAssay(sample6.combined_withnonoise_nolevel4_nocontrol)<-"SCT"
# violin_data_fetch<-data.frame(category=sample6.combined_withnonoise_nolevel4$pathology,
#                               expression_PLCG2=FetchData(sample6.combined_withnonoise_nolevel4, vars = "PLCG2", slot = "data")[,1])

#identical(violin_data_fetch$expression_PLCG2,violin_data$expression_PLCG2)

violin_data$category <- gsub("β", "beta", violin_data$category)
violin_data$category<-factor(violin_data$category,levels = c("Abeta","non_Abeta_level1","non_Abeta_level2","non_Abeta_level3"))
p <-
  ggplot(violin_data, aes(x = category, y = expression_PLCG2, fill = category)) +
  geom_violin(trim = T, scale = "width") +
  # stat_summary(fun.data = data_summary) +
  #scale_fill_manual(values=c("#56B4E9","#97f7a7","#f58e62"))+
  theme_classic() +
  #stat_compare_means(label = "p.signif", method = "anova") +
  theme(text = element_text(size = 20),axis.text.x = element_text(angle = 45, vjust = 0.1, hjust=0.05))
#pdf_version
setwd("/fs/ess/PCON0022/guoqi/AD/Tae_yeon/Results")
# Customizing the output
pdf(
  "violin_PLCG2_nonoise_nolevel4_nocontrol_ab.pdf",
  # File name
  width = 10,
  height = 9,
  # Width and height in inches
  bg = "white" ,
  # Background color
  colormodel = "cmyk"
)  # Color model (cmyk is required for most publications)
#paper = "A4")          # Paper size
p
# Closing the graphical device
dev.off()

#statistical summary
statistical_compare2_funtion <-
  function(violin_data,
           gene_index,
           compare1,
           compare2,
           outputname) {
    violin_data_compare1<-violin_data[which(violin_data$category==compare1),c(1,gene_index)]
    violin_data_compare2<-violin_data[which(violin_data$category==compare2),c(1,gene_index)]
    violin_data_compare12<-rbind(violin_data_compare1,violin_data_compare2)
    colnames(violin_data_compare12)[2]<-"expression"
    p <-
      as.data.frame(compare_means(expression ~ category,  data = violin_data_compare12, method = "wilcox.test"))
    group_df <- data.frame()
    category <- c(compare1,compare2)
    for (i in unique(category)) {
      group <- violin_data_compare12[which(violin_data_compare12$category == i), ]
      group_df_temp <- data.frame(
        group = outputname,
        category = i,
        mean = mean(group$expression),
        sd = sd(group$expression),
        p = p$p,
        p.adj = p$p.adj,
        p.signif = p$p.signif,
        method = p$method
      )
      group_df <- rbind(group_df, group_df_temp)
    }
    group_df$comparison <- paste(compare1, "vs", compare2, sep = " ")
    return(group_df)
  }

#---compare each two groups(wilcoxen test)
colnames(violin_data)[2]<-"geneexp"
sta_df <- data.frame()
outputname<-c("PLCG2_pairwise_comparison")
for (i in 1) {
  conditions <- levels(violin_data$category)
  for (j in 1:(length(conditions)-1)) {
    for (t in (j + 1):length(conditions)) {
      temp <-
        statistical_compare2_funtion(
          violin_data,
          i+1,
          conditions[j],
          conditions[t],
          outputname[i]
        )
      sta_df <- rbind(sta_df, temp)
    }
  }
}
write.csv(sta_df, "/fs/ess/PCON0022/guoqi/AD/Tae_yeon/Results/07012025/Statisticalsummary_pairwisecompare_ab_PLCG2_nocontrol.csv")

#five
statistical_funtion <-
  function(violin_data,outputname) {
    p <-
      as.data.frame(compare_means(geneexp ~ category,  data = violin_data, method = "anova"))
    group_df <- data.frame()
    category<-violin_data$category
    for (i in levels(category)) {
      group <- violin_data[which(violin_data$category == i), ]
      group_df_temp <- data.frame(
        group = outputname,
        category = i,
        mean = mean(group$geneexp),
        sd = sd(group$geneexp),
        p = p$p,
        p.adj = p$p.adj,
        p.signif = p$p.signif,
        method = p$method
      )
      group_df <- rbind(group_df, group_df_temp)
    }
    return(group_df)
  }
colnames(violin_data)[2]<-"geneexp"
plcg2_stat_all<-statistical_funtion(violin_data,"PLCG2_allcomparison")
write.csv(plcg2_stat_all,"/fs/ess/PCON0022/guoqi/AD/Tae_yeon/Results/07012025/Statisticalsummary_allfourlevelscompare_ab_PLCG2_nocontrol_ab.csv")
