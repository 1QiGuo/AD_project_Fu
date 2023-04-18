#######################################load packages
library(Seurat)
library(qs)

#######################################input
#--object
load("/fs/ess/PCON0022/guoqi/AD/Spatial_6samples(original name:NC)/output/sample6_clusters_n.RData")

#-----------------read and process annotation data
#--------using pathological annotation with no noise
setwd("/fs/ess/PCON0022/guoqi/AD/Spatial_6samples(original name:NC)/output/pathology(withoutnoise)")
spot_2_3_nonoise <- read.csv("./2-3 AT8 annotation.csv")
spot_2_3_nonoise[, 1] <- sub("1", "2_3", spot_2_3_nonoise[, 1])
colnames(spot_2_3_nonoise)<-c("Barcode","AT8")
spot_2_3_specifc <- read.csv("./2-3 AT8 annotation_specific.csv")
spot_2_3_specifc[, 1] <- sub("1", "2_3", spot_2_3_specifc[, 1])
colnames(spot_2_3_specifc)<-colnames(spot_2_3_nonoise)
table(sample6.combined_withnonoise$patientID)
length(intersect(spot_2_3_nonoise$Barcode, colnames(sample6.combined_withnonoise)))#4289 correct
#2_8
spot_2_8 <- read.csv("./2-8 AT8 annotation.csv")
spot_2_8[, 1] <- sub("1", "2_8", spot_2_8[, 1])
colnames(spot_2_8)<-colnames(spot_2_3_nonoise)
spot_2_8_specific <- read.csv("./2-8 AT8 annotation_specific.csv")
spot_2_8_specific[, 1] <- sub("1", "2_8", spot_2_8_specific[, 1])
colnames(spot_2_8_specific)<-colnames(spot_2_3_nonoise)
length(intersect(spot_2_8_specific$Barcode, colnames(sample6.combined_withnonoise)))#3445 correct
#T4857
spot_T4857 <- read.csv("./T4857 AT8 annotation.csv")
spot_T4857[, 1] <- sub("1", "T4857", spot_T4857[, 1])
colnames(spot_T4857)<-colnames(spot_2_3_nonoise)
spot_T4857_specific <- read.csv("T4857 AT8 annotation_specific.csv")
spot_T4857_specific[, 1] <-
  sub("1", "T4857", spot_T4857_specific[, 1])
colnames(spot_T4857_specific)<-colnames(spot_2_3_nonoise)
length(intersect(spot_T4857$Barcode, colnames(sample6.combined_withnonoise)))#4823 correct
spot_annotation <- rbind(spot_2_3_nonoise, spot_2_8, spot_T4857)
#change blank into level4
spot_annotation[, 2][which(spot_annotation[, 2] == "")] <-
  "non_AT8_level4"
spot_specifc <-
  rbind(spot_2_3_specifc, spot_2_8_specific, spot_T4857_specific)
spot_specifc <- spot_specifc[-which(spot_specifc$AT8 == ""), ]#blank spots can be deleted directly in specific files?

#--------------geneset
geneset1 <- c("KCNLP4", "CSMD1", "NRXN1", "DLGAP1", "EPHA6")
# no KCNLP4
geneset1 <- list(geneset1[-1])
geneset2 <- list(
  c(
    "NRG3",
    "DPP10",
    "MEF2C",
    "RBFOX1",
    "TENM2",
    "SYT1",
    "CNTNAP2" ,
    "GRIN2A",
    "GRIN1",
    "NRXN3"
  )
)
geneset3 <- c(
  "KCNLP4",
  "CSMD1",
  "NRXN1",
  "DLGAP1",
  "EPHA6",
  "NRG3",
  "DPP10" ,
  "MEF2C",
  "RBFOX1",
  "TENM2",
  "SYT1",
  "CNTNAP2",
  "GRIN2A",
  "GRIN1",
  "NRXN3"
)
# no KCNLP4
geneset3 <- list(geneset3[-1])

################################################function
#--statistical summary
data_summary <- function(x) {
  m <- mean(x)
  ymin <- m - sd(x)
  ymax <- m + sd(x)
  return(c(y = m, ymin = ymin, ymax = ymax))
}
#--violin plot
violin_pdf <- function(object, metadata, geneset) {
  object$pathology <- "control"
  ct_object <- subset(object, subset = stage == "control")
  ct_object <- AddModuleScore(object = ct_object,
                              features = geneset,
                              name = 'CD_Features')
  violin_data_control <-
    data.frame(
      barcode = colnames(ct_object),
      category = ct_object$pathology,
      modulescore = ct_object$CD_Features1
    )
  #midad
  midad_object <- subset(object, cells = metadata$Barcode)
  midad_object$pathology[match(metadata[, 1], colnames(midad_object))] <-
    metadata[, 2]
  violin_data <- data.frame()
  for (i in unique(midad_object$pathology)) {
    level_object <- subset(midad_object, subset = pathology == i)
    level_object <- AddModuleScore(object = level_object,
                                   features = geneset,
                                   name = 'CD_Features')
    violin_data_temp <-
      data.frame(
        barcode = colnames(level_object),
        category = level_object$pathology,
        modulescore = level_object$CD_Features1
      )
    violin_data <- rbind(violin_data, violin_data_temp)
  }
  violin_data <- rbind(violin_data, violin_data_control)
  category <- sort(unique(midad_object$pathology))
  category <- c(category, "control")
  violin_data$category <-
    factor(violin_data$category, levels = category)
  library(ggpubr)
  library(ggplot2)
  p <-
    ggplot(violin_data, aes(x = category, y = modulescore, fill = category)) +
    geom_violin(trim = FALSE) +
    stat_summary(fun.data = data_summary) +
    #scale_fill_manual(values=c("#56B4E9","#97f7a7","#f58e62"))+
    theme_classic() +
    stat_compare_means(label = "p.signif", method = "anova") +
    theme(text = element_text(size = 15))
  return(p)
}
#--statistical summary
statistical_funtion <-
  function(object, metadata, geneset, outputname) {
    object$pathology <- "control"
    ct_object <- subset(object, subset = stage == "control")
    ct_object <- AddModuleScore(object = ct_object,
                                features = geneset,
                                name = 'CD_Features')
    violin_data_control <-
      data.frame(
        barcode = colnames(ct_object),
        category = ct_object$pathology,
        modulescore = ct_object$CD_Features1
      )
    #midad
    midad_object <- subset(object, cells = metadata$Barcode)
    midad_object$pathology[match(metadata[, 1], colnames(midad_object))] <-
      metadata[, 2]
    violin_data <- data.frame()
    for (i in unique(midad_object$pathology)) {
      level_object <- subset(midad_object, subset = pathology == i)
      level_object <- AddModuleScore(object = level_object,
                                     features = geneset,
                                     name = 'CD_Features')
      violin_data_temp <-
        data.frame(
          barcode = colnames(level_object),
          category = level_object$pathology,
          modulescore = level_object$CD_Features1
        )
      violin_data <- rbind(violin_data, violin_data_temp)
    }
    violin_data <- rbind(violin_data, violin_data_control)
    category <- sort(unique(midad_object$pathology))
    category <- c(category, "control")
    p <-
      as.data.frame(compare_means(modulescore ~ category,  data = violin_data, method = "anova"))
    group_df <- data.frame()
    for (i in unique(category)) {
      group <- violin_data[which(violin_data$category == i), ]
      group_df_temp <- data.frame(
        group = outputname,
        category = i,
        mean = mean(group$modulescore),
        sd = sd(group$modulescore),
        p = p$p,
        p.adj = p$p.adj,
        p.signif = p$p.signif,
        method = p$method
      )
      group_df <- rbind(group_df, group_df_temp)
    }
    return(group_df)
  }
#---pairwise testing
statistical_compare2_funtion <-
  function(object,
           metadata,
           geneset,
           compare1,
           compare2,
           outputname) {
    object$pathology <- "control"
    violin_data <- data.frame()
    if (compare2 == "control") {
      ct_object <- subset(object, subset = stage == "control")
      ct_object <- AddModuleScore(object = ct_object,
                                  features = geneset,
                                  name = 'CD_Features')
      violin_data_control <-
        data.frame(
          barcode = colnames(ct_object),
          category = ct_object$pathology,
          modulescore = ct_object$CD_Features1
        )
      #midad
      metadata <- metadata[which(metadata[, 2] == compare1), ]
      midad_object <- subset(object, cells = metadata$Barcode)
      midad_object$pathology[match(metadata[, 1], colnames(midad_object))] <-
        metadata[, 2]
      midad_object <- AddModuleScore(object = midad_object,
                                     features = geneset,
                                     name = 'CD_Features')
      violin_data <-
        data.frame(
          barcode = colnames(midad_object),
          category = midad_object$pathology,
          modulescore = midad_object$CD_Features1
        )
      violin_data <- rbind(violin_data, violin_data_control)
      category <- c(unique(midad_object$pathology), "control")
    }
    else{
      metadata1 <- metadata[which(metadata[, 2] == compare1), ]
      metadata2 <- metadata[which(metadata[, 2] == compare2), ]
      metadata <- rbind(metadata1, metadata2)
      midad_object <- subset(object, cells = metadata$Barcode)
      midad_object$pathology[match(metadata[, 1], colnames(midad_object))] <-
        metadata[, 2]
      for (a in c(compare1, compare2)) {
        level_object <- subset(midad_object, subset = pathology == a)
        level_object <- AddModuleScore(object = level_object,
                                       features = geneset,
                                       name = 'CD_Features')
        violin_data_temp <-
          data.frame(
            barcode = colnames(level_object),
            category = level_object$pathology,
            modulescore = level_object$CD_Features1
          )
        violin_data <- rbind(violin_data, violin_data_temp)
      }
      category <- sort(unique(midad_object$pathology))
    }
    p <-
      as.data.frame(compare_means(modulescore ~ category,  data = violin_data, method = "wilcox.test"))
    group_df <- data.frame()
    for (i in unique(category)) {
      group <- violin_data[which(violin_data$category == i), ]
      group_df_temp <- data.frame(
        group = outputname,
        category = i,
        mean = mean(group$modulescore),
        sd = sd(group$modulescore),
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

#############################output
#-----------remove level4 spot
#change blank into level4
spot_annotation_new<-spot_annotation[-which(spot_annotation[, 2] == ""), ] 
spot_specifc <-
  rbind(spot_2_3_specifc, spot_2_8_specific, spot_T4857_specific)
spot_specifc_new <- spot_specifc[-which(spot_specifc[, 2] == ""), ] 
violin(sample6.combined_withnonoise,
       spot_annotation_new,
       geneset1,
       "annotation_sixgroups_geneset1_OLIG_nonoise_nolevel4")
violin(sample6.combined_withnonoise,
       spot_annotation_new,
       geneset2,
       "annotation_sixgroups_geneset2_OLIG_nonoise_nolevel4")
violin(sample6.combined_withnonoise,
       spot_annotation_new,
       geneset3,
       "annotation_sixgroups_geneset3_OLIG_nonoise_nolevel4")
violin(sample6.combined_withnonoise,
       spot_specifc_new,
       geneset1,
       "specifc_sixgroups_geneset1_OLIG_nonoise_nolevel4")
violin(sample6.combined_withnonoise,
       spot_specifc_new,
       geneset2,
       "specific_sixgroups_geneset2_OLIG_nonoise_nolevel4")
violin(sample6.combined_withnonoise,
       spot_specifc_new,
       geneset3,
       "specifc_sixgroups_geneset3_OLIG_nonoise_nolevel4")

outputname_all <-
  c(
    "annotation_sixgroups_geneset1_OLIG_nonoise_nolevel4",
    "annotation_sixgroups_geneset2_OLIG_nonoise_nolevel4",
    "annotation_sixgroups_geneset3_OLIG_nonoise_nolevel4",
    "specific_sixgroups_geneset1_OLIG_nonoise_nolevel4",
    "specific_sixgroups_geneset2_OLIG_nonoise_nolevel4",
    "specific_sixgroups_geneset3_OLIG_nonoise_nolevel4"
  )
geneset_all <- c(geneset1, geneset2, geneset3)
geneset_all <- rep(geneset_all, times = 2)
sta_df <- data.frame()
spot_all <- rep(c(list(spot_annotation_new), list(spot_specifc_new)), each = 3)
for (i in 1:6) {
  temp <-
    statistical_funtion(sample6.combined_withnonoise, spot_all[[i]], geneset_all[i], outputname_all[i])
  sta_df <- rbind(sta_df, temp)
}
write.csv(sta_df, "/fs/ess/PCON0022/guoqi/AD/Wenjie/results/statistical_summary/statistical_summary_nonoise_nolevel4_addolig.csv")
#---compare each two groups(wilcoxen test)
sta_df <- data.frame()
for (i in 1:6) {
  for (j in 1:length(unique(spot_all[[i]][, 2]))) {
    conditions <- c(sort(unique(spot_all[[i]][, 2])), "control")
    for (t in (j + 1):length(conditions)) {
      temp <-
        statistical_compare2_funtion(
          sample6.combined_withnonoise,
          spot_all[[i]],
          geneset_all[i],
          conditions[j],
          conditions[t],
          outputname_all[i]
        )
      sta_df <- rbind(sta_df, temp)
    }
  }
}
write.csv(sta_df, "/fs/ess/PCON0022/guoqi/AD/Wenjie/results/statistical_summary/statistical_twocomparison_summary_nonoise_nolevel4_addolig.csv")

#pdf_version
setwd("/fs/ess/PCON0022/guoqi/AD/Wenjie/results/addoligmarker/nolevel4_nonoise_pdf")
# Customizing the output
pdf("annotation_sixgroups_geneset1_OLIG_nonoise_nolevel4_violin.pdf",         # File name
    width = 15, height = 10, # Width and height in inches
    bg = "white" ,         # Background color
    colormodel = "cmyk"  )  # Color model (cmyk is required for most publications)
#paper = "A4")          # Paper size

# Creating a plot
p<-violin_pdf(sample6.combined_withnonoise,
              spot_annotation_new,
              geneset1)
p
# Closing the graphical device
dev.off() 
# Customizing the output
pdf("annotation_sixgroups_geneset2_OLIG_nonoise_nolevel4_violin.pdf",         # File name
    width = 15, height = 10, # Width and height in inches
    bg = "white" ,         # Background color
    colormodel = "cmyk"  )  # Color model (cmyk is required for most publications)
#paper = "A4")          # Paper size
# Creating a plot
g<-violin_pdf(sample6.combined_withnonoise,
              spot_annotation_new,
              geneset2)
g
# Closing the graphical device
dev.off() 

for(i in 1:6){
  pdf(paste0(outputname_all[i],"_violin.pdf"),         # File name
      width = 15, height = 10, # Width and height in inches
      bg = "white" ,         # Background color
      colormodel = "cmyk"  )
  g<-violin_pdf(sample6.combined_withnonoise,
                spot_all[[i]],
                geneset_all[i])
  g
  dev.off() 
}
