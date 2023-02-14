# Violin plots

```{r}
data_summary <- function(x) {
  m <- mean(x)
  ymin <- m - sd(x)
  ymax <- m + sd(x)
  return(c(y = m, ymin = ymin, ymax = ymax))
}
#--violin plot
violin <- function(object, metadata, geneset, outputname) {
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
  ggsave(
    plot = p,
    filename = paste0(
      outputname,
      "_violin.tiff"
    ),
    device = "tiff",
    dpi = 150,
    width = 15,
    height = 10,
    units = "in"
  )
}
```

# Statistical summary for three comparison using annova
```{r}
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
```

## Pairwise testing
```{r}
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
```

# output
```{r}
# violin plot
DefaultAssay(sample6.combined) <- "SCT"
violin(sample6.combined,
       spot_annotation,
       geneset1,
       "annotation_sixgroups_geneset1")

#------statistical summary (anov)
outputname_all <-
  c(
    "annotation_sixgroups_geneset1",
    "annotation_sixgroups_geneset2",
    "annotation_sixgroups_geneset3",
    "specific_sixgroups_geneset1",
    "specific_sixgroups_geneset2",
    "specific_sixgroups_geneset3"
  )
geneset_all <- c(geneset1, geneset2, geneset3)
geneset_all <- rep(geneset_all, times = 2)
sta_df <- data.frame()
spot_all <- rep(c(list(spot_annotation), list(spot_specifc)), each = 3)
for (i in 1:6) {
  temp <-
    statistical_funtion(sample6.combined, spot_all[[i]], geneset_all[i], outputname_all[i])
  sta_df <- rbind(sta_df, temp)
}
write.csv(sta_df, "../results/statistical_summary.csv")
#---compare each two groups(wilcoxen test)
sta_df <- data.frame()
for (i in 1:6) {
  for (j in 1:length(unique(spot_all[[i]][, 2]))) {
    conditions <- c(sort(unique(spot_all[[i]][, 2])), "control")
    for (t in (j + 1):length(conditions)) {
      temp <-
        statistical_compare2_funtion(
          sample6.combined,
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
write.csv(sta_df, "../results/statistical_twocomparison_summary.csv")

```
