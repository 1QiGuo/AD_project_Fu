# load packages
library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library("dplyr")
library(BSgenome.Hsapiens.UCSC.hg38)
dyn.load("/apps/hdf5-serial/gnu/9.1/1.12.0/lib/libhdf5_hl.so.200")
library(hdf5r)
library(GenomicRanges)
library(GenomeInfoDb)
library(qs)
library(sctransform)
set.seed(123)
library(ggplot2)


# Load data
#_________We choose 18-64 as the reference sample to ensure that there are common features across the eight datasets since it has the most features. So we can identify as much ATAC features as possible
setwd("/fs/ess/PCON0022/Yuzhou/Shane_ATAC/Primary_data_processing")
counts <- Read10X_h5("18-64_results/filtered_feature_bc_matrix.h5")
fragpath <- "18-64_results/atac_fragments.tsv.gz"

# create a Seurat object containing the RNA adata
multi <- CreateSeuratObject(counts = counts$`Gene Expression`, assay = "RNA")

# Get annotation information
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotation) <- "UCSC"

# Create ATAC assay and add it to the object
  multi[["peaks"]] <- CreateChromatinAssay(
    counts = counts$Peaks,
    sep = c(":", "-"),
    fragments = fragpath,
    annotation = annotation
  )
                            
# QC
DefaultAssay(multi) <- "peaks"
multi <- NucleosomeSignal(multi)
multi <- TSSEnrichment(multi)
p <- VlnPlot(
  object = multi,
  features = c(
    "nCount_RNA",
    "nCount_peaks",
    "TSS.enrichment",
    "nucleosome_signal"
    ),
  ncol = 4,
  pt.size = 0
)
setwd("/fs/ess/PCON0022/guoqi/NC-snrna/atac_output")
 ggsave(
   plot = p,
   filename = paste0(names(object)[i], "-qc-vlnplot.png"),
   device = "png",
   dpi = 150,
   width = 12,
   height = 10,
   units = "in"
 )
 
  multi <- subset(
    x = multi,
    subset = nCount_peaks < 40000 &
      nCount_RNA < 75000 &
      nCount_peaks > 1000 &
      nCount_RNA > 1000 &
      nucleosome_signal < 2.5 &
      TSS.enrichment > 1
  )
  p1 <- VlnPlot(
    object = multi,
    features = c(
      "nCount_RNA",
      "nCount_peaks",
      "TSS.enrichment",
      "nucleosome_signal"
    ),
    ncol = 4,
    pt.size = 0
  )

  ggsave(
    plot = p1,
    filename = "18-64_results-qc_ed_-vlnplot.png",
    device = "png",
    dpi = 150,
    width = 12,
    height = 10,
    units = "in"
  )
  # Normalization and 
  DefaultAssay(multi) <- "RNA"
  multi <- SCTransform(multi)
  multi <- RunPCA(multi)
  DefaultAssay(multi) <- "peaks"
  multi <- RunTFIDF(multi)
  multi <- FindTopFeatures(multi, min.cutoff = 5)
  multi <- RunSVD(multi)
  qsave(multi, "/fs/ess/PCON0022/guoqi/NC-snrna/atac_output/reference-18_64.qs")
  
#_________________Integrate 18-64 with other seven datasets
#-------Preprocess other seven data
setwd("/fs/ess/PCON0022/Yuzhou/Shane_ATAC/Primary_data_processing")
## pre integration
# function
pre_integration <- function(fragpath) {
  fragcounts <- CountFragments(fragments = fragpath)
  atac.cells <- fragcounts[fragcounts$frequency_count > 2000, "CB"]
  # create the fragment object
  atac.frags <-
    CreateFragmentObject(path = fragpath, cells = atac.cells)
  # quantify multiome peaks in the scATAC-seq dataset
  counts <- FeatureMatrix(fragments = atac.frags,
                          features = granges(multi),
                          cells = atac.cells)
  # create object
  atac.assay <- CreateChromatinAssay(counts = counts,
                                     min.features = 1000,
                                     fragments = atac.frags)
  pbmc.atac <-
    CreateSeuratObject(counts = atac.assay, assay = "peaks")
  return(pbmc.atac)
}
# input fragment path, output-target object before filtering and LSI
setwd("/fs/ess/PCON0022/Yuzhou/Shane_ATAC/Primary_data_processing")
target_list<-list()
target_file<-all_file_name[-c(3)]
for(i in 1:length(target_file)){
  fragpath<-paste(target_file[i],"atac_fragments.tsv.gz",sep = "/")
  target_list[[i]]<-pre_integration(fragpath)
}
names(target_list)<-target_file
setwd("/fs/ess/PCON0022/guoqi/NC-snrna/atac_output")
qsave(target_list,"preqc_7samples.qs")#note, first two samples are filtered, so should rerun
#QC
p<-VlnPlot(
  object = target_list[[1]],
  features = c("nCount_peaks", "nFeature_peaks"),
  ncol = 2,
  pt.size = 0
)
target_list[[1]] <- subset(target_list[[1]], nCount_peaks > 1000 & nCount_peaks < 30000)
p<-VlnPlot(
  object = target_list[[2]],
  features = c("nCount_peaks", "nFeature_peaks"),
  ncol = 2,
  pt.size = 0
)
 target_list[[2]] <- subset(target_list[[2]], nCount_peaks > 1000 & nCount_peaks < 30000)

p<-VlnPlot(
  object = target_list[[3]],
  features = c("nCount_peaks", "nFeature_peaks"),
  ncol = 2,
  pt.size = 0
)
target_list[[3]] <- subset(target_list[[3]], nCount_peaks > 1000 & nCount_peaks < 40000)

p<-VlnPlot(
  object = target_list[[4]],
  features = c("nCount_peaks", "nFeature_peaks"),
  ncol = 2,
  pt.size = 0
)
target_list[[4]] <- subset(target_list[[4]], nCount_peaks > 1000 & nCount_peaks < 40000)
p<-VlnPlot(
  object = target_list[[5]],
  features = c("nCount_peaks", "nFeature_peaks"),
  ncol = 2,
  pt.size = 0
)
target_list[[5]] <- subset(target_list[[5]], nCount_peaks > 1000 & nCount_peaks < 30000)

p<-VlnPlot(
  object = target_list[[6]],
  features = c("nCount_peaks", "nFeature_peaks"),
  ncol = 2,
  pt.size = 0
)
target_list[[6]] <- subset(target_list[[6]], nCount_peaks > 1000 & nCount_peaks < 30000)
p<-VlnPlot(
  object = target_list[[7]],
  features = c("nCount_peaks", "nFeature_peaks"),
  ncol = 2,
  pt.size = 0
)
target_list[[7]] <- subset(target_list[[7]], nCount_peaks > 1000 & nCount_peaks < 30000)

#Normalization and dimensional reduction
for(i in 1:7){
  target_list[[i]]<-FindTopFeatures(target_list[[i]], min.cutoff = 10)
  target_list[[i]] <- RunTFIDF(target_list[[i]])
  target_list[[i]] <- RunSVD(target_list[[i]])
}
setwd("/fs/ess/PCON0022/guoqi/NC-snrna/atac_output/New_results_qc/")
qsave(target_list,"qc_7samples.qs")
target_list<-qread("qc_7samples.qs")
#annotation
#id
target_list[[1]]$orig.ident<-"1-1"
target_list[[2]]$orig.ident<-"1-7"
target_list[[3]]$orig.ident<-"2-10"
target_list[[4]]$orig.ident<-"2-3"
target_list[[5]]$orig.ident<-"2-5"
target_list[[6]]$orig.ident<-"2-8"
target_list[[7]]$orig.ident<-"T4857"
#stage
target_list[[1]]$stage<-"Control"
target_list[[2]]$stage<-"Late_AD"
target_list[[3]]$stage<-"Late_AD"
target_list[[4]]$stage<-"Mid-AD"
target_list[[5]]$stage<-"Control"
target_list[[6]]$stage<-"Mid-AD"
target_list[[7]]$stage<-"Mid-AD"
#condition
target_list[[1]]$condition<-"Control"
target_list[[2]]$condition<-"AD"
target_list[[3]]$condition<-"AD"
target_list[[4]]$condition<-"AD"
target_list[[5]]$condition<-"Control"
target_list[[6]]$condition<-"AD"
target_list[[7]]$condition<-"AD"
#18-64
multi<-qread("newreference-18_64.qs")
multi$condition<-"Control"
multi$stage<-"Control"
multi$orig.ident<-"18-64"

#Integration:
pbmc.combined<-merge(multi, c(target_list[[1]],target_list[[2]],
                               target_list[[3]],target_list[[4]],
                               target_list[[5]],target_list[[6]],target_list[[7]]),
                      add.cell.ids=c("18_64","1_1","1_7","2_10","2_3","2_5","2_8","T4857"))
# process the combined dataset
pbmc.combined <- FindTopFeatures(pbmc.combined, min.cutoff = 10)
pbmc.combined <- RunTFIDF(pbmc.combined)
pbmc.combined <- RunSVD(pbmc.combined)
pbmc.combined <- RunUMAP(pbmc.combined, reduction = "lsi", dims = 2:30)
p1 <- DimPlot(pbmc.combined, group.by = "orig.ident")

# find integration anchors
resul_8<-append(multi,target_list)
add.cell.ids=c("18_64","1_1","1_7","2_10","2_3","2_5","2_8","T4857")

for(i in 1:8){
  resul_8[[i]] <- RenameCells(
    resul_8[[i]],
    add.cell.id =add.cell.ids[i]
  )
}

integration.anchors <- FindIntegrationAnchors(
  object.list = resul_8,
  anchor.features = rownames(multi),
  reduction = "rlsi",
  dims = 2:30
)

# integrate LSI embeddings
integrated <- IntegrateEmbeddings(
  anchorset = integration.anchors,
  reductions = pbmc.combined[["lsi"]],
  new.reduction.name = "integrated_lsi",
  dims.to.integrate = 1:30
)
# create a new UMAP using the integrated embeddings
integrated <- RunUMAP(integrated, reduction = "integrated_lsi", dims = 2:30)
p2 <- DimPlot(integrated, group.by = "orig.ident")
qsave(integrated,"atac_8_integration.qs")
