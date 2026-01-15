# Multiome snRNA–snATAC Integration Pipeline

## Introduction
This pipeline integrates single-nucleus RNA-seq (snRNA-seq) and single-nucleus ATAC-seq (snATAC-seq) data from multiple human brain samples.  
It performs quality control, dataset integration, multimodal (RNA + ATAC) analysis, clustering, and cell-type annotation.

The workflow is organized into four R scripts executed in order.

---

## Input

### snRNA-seq (MTX format)
Each sample directory must contain:
- `barcodes.tsv.gz`
- `features.tsv.gz`
- `matrix.mtx.gz`

### snATAC-seq
- `filtered_feature_bc_matrix.h5`
- `atac_fragments.tsv.gz`

### Reference and metadata
- Human genome: hg38
- Sample metadata: `orig.ident`, `stage`, `condition`
- Marker gene list (Excel file) for cell-type annotation

---

## Output

- Integrated snATAC-seq object (`atac_8_integration.qs`)
- Integrated snRNA-seq object (`rna8_integration.qs`)
- Multimodal WNN-integrated object (`wnn.qs`)
- UMAP plots (clusters, sample ID, condition, stage)
- Heatmap of marker gene expression
- Final cell-type labels stored in `celltype`

---

## Steps

### 1. snATAC-seq Integration
- Perform QC on ATAC data
- Normalize using TF-IDF
- Reduce dimensions with LSI
- Integrate all ATAC datasets using a reference sample
- Generate ATAC UMAP

---

### 2. snRNA-seq Integration
- Load RNA data from MTX files
- Perform QC filtering
- Normalize using SCTransform
- Integrate datasets across samples
- Run PCA and UMAP

---

### 3. Multimodal Integration (RNA + ATAC)
- Match cells shared between RNA and ATAC
- Combine RNA (PCA) and ATAC (LSI) information
- Build weighted nearest neighbor (WNN) graph
- Run WNN UMAP and clustering

---

### 4. Clustering and Cell-Type Annotation
- Tune clustering resolution
- Visualize UMAP by cluster and metadata
- Generate marker gene heatmap
- Assign cell-type labels based on marker expression

---

## Notes

- Sample **18–64** is used as the ATAC reference
- RNA data are read from MTX format, not HDF5
- QC thresholds may vary between samples
- Cell-type annotation is marker-based and manually curated
