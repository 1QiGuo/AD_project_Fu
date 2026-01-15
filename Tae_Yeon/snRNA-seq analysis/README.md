# snRNA analysis Pipeline

author: "Qi Guo"

date: "2026-01-14"

## Introduction
This pipeline integrates single-nucleus RNA-seq (snRNA-seq) data from multiple human brain samples.  
It performs quality control, dataset integration, clustering, and cell-type annotation.

The workflow is organized into four R scripts executed in order.

---

## Input

Raw count data can be found at: https://zenodo.org/records/18252374?preview=1&token=eyJhbGciOiJIUzUxMiIsImlhdCI6MTc2ODQ5MDM1OCwiZXhwIjoxNzk5OTcxMTk5fQ.eyJpZCI6ImFhZjE0MzY1LTA0NWEtNDQ3ZS04MGUzLTk1MWE1YzFmMjgxYiIsImRhdGEiOnt9LCJyYW5kb20iOiIxYTg3MjlhZGQ1OGJjY2QzOGY4NWEzMjdhMDE5NmZmNiJ9.YsGq0rtZhzGu9ha8f_Oe8n4K7SWfrFkRmIDHDj7xuCNNN-QgfdzS4vVTFFPzpIVMad1d4_P9XtJSWYIFZWbYHg

### snRNA-seq (MTX format)
Each sample directory must contain:
- `barcodes.tsv.gz`
- `features.tsv.gz`
- `matrix.mtx.gz`


### Reference and metadata
- Human genome: hg38
- Sample metadata: `orig.ident`, `stage`, `condition`
- Marker gene list (Excel file) for cell-type annotation

---

## Output

- Integrated snRNA-seq object (`rna6_integration.qs`)
- UMAP plots (clusters, sample ID, condition, stage)
- Heatmap of marker gene expression
- Final cell-type labels stored in `celltype`

---

## Steps


### 1. snRNA-seq Integration
- Load RNA data from MTX files
- Perform QC filtering
- Normalize using SCTransform
- Integrate datasets across samples
- Run PCA and UMAP

---

### 2. Clustering and Cell-Type Annotation
- Tune clustering resolution
- Visualize UMAP by cluster and metadata
- Generate marker gene heatmap
- Assign cell-type labels based on marker expression

---

## Notes

- QC thresholds may vary between samples
- Cell-type annotation is marker-based and manually curated
