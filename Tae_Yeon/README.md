---
title: "README: PLCG2 Spatial Transcriptomics Analysis"
author: "Qi Guo"
date: "2025-06-30"
output: html_document
---

## Overview

This analysis investigates **PLCG2 gene expression** across spatial transcriptomics data from Alzheimer's disease samples, excluding control spots. It includes integration of pathology annotations, data scaling, violin plot visualization, and statistical testing.

---

## Input

- `sample6_clusters.RData`: Seurat object containing spatial transcriptomics data for 6 AD samples.
- AT8 pathology annotation CSVs:
  - `2-3 AT8 annotation.csv`
  - `2-8 AT8 annotation.csv`
  - `T4857 AT8 annotation.csv`
- Abeta pathology annotation CSVs:
  - `2-3 Abeta annotation.csv`
  - `2-8 Abeta annotation.csv`
  - `T4857 Abeta annotation.csv`

---

## Output

- Violin plots of PLCG2 expression:
- Statistical comparison results:
  - Statisticalsummary_pairwisecompare
  - Statisticalsummary_allfour
---

## Steps

1. **Load the Seurat object** (`sample6_clusters.RData`) and AT8 annotation files.
2. **Merge annotation data** across three samples and assign pathology levels.
3. **Filter out control spots** and generate a merged Seurat object with updated `pathology` metadata.
4. **Visualize PLCG2 expression** using `ggplot2::geom_violin()` and `Seurat::VlnPlot()`.
5. **Perform statistical testing**:
   - **ANOVA** across all levels.
   - **Wilcoxon rank-sum test** for all pairwise comparisons.
6. **Export results** to CSV and plot to PDF.

---

## Notes

- Control spots are excluded from the analysis.
- Only non-empty AT8, and Abeta levels (`level1`, `level2`, `level3`, `AT8+`) are retained.
- Ensure that the working directories and file paths are updated based on your system.

---
