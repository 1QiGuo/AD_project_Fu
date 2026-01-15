# PLCG2 computational Analysis

author: "Qi Guo"

date: "2026-01-14"


## Overview

In the folder "snRNA-seq analysis", we preprocessed, integrated, and annotated single-cell multiome (RNA+ATAC) data from healthy and Alzheimer's disease samples. 

In the folder "Spatial_transcriptomics_analysis", we investigated **PLCG2 gene expression** across spatial transcriptomics data from matched Alzheimer's disease samples, excluding control spots. It includes integration of pathology annotations, data scaling, violin plot visualization, and statistical testing.

## Notes on usage and reproducibility
- All dependencies are standard R packages.
- Installation typically takes 10–20 minutes.
- Scripts should be executed sequentially.
- Expected outputs include integrated Seurat objects and UMAP visualizations.
- All manuscript results can be reproduced using the full dataset.
