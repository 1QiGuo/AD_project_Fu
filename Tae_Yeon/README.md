# PLCG2 computational Analysis

author: "Qi Guo"

date: "2026-01-14"


## Overview

In the folder "snRNA-seq analysis", we preprocessed, integrated, and annotated single-nuclei RNA data from healthy and Alzheimer's disease samples. 

In the folder "Spatial_transcriptomics_analysis", we investigated **PLCG2 gene expression** across spatial transcriptomics data from matched Alzheimer's disease samples, excluding control spots. It includes integration of pathology annotations, data scaling, violin plot visualization, and statistical testing.



## Hardware requirements
No non-standard hardware is required.  
The pipeline runs on standard CPU-based systems. Large datasets may benefit
from increased memory, but this is not mandatory.

---

## Installation
All dependencies are standard R packages available from CRAN and Bioconductor.
No compilation or special installation steps are required.

Typical installation time on a standard desktop computer is approximately
10–20 minutes, depending on internet speed and existing package installations.

---

## Instructions for use
The analysis is performed by running the provided R scripts sequentially.
Users can apply the pipeline to their own data by replacing the input directories
with their snRNA-seq (MTX format) and snATAC-seq data while keeping the same
data structure and metadata fields.

---

## Reproducibility
All quantitative results reported in the manuscript can be reproduced by
running the provided scripts on the full dataset described in the Methods.

## License
This code is released under the MIT License.
