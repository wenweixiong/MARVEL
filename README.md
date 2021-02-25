# MARVEL
MARVEL is an R package developed for alternative splicing analysis at single-cell resolution. The main functionalities are:
1. Compute PSI values for all five main exon-level splicing events, i.e. skipped-exon (SE),  mutually-exclusive exons (MXE), retained-intron (RI), alternative 5' splice site (A5SS), and alternative 3' splice site (A3SS).
2. Stratify PSI distribution for each splicing event into the five main modalities, i.e. included, excluded, bimodal, middle, and multimodal. Further stratify included and excluded into primary and dispersed sub-modalities. 
3. Perform differential splicing analysis and identify network of genes which are coordinately spliced.
4. Integrate both splicing and gene expression data to compare and contrast splicing and gene expression profiles.

# Installation
```
# Install required Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("AnnotationDbi")
BiocManager::install("GO.db")
BiocManager::install("GOstats")
BiocManager::install("org.Hs.eg.db")

# Install MARVEL package
library(devtools)
install_github("wenweixiong/MARVEL")

# Load package
library(MARVEL)

# Launch vignette
??MARVEL
```

# Tutorial
A comprehensive tutorial for using MARVEL to extract biological insights can be found here: https://wenweixiong.github.io/MARVEL
