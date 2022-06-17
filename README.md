# MARVEL
MARVEL is an R package developed for alternative splicing analysis at single-cell resolution. MARVEL complements published single-cell splicing softwares with the following features:
1. Percent spliced-in (PSI) quantification for all seven main exon-level splicing events, i.e. skipped-exon (SE), mutually-exclusive exons (MXE), retained-intron (RI), alternative 5' and 3' splice sites (A5SS, A3SS), and alternative first and last exons (AFE, ALE).
2. Stratify PSI distribution for each splicing event into the modalities (discrete splicing patterns), and adjust for technical biases during this assignment.
3. Integrated differential splicing and gene expression analysis to reveal gene-splicing dynamics.
4. Dimension reduction analysis.
5. Pathway enrichment analysis.
6. Splicing-associated nonsense-mediated decay (NMD) prediction.
7. Multiple visualisation functions for exploring splicing and gene expression across cell populations.
8. Supports both plate-based (e.g., Smart-seq2) and droplet-based (e.g., 10x Genomics) single-cell RNA-sequencing data analysis. 
9. In principle,  also applicable to bulk RNA-sequencing data analysis.

# Workflow for plate-based splicing analysis
![](inst/extdata/Cover_Figure.png)

# Installation
Please install the following pre-requisite R packages from Bioconductor prior to installing MARVEL.
```
# Install required Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("AnnotationDbi")
BiocManager::install("Biostrings")
BiocManager::install("BSgenome")
BiocManager::install("BSgenome.Hsapiens.NCBI.GRCh38")
BiocManager::install("clusterProfiler")
BiocManager::install("GenomicRanges")
BiocManager::install("GenomicScores")
BiocManager::install("GO.db")
BiocManager::install("IRanges")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("org.Mm.eg.db")
BiocManager::install("phastCons100way.UCSC.hg38")
BiocManager::install("waddR")
```
 
Please install the following pre-requisite R packages from CRAN prior to installing MARVEL.
```
install.packages("factoextra")
install.packages("FactoMineR")
install.packages("fitdistrplus")
install.packages("ggplot2")
install.packages("ggrepel")
install.packages("gtools")
install.packages("kableExtra")
install.packages("knitr")
install.packages("kSamples")
install.packages("kuiper.2samp")
install.packages("markdown")
install.packages("Matrix")
install.packages("pheatmap")
install.packages("plyr")
install.packages("reshape2")
install.packages("rmarkdown")
install.packages("S4Vectors")
install.packages("scales")
install.packages("scran")
install.packages("stringr")
install.packages("textclean")
install.packages("twosamples")
```

Please install the modified wiggleplotr R package from here: http://datashare.molbiol.ox.ac.uk/public/wwen/wiggleplotr_1.18.0_master.tar.gz
```
install.packages("wiggleplotr_1.18.0.tar.gz", repos=NULL, type="source")
```

Finally, please install MARVEL from this repository.
```
# Install MARVEL package
library(devtools)
install_github("wenweixiong/MARVEL")

# Load package
library(MARVEL)
```


# Tutorial
Single-cell plate-based alternative splicing analysis: https://wenweixiong.github.io/MARVEL_Plate.html  
Single-cell droplet-based alternative splicing analysis: https://wenweixiong.github.io/MARVEL_Droplet.html

# Preprint

# Contact
Sean Wen <wei.wen@imm.ox.ac.uk>
