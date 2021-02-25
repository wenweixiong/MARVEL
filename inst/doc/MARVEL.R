## ----message=FALSE, echo=FALSE------------------------------------------------
# Check plot
output <- system.file("extdata", "Vignette_Fig_0_Software_Comparison.png", package="MARVEL")
knitr::include_graphics(output)

## ----eval=FALSE---------------------------------------------------------------
#  library(devtools)
#  install_github("wenweixiong/MARVEL")

## ----eval=FALSE---------------------------------------------------------------
#  install.packages("MARVEL")

## ----message=FALSE------------------------------------------------------------
# Load package
library(MARVEL)

## ----size="small"-------------------------------------------------------------
# Read splicing files
  # Sample metadata
  path_to_file <- system.file("extdata", "SE_phenoData.txt", package="MARVEL")
  df.pheno <- read.table(path_to_file, sep="\t", header=TRUE, stringsAsFactors=FALSE, na.strings="NA")

  # Subset samples that passed sequencing QC
  df.pheno <- df.pheno[which(df.pheno$qc.seq=="pass"), ]
  df.pheno[1:5,]
  
  # Splice junction file
  path_to_file <- system.file("extdata", "SJ.txt", package="MARVEL")
  sj <- read.table(path_to_file, sep="\t", header=TRUE, stringsAsFactors=FALSE, na.strings="NA")
  sj[1:5,1:5]
  
  # Splicing event metadata
  df.feature.list <- list()
  path_to_file <- system.file("extdata", "SE_featureData.txt", package="MARVEL")
  df.feature.list[[1]] <- read.table(path_to_file, sep="\t", header=TRUE, stringsAsFactors=FALSE, na.strings="NA")
  names(df.feature.list) <- "SE"
  df.feature.list[["SE"]][1:5,]

# Read gene files
  # featureData
  path_to_file <- system.file("extdata", "TPM_featureData.txt", package="MARVEL")
  df.tpm.feature <- read.table(path_to_file, sep="\t", header=TRUE, stringsAsFactors=FALSE, na.strings="NA")
  df.tpm.feature[1:5,]

  # phenoData
  path_to_file <- system.file("extdata", "SE_phenoData.txt", package="MARVEL")
  df.tpm.pheno <- read.table(path_to_file, sep="\t", header=TRUE, stringsAsFactors=FALSE, na.strings="NA")

  # Normalised expression matrix
  path_to_file <- system.file("extdata", "TPM.txt", package="MARVEL")
  df.tpm <- read.table(path_to_file, sep="\t", header=TRUE, stringsAsFactors=FALSE, na.strings="NA")
  
  # Log-transform values
  df.tpm[,-1] <- log2(df.tpm[,-1])
  df.tpm[,-1][df.tpm[,-1] < 1] <- 0

  df.tpm[1:5,1:5]
  
  # Subset samples that passed sequencing QC
  df.tpm.pheno <- df.tpm.pheno[which(df.pheno$qc.seq=="pass"), ]
  df.tpm.pheno[1:5,]
  df.tpm <- df.tpm[, which(names(df.tpm) %in% c("gene_id", df.tpm.pheno$sample.id))]
  
# Create Marvel object
marvel <- CreateMarvelObject(
            SplicePheno=df.pheno,          # Sample metadata
            SpliceJunction=sj,             # Splice junction counts 
            SpliceFeature=df.feature.list, # Splicing event metadata
            GenePheno=df.tpm.pheno,        # Sample metadata
            GeneFeature=df.tpm.feature,    # Gene metadata
            Exp=df.tpm                     # Gene expression matrix
            )

## ----message=FALSE, echo=FALSE------------------------------------------------
# Check plot
output <- system.file("extdata", "Vignette_Fig_1_PSI_Formula-min.png", package="MARVEL")
knitr::include_graphics(output)

## ----size="small"-------------------------------------------------------------
# Validate and filter splicing events, compute PSI
marvel <- ComputePSI.SE(MarvelObject=marvel, CoverageThreshold=10)

# Check validated splicing events
marvel$SpliceFeatureValidated$SE[1:5,]

# Check computed PSI values
marvel$PSI$SE[1:5,1:5]

## ----size="small"-------------------------------------------------------------
# Read sample metadata
path_to_file <- system.file("extdata", "SE_phenoData.txt", package="MARVEL")
df.pheno <- read.table(path_to_file, sep="\t", header=TRUE, stringsAsFactors=FALSE, na.strings="NA")

# Subset samples that passed sequencing QC
df.pheno <- df.pheno[which(df.pheno$qc.seq=="pass"), ]
df.pheno[1:5,]
  
# Read pre-validated splicing event metadata
df.feature.list <- list()
path_to_file <- system.file("extdata", "SE_featureDataValidated.txt", package="MARVEL")
df.feature.list[[1]] <- read.table(path_to_file, sep="\t", header=TRUE, stringsAsFactors=FALSE, na.strings="NA")
names(df.feature.list) <- "SE"
df.feature.list[["SE"]][1:5, ]

# Read PSI file (pre-computed)
df.list <- list()
path_to_file <- system.file("extdata", "SE.txt", package="MARVEL")
df.list[[1]] <- read.table(path_to_file, sep="\t", header=TRUE, stringsAsFactors=FALSE, na.strings="NA")
names(df.list) <- "SE"
df.list[["SE"]][1:5,1:5]

# Create Marvel object
marvel.temp <- CreateMarvelObject(
            SplicePheno=df.pheno,                    # Sample metadata
            SpliceFeatureValidated=df.feature.list,  # Validated splicing event metadata
            PSI=df.list,                             # Pre-computed PSI matrices
            GenePheno=df.tpm.pheno,                  # Sample metadata
            GeneFeature=df.tpm.feature,              # Gene metadata
            Exp=df.tpm                               # Gene expression matrix
            )

## ----message=FALSE, echo=FALSE------------------------------------------------
# Check plot
output <- system.file("extdata", "Vignette_Fig_2_PSI_Modalities-min.png", package="MARVEL")
knitr::include_graphics(output)

## ----message=FALSE, echo=FALSE------------------------------------------------
# Check plot
output <- system.file("extdata", "Vignette_Fig_3_PSI_Modalities_Included-min.png", package="MARVEL")
knitr::include_graphics(output)

## ----message=FALSE, echo=FALSE------------------------------------------------
# Check plot
output <- system.file("extdata", "Vignette_Fig_4_PSI_Modalities_Excluded-min.png", package="MARVEL")
knitr::include_graphics(output)

## ----message=FALSE, echo=FALSE------------------------------------------------
# Check plot
output <- system.file("extdata", "Vignette_Fig_5_PSI_True_vs_False_Bimodal-min.png", package="MARVEL")
knitr::include_graphics(output)

## ----size="small"-------------------------------------------------------------
# Assign modality
marvel <- AssignModality(
            MarvelObject=marvel, # Marvel object
            cell.type="iPSC",    # Cell type to include for analysis
            n.cells=25,          # Min. no. of cells PSI != NA required 
            sigma.sq=0.001,      # Variance value below which sub-modality is primary,
                                                # above which sub-modality is dispersed
            bimodal.adjust=TRUE, # Detect and rectify false bimodals
            seed=1               # Ensures MLE model returns reproducible parameter values
            )

marvel$Modality$Results[1:5,]

## ----size="small", warning=FALSE, message=FALSE, echo=FALSE, eval=FALSE-------
#  # Differential splicing analysis
#  marvel <- CompareValues(
#              MarvelObject=marvel,               # Marvel object
#              cell.types=c("iPSC", "Endoderm"),  # Cell types to analyse
#              n.cells=25,                        # Min. no. of cells PSI != NA required
#              method="ks",                       # "ks"/"wilcox"/"t.test"
#              method.adj="fdr",                  # Adjust for multiple testing as per p.adjust
#              level="splicing"                   # "gene"/"splicing" data to analyse
#              )
#  
#  marvel$DE$PSI[1:5,]
#  
#  # Differential gene expression analysis
#  marvel <- CompareValues(
#              MarvelObject=marvel,              # Marvel object
#              cell.types=c("iPSC", "Endoderm"), # Cell types to include for analysis
#              n.cells=3,                        # Min. no. of cells expression value > 1 required
#              method="wilcox",                  # "wilcox"/"t.test"
#              method.adj="fdr",                 # Adjust for multiple testing as per p.adjust
#              level="gene"                      # "gene"/"splicing" data for DE analysis
#              )
#  
#  marvel$DE$Gene[1:5,]

## ----message=FALSE, echo=FALSE------------------------------------------------
# Check plot
output <- system.file("extdata", "Vignette_Fig_6_Appendix_Single-Cell_Correlation-min.png", package="MARVEL")
knitr::include_graphics(output)

## ----message=FALSE, echo=FALSE------------------------------------------------
# Check plot
output <- system.file("extdata", "Vignette_Fig_7_Appendix_Single-Cell-Bulk_Correlation-min.png", package="MARVEL")
knitr::include_graphics(output)

## ----message=FALSE, echo=FALSE------------------------------------------------
# Check plot
output <- system.file("extdata", "Vignette_Fig_8_Appendix_Genes_vs_mRNA_Counts-min.png", package="MARVEL")
knitr::include_graphics(output)

## ----message=FALSE, echo=FALSE------------------------------------------------
# Check plot
output <- system.file("extdata", "Vignette_Fig_9_Appendix_Bimodal_Features-min.png", package="MARVEL")
knitr::include_graphics(output)

## ----message=FALSE, echo=FALSE------------------------------------------------
# Check plot
output <- system.file("extdata", "Vignette_Fig_10_Appendix_Modality_Proportions-min.png", package="MARVEL")
knitr::include_graphics(output)

