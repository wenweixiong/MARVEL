## ----echo=FALSE---------------------------------------------------------------
library(knitr)

## ---- message=FALSE, warning=FALSE--------------------------------------------
# Load MARVEL package
library(MARVEL)

# Load adjunct packages for this tutorial
library(ggplot2)
library(gridExtra)

## ---- eval = FALSE------------------------------------------------------------
#  # Load adjunct packages to support additional functionalities
#  library(AnnotationDbi) # GO analysis
#  library(clusterProfiler)
#  library(org.Hs.eg.db)
#  library(org.Mm.eg.db)

## ---- message=FALSE, warning=FALSE--------------------------------------------
# Load adjunct packages to support additional functionalities
library(plyr) # General data processing
library(ggrepel) # General plotting
library(parallel) # To enable multi-thread during RI PSI computation
library(textclean) # AFE, ALE detection
library(fitdistrplus) # Modality analysis: Fit beta distribution
library(FactoMineR) # PCA: Reduce dimension
library(factoextra) # PCA: Retrieve eigenvalues
library(kSamples) # Anderson-Darling (AD) statistical test
library(twosamples) # D Test Statistic (DTS) statistical test
library(stringr) # Plot GO results

## ---- message=FALSE, warning=FALSE--------------------------------------------
# Load saved MARVEL object
marvel.demo <- readRDS(system.file("extdata/data", "marvel.demo.rds", package="MARVEL"))
class(marvel.demo)

## ---- message=FALSE, warning=FALSE--------------------------------------------
SplicePheno <- marvel.demo$SplicePheno
head(SplicePheno)

## ---- eval = FALSE------------------------------------------------------------
#  # STAR in 1st pass mode
#  STAR --runThreadN 16 \
#       --genomeDir GRCh38_GENCODE_genome_STAR_indexed \
#       --readFilesCommand zcat \
#       --readFilesIn ERR1562083_1_val_1.fq.gz ERR1562083_2_val_2.fq.gz \
#       --outFileNamePrefix SJ/ERR1562083. \
#       --outSAMtype None
#  
#  # STAR in 2nd pass mode
#  STAR --runThreadN 16 \
#       --genomeDir GRCh38_GENCODE_genome_STAR_indexed \
#       --readFilesCommand zcat \
#       --readFilesIn ERR1562083_1_val_1.fq.gz ERR1562083_2_val_2.fq.gz \
#       --outFileNamePrefix ERR1562083. \
#       --sjdbFileChrStartEnd SJ/*SJ.out.tab \
#       --outSAMtype BAM SortedByCoordinate \
#       --outSAMattributes NH HI AS nM XS \
#       --quantMode TranscriptomeSAM

## ---- message=FALSE, warning=FALSE--------------------------------------------
SpliceJunction <- marvel.demo$SpliceJunction
SpliceJunction[1:5,1:5]

## ---- eval = FALSE------------------------------------------------------------
#  rmats \
#      --b1 path_to_BAM_sample_1.txt \
#      --b2 path_to_BAM_sample_2.txt \
#      --gtf gencode.v31.annotation.gtf \
#      --od rMATS/ \
#      --tmp rMATS/ \
#      -t paired \
#      --readLength 125 \
#      --variable-read-length \
#      --nthread 8 \
#      --statoff

## ---- message=FALSE, warning=FALSE--------------------------------------------
SpliceFeature <-marvel.demo$SpliceFeature
lapply(SpliceFeature, head)

## ---- eval = FALSE------------------------------------------------------------
#  bedtools coverage \
#                 -g GRCh38.primary_assembly.genome_bedtools.txt \
#                 -split \
#                 -sorted \
#                 -a RI_Coordinates.bed \
#                 -b ERR1562083.Aligned.sortedByCoord.out.bam > \
#                    ERR1562083.txt \
#                 -d

## ---- message=FALSE, warning=FALSE--------------------------------------------
IntronCounts <- marvel.demo$IntronCounts
IntronCounts[1:5,1:5]

## ---- eval = FALSE------------------------------------------------------------
#  rsem-calculate-expression --bam \
#                            --paired-end \
#                            -p 8 \
#                            ERR1562083.Aligned.toTranscriptome.out.bam \
#                            GRCh38_GENCODE_genome_RSEM_indexed/gencode.v31 \
#                            ERR1562083

## ---- message=FALSE, warning=FALSE--------------------------------------------
Exp <- marvel.demo$Exp
Exp[1:5,1:5]

## ---- message=FALSE, warning=FALSE--------------------------------------------
GeneFeature <- marvel.demo$GeneFeature
head(GeneFeature)

## ---- message=FALSE, warning=FALSE--------------------------------------------
marvel <- CreateMarvelObject(SpliceJunction=SpliceJunction,
                             SplicePheno=SplicePheno,
                             SpliceFeature=SpliceFeature,
                             IntronCounts=IntronCounts,
                             GeneFeature=GeneFeature,
                             Exp=Exp
                             )

## ---- echo=FALSE--------------------------------------------------------------
include_graphics(system.file("extdata/figures", "PSI_Validation.png", package="MARVEL"))

## ---- echo=FALSE--------------------------------------------------------------
include_graphics(system.file("extdata/figures", "PSI_Computation.png", package="MARVEL"))

## ---- message=FALSE, warning=FALSE--------------------------------------------
# Check splicing junction data
marvel.demo <- CheckAlignment(MarvelObject=marvel.demo, level="SJ")

## ---- eval = FALSE------------------------------------------------------------
#  # Validate, filter, compute SE splicing events
#  marvel.demo <- ComputePSI(MarvelObject=marvel.demo,
#                            CoverageThreshold=10,
#                            UnevenCoverageMultiplier=10,
#                            EventType="SE"
#                            )
#  
#  # Validate, filter, compute MXE splicing events
#  marvel.demo <- ComputePSI(MarvelObject=marvel.demo,
#                            CoverageThreshold=10,
#                            UnevenCoverageMultiplier=10,
#                            EventType="MXE"
#                            )
#  
#  # Validate, filter, compute RI splicing events
#  marvel.demo <- ComputePSI(MarvelObject=marvel.demo,
#                            CoverageThreshold=10,
#                            EventType="RI",
#                            thread=4
#                            )
#  
#  # Validate, filter, compute A5SS splicing events
#  marvel.demo <- ComputePSI(MarvelObject=marvel.demo,
#                            CoverageThreshold=10,
#                            EventType="A5SS"
#                            )
#  
#  # Validate, filter, compute A3SS splicing events
#  marvel.demo <- ComputePSI(MarvelObject=marvel.demo,
#                            CoverageThreshold=10,
#                            EventType="A3SS"
#                            )
#  
#  # Validate, filter, compute AFE splicing events
#  marvel.demo <- ComputePSI(MarvelObject=marvel.demo,
#                            CoverageThreshold=10,
#                            EventType="AFE"
#                            )
#  
#  # Validate, filter, compute ALE splicing events
#  marvel.demo <- ComputePSI(MarvelObject=marvel.demo,
#                            CoverageThreshold=10,
#                            EventType="ALE"
#                            )

## ---- message=FALSE, warning=FALSE--------------------------------------------
marvel.demo <- TransformExpValues(MarvelObject=marvel.demo,
                                  offset=1,
                                  transformation="log2",
                                  threshold.lower=1
                                  )

## ---- message=FALSE, warning=FALSE--------------------------------------------
# Check splicing data
marvel.demo <- CheckAlignment(MarvelObject=marvel.demo, level="splicing")

# Check gene data
marvel.demo <- CheckAlignment(MarvelObject=marvel.demo, level="gene")

# Cross-check splicing and gene data
marvel.demo <- CheckAlignment(MarvelObject=marvel.demo, level="splicing and gene")

## ---- message=FALSE, warning=FALSE, fig.width=4, fig.height=4, fig.align="center"----
# Retrieve sample metadata
df.pheno <- marvel.demo$SplicePheno

# Define sample ids
sample.ids <- df.pheno[which(df.pheno$cell.type=="iPSC"), "sample.id"]

# Tabulate expressed events
marvel.demo <- CountEvents(MarvelObject=marvel.demo,
                           sample.ids=sample.ids,
                           min.cells=5
                           )

# Output (1): Plot
marvel.demo$N.Events$Plot

# Output (2): Table
marvel.demo$N.Events$Table

## ---- message=FALSE, warning=FALSE, fig.width=4, fig.height=4, fig.align="center"----
# Retrieve sample metadata
df.pheno <- marvel.demo$SplicePheno

# Define sample ids
sample.ids <- df.pheno[which(df.pheno$cell.type=="Endoderm"), "sample.id"]

# Tabulate expressed events
marvel.demo <- CountEvents(MarvelObject=marvel.demo,
                           sample.ids=sample.ids,
                           min.cells=5
                           )

# Output (1): Plot
marvel.demo$N.Events$Plot

# Output (2): Table
marvel.demo$N.Events$Table

## ---- echo=FALSE--------------------------------------------------------------
include_graphics(system.file("extdata/figures", "Modality.png", package="MARVEL"))

## ---- message=FALSE, warning=FALSE, fig.width=4, fig.height=4, fig.align="center"----
# Retrieve sample metadata
df.pheno <- marvel.demo$SplicePheno

# Define sample IDs
sample.ids <- df.pheno[which(df.pheno$cell.type=="iPSC"), "sample.id"]

# Assign modality
marvel.demo <- AssignModality(MarvelObject=marvel.demo,
                              sample.ids=sample.ids,
                              min.cells=5,
                              seed=1
                              )
                         
marvel.demo$Modality$Results[1:5, c("tran_id", "event_type", "gene_id", "gene_short_name", "modality.bimodal.adj")]

# Tabulate modality proportion (overall)
marvel.demo <- PropModality(MarvelObject=marvel.demo,
                            modality.column="modality.bimodal.adj",
                            modality.type="extended",
                            event.type=c("SE", "MXE", "RI", "A5SS", "A3SS", "AFE", "ALE"),
                            across.event.type=FALSE
                            )

marvel.demo$Modality$Prop$DoughnutChart$Plot
marvel.demo$Modality$Prop$DoughnutChart$Table

## ---- message=FALSE, warning=FALSE, fig.width=6, fig.height=4, fig.align="center"----
# Tabulate modality proportion (by event type)
marvel.demo <- PropModality(MarvelObject=marvel.demo,
                            modality.column="modality.bimodal.adj",
                            modality.type="extended",
                            event.type=c("SE", "MXE", "RI", "A5SS", "A3SS", "AFE", "ALE"),
                            across.event.type=TRUE,
                            prop.test="fisher",
                            prop.adj="fdr",
                            xlabels.size=8
                            )

marvel.demo$Modality$Prop$BarChart$Plot
head(marvel.demo$Modality$Prop$BarChart$Table)

## ---- message=FALSE, warning=FALSE, fig.width=4, fig.height=4, fig.align="center"----
# Retrieve sample metadata
df.pheno <- marvel.demo$SplicePheno

# Define sample IDs
sample.ids <- df.pheno[which(df.pheno$cell.type=="Endoderm"), "sample.id"]

# Assign modality
marvel.demo <- AssignModality(MarvelObject=marvel.demo,
                              sample.ids=sample.ids,
                              min.cells=5,
                              seed=1
                              )
                         
marvel.demo$Modality$Results[1:5, c("tran_id", "event_type", "gene_id", "gene_short_name", "modality.bimodal.adj")]

# Tabulate modality proportion (overall)
marvel.demo <- PropModality(MarvelObject=marvel.demo,
                            modality.column="modality.bimodal.adj",
                            modality.type="extended",
                            event.type=c("SE", "MXE", "RI", "A5SS", "A3SS", "AFE", "ALE"),
                            across.event.type=FALSE
                            )

marvel.demo$Modality$Prop$DoughnutChart$Plot
marvel.demo$Modality$Prop$DoughnutChart$Table

## ---- message=FALSE, warning=FALSE, fig.width=6, fig.height=4, fig.align="center"----
# Tabulate modality proportion (by event type)
marvel.demo <- PropModality(MarvelObject=marvel.demo,
                            modality.column="modality.bimodal.adj",
                            modality.type="extended",
                            event.type=c("SE", "MXE", "RI", "A5SS", "A3SS", "AFE", "ALE"),
                            across.event.type=TRUE,
                            prop.test="fisher",
                            prop.adj="fdr",
                            xlabels.size=8
                           )

marvel.demo$Modality$Prop$BarChart$Plot
head(marvel.demo$Modality$Prop$BarChart$Table)

## ---- echo=FALSE--------------------------------------------------------------
include_graphics(system.file("extdata/figures", "DE.png", package="MARVEL"))

## ---- message=FALSE, warning=FALSE, fig.width=4, fig.height=4, fig.align="center"----
# Define cell groups
    # Retrieve sample metadata
    df.pheno <- marvel.demo$SplicePheno
    
    # Cell group 1 (reference)
    cell.group.g1 <- df.pheno[which(df.pheno$cell.type=="iPSC"), "sample.id"]
    
    # Cell group 2
    cell.group.g2 <- df.pheno[which(df.pheno$cell.type=="Endoderm"), "sample.id"]

# DE analysis
marvel.demo <- CompareValues(MarvelObject=marvel.demo,
                             cell.group.g1=cell.group.g1,
                             cell.group.g2=cell.group.g2,
                             min.cells=3,
                             method="t.test",
                             method.adjust="fdr",
                             level="gene",
                             show.progress=FALSE
                             )

marvel.demo$DE$Exp$Table[1:5, ]

## ---- message=FALSE, warning=FALSE, fig.width=4, fig.height=4, fig.align="center"----
# Plot DE results
marvel.demo <- PlotDEValues(MarvelObject=marvel.demo,
                            pval=0.10,
                            log2fc=0.5,
                            point.size=0.1,
                            level="gene.global",
                            anno=FALSE
                            )

marvel.demo$DE$Exp.Global$Plot
marvel.demo$DE$Exp.Global$Summary
head(marvel.demo$DE$Exp.Global$Table[,c("gene_id", "gene_short_name", "sig")])

## ---- message=FALSE, warning=FALSE, fig.width=4, fig.height=4, fig.align="center"----
# Plot DE results with annotation of selected genes
    # Retrieve DE output table
    results <- marvel.demo$DE$Exp$Table
    
    # Retrieve top genes
    index <- which(results$log2fc > 2 | results$log2fc < -2)
    gene_short_names <- results[index, "gene_short_name"]

    # Plot
    marvel.demo <- PlotDEValues(MarvelObject=marvel.demo,
                                pval=0.10,
                                log2fc=0.5,
                                point.size=0.1,
                                xlabel.size=10,
                                level="gene.global",
                                anno=TRUE,
                                anno.gene_short_name=gene_short_names
                                )

    marvel.demo$DE$Exp.Global$Plot

## ---- message=FALSE, warning=FALSE, fig.width=4, fig.height=4, fig.align="center"----
marvel.demo <- CompareValues(MarvelObject=marvel.demo,
                             cell.group.g1=cell.group.g1,
                             cell.group.g2=cell.group.g2,
                             min.cells=5,
                             method=c("ad", "dts"),
                             method.adjust="fdr",
                             level="splicing",
                             event.type=c("SE", "MXE", "RI", "A5SS", "A3SS", "ALE", "AFE"),
                             show.progress=FALSE
                             )

head(marvel.demo$DE$PSI$Table[["ad"]])
head(marvel.demo$DE$PSI$Table[["dts"]])

## ---- message=FALSE, warning=FALSE, fig.width=4, fig.height=4, fig.align="center"----

marvel.demo <- PlotDEValues(MarvelObject=marvel.demo,
                       method="ad",
                       pval=0.10,
                       level="splicing.distance",
                       anno=TRUE,
                       anno.tran_id=marvel.demo$DE$PSI$Table[["ad"]]$tran_id[c(1:10)]
                       )

marvel.demo$DE$PSI$Plot[["ad"]]

## ---- message=FALSE, warning=FALSE, fig.width=4, fig.height=4, fig.align="center"----
marvel.demo <- CompareValues(MarvelObject=marvel.demo,
                             cell.group.g1=cell.group.g1,
                             cell.group.g2=cell.group.g2,
                             psi.method=c("ad", "dts"),
                             psi.pval=c(0.10, 0.10),
                             psi.delta=0,
                             method.de.gene="t.test",
                             method.adjust.de.gene="fdr",
                             downsample=FALSE,
                             show.progress=FALSE,
                             level="gene.spliced"
                             )

head(marvel.demo$DE$Exp.Spliced$Table)

## ---- message=FALSE, warning=FALSE, fig.width=4, fig.height=4, fig.align="center"----
# Plot: No annotation
marvel.demo <- PlotDEValues(MarvelObject=marvel.demo,
                            method=c("ad", "dts"),
                            psi.pval=c(0.10, 0.10),
                            psi.delta=0,
                            gene.pval=0.10,
                            gene.log2fc=0.5,
                            point.size=0.1,
                            xlabel.size=8,
                            level="gene.spliced",
                            anno=FALSE
                            )

marvel.demo$DE$Exp.Spliced$Plot
marvel.demo$DE$Exp.Spliced$Summary

## ---- message=FALSE, warning=FALSE, fig.width=4, fig.height=4, fig.align="center"----
# Plot: Annotate top genes
results <- marvel.demo$DE$Exp.Spliced$Table

index <- which((results$log2fc > 2 | results$log2fc < -2) & -log10(results$p.val.adj) > 1)
gene_short_names <- results[index, "gene_short_name"]

marvel.demo <- PlotDEValues(MarvelObject=marvel.demo,
                            method=c("ad", "dts"),
                            psi.pval=c(0.10, 0.10),
                            psi.delta=0,
                            gene.pval=0.10,
                            gene.log2fc=0.5,
                            point.size=0.1,
                            xlabel.size=8,
                            level="gene.spliced",
                            anno=TRUE,
                            anno.gene_short_name=gene_short_names
                            )

marvel.demo$DE$Exp.Spliced$Plot

## ---- message=FALSE, warning=FALSE, fig.width=4, fig.height=3, fig.align="center"----
# Define sample groups
    # Retrieve sample metadata
    df.pheno <- marvel.demo$SplicePheno

    # Group 1
    sample.ids.1 <- df.pheno[which(df.pheno$cell.type=="iPSC"), "sample.id"]
    
    # Group 2
    sample.ids.2 <- df.pheno[which(df.pheno$cell.type=="Endoderm"), "sample.id"]

    # Merge
    cell.group.list <- list("iPSC"=sample.ids.1,
                            "Endoderm"=sample.ids.2
                            )

# Retrieve DE genes
  # Retrieve DE result table
  results.de.exp <- marvel.demo$DE$Exp$Table    
  
  # Retrieve relevant gene_ids
  index <- which(results.de.exp$p.val.adj < 0.10 & abs(results.de.exp$log2fc) > 0.5)
  gene_ids <- results.de.exp[index, "gene_id"]

# Reduce dimension
marvel.demo <- RunPCA(MarvelObject=marvel.demo,
                      cell.group.column="cell.type",
                      cell.group.order=c("iPSC", "Endoderm"),
                      cell.group.colors=NULL,
                      min.cells=5,
                      features=gene_ids,
                      point.size=2.5,
                      level="gene"
                      )

marvel.demo$PCA$Exp$Plot

## ---- message=FALSE, warning=FALSE, fig.width=4, fig.height=3, fig.align="center"----
# Retrieve DE tran_ids
method <- c("ad", "dts")

tran_ids.list <- list()
    
for(i in 1:length(method)) {

    results.de.psi <- marvel.demo$DE$PSI$Table[[method[i]]]
    index <- which(results.de.psi$p.val.adj < 0.10 & results.de.psi$outlier==FALSE)
    tran_ids <- results.de.psi[index, "tran_id"]
    tran_ids.list[[i]] <- tran_ids

}

tran_ids <- unique(unlist(tran_ids.list))

# Reduce dimension
marvel.demo <- RunPCA(MarvelObject=marvel.demo,
                      cell.group.column="cell.type",
                      cell.group.order=c("iPSC", "Endoderm"),
                      cell.group.colors=NULL,
                      min.cells=5,
                      features=tran_ids,
                      point.size=2.5,
                      level="splicing",
                      method.impute="random",
                      seed=1
                      )

marvel.demo$PCA$PSI$Plot

## ---- message=FALSE, warning=FALSE, fig.width=4, fig.height=3, fig.align="center"----
# Retrieve relevant gene_ids
results.de.exp <- marvel.demo$DE$Exp$Table
index <- which(results.de.exp$p.val.adj < 0.10 & abs(results.de.exp$log2fc) > 0.5)
gene_ids <- results.de.exp[-index, "gene_id"]

# Reduce dimension
marvel.demo <- RunPCA(MarvelObject=marvel.demo,
                      cell.group.column="cell.type",
                      cell.group.order=c("iPSC", "Endoderm"),
                      cell.group.colors=NULL,
                      min.cells=5,
                      features=gene_ids,
                      point.size=2.5,
                      level="gene"
                      )

marvel.demo$PCA$Exp$Plot

## ---- message=FALSE, warning=FALSE, fig.width=7, fig.height=8, fig.align="center"----
# Retrieve non-DE gene_ids
results.de.exp <- marvel.demo$DE$Exp$Table
index <- which(results.de.exp$p.val.adj > 0.10 )
gene_ids <- results.de.exp[, "gene_id"]

# Retrieve tran_ids
df.feature <- do.call(rbind.data.frame, marvel.demo$SpliceFeatureValidated)
df.feature <- df.feature[which(df.feature$gene_id %in% gene_ids), ]

# Reduce dimension: All DE splicing events
tran_ids <- df.feature$tran_id

marvel.demo <- RunPCA(MarvelObject=marvel.demo,
                      cell.group.column="cell.type",
                      cell.group.order=c("iPSC", "Endoderm"),
                      cell.group.colors=NULL,
                      min.cells=5,
                      features=tran_ids,
                      point.size=2.5,
                      level="splicing",
                      method.impute="random",
                      seed=1
                      )

plot.all <- marvel.demo$PCA$PSI$Plot

# Reduce dimension: SE
tran_ids <- df.feature[which(df.feature$event_type=="SE"), "tran_id"]

marvel.demo <- RunPCA(MarvelObject=marvel.demo,
                      cell.group.column="cell.type",
                      cell.group.order=c("iPSC", "Endoderm"),
                      cell.group.colors=NULL,
                      min.cells=5,
                      features=tran_ids,
                      point.size=2.5,
                      level="splicing",
                      method.impute="random",
                      seed=1
                      )

plot.se <- marvel.demo$PCA$PSI$Plot

# Reduce dimension: MXE
tran_ids <- df.feature[which(df.feature$event_type=="MXE"), "tran_id"]

marvel.demo <- RunPCA(MarvelObject=marvel.demo,
                      cell.group.column="cell.type",
                      cell.group.order=c("iPSC", "Endoderm"),
                      cell.group.colors=NULL,
                      min.cells=5,
                      features=tran_ids,
                      point.size=2.5,
                      level="splicing",
                      method.impute="random",
                      seed=1
                      )

plot.mxe <- marvel.demo$PCA$PSI$Plot

# Reduce dimension: RI
tran_ids <- df.feature[which(df.feature$event_type=="RI"), "tran_id"]

marvel.demo <- RunPCA(MarvelObject=marvel.demo,
                      cell.group.column="cell.type",
                      cell.group.order=c("iPSC", "Endoderm"),
                      cell.group.colors=NULL,
                      min.cells=5,
                      features=tran_ids,
                      point.size=2.5,
                      level="splicing",
                      method.impute="random",
                      seed=1
                      )

plot.ri <- marvel.demo$PCA$PSI$Plot

# Reduce dimension: A5SS
tran_ids <- df.feature[which(df.feature$event_type=="A5SS"), "tran_id"]

marvel.demo <- RunPCA(MarvelObject=marvel.demo,
                      cell.group.column="cell.type",
                      cell.group.order=c("iPSC", "Endoderm"),
                      cell.group.colors=NULL,
                      min.cells=5,
                      features=tran_ids,
                      point.size=2.5,
                      level="splicing",
                      method.impute="random",
                      seed=1
                      )

plot.a5ss <- marvel.demo$PCA$PSI$Plot

# Reduce dimension: A3SS
tran_ids <- df.feature[which(df.feature$event_type=="A3SS"), "tran_id"]

marvel.demo <- RunPCA(MarvelObject=marvel.demo,
                      cell.group.column="cell.type",
                      cell.group.order=c("iPSC", "Endoderm"),
                      cell.group.colors=NULL,
                      min.cells=5,
                      features=tran_ids,
                      point.size=2.5,
                      level="splicing",
                      method.impute="random",
                      seed=1
                      )

plot.a3ss <- marvel.demo$PCA$PSI$Plot

# Reduce dimension: AFE
tran_ids <- df.feature[which(df.feature$event_type=="AFE"), "tran_id"]

marvel.demo <- RunPCA(MarvelObject=marvel.demo,
                      cell.group.column="cell.type",
                      cell.group.order=c("iPSC", "Endoderm"),
                      cell.group.colors=NULL,
                      min.cells=5,
                      features=tran_ids,
                      point.size=2.5,
                      level="splicing",
                      method.impute="random",
                      seed=1
                      )

plot.afe <- marvel.demo$PCA$PSI$Plot

# Reduce dimension: 
tran_ids <- df.feature[which(df.feature$event_type=="ALE"), "tran_id"]

marvel.demo <- RunPCA(MarvelObject=marvel.demo,
                      cell.group.column="cell.type",
                      cell.group.order=c("iPSC", "Endoderm"),
                      cell.group.colors=NULL,
                      min.cells=5,
                      features=tran_ids,
                      point.size=2.5,
                      level="splicing",
                      method.impute="random",
                      seed=1
                      )

plot.ale <- marvel.demo$PCA$PSI$Plot

# Arrange and view plots
# Read plots from right to left for each row
grid.arrange(plot.all, plot.se, 
             plot.mxe, plot.ri, 
             plot.a5ss, plot.a3ss, 
             plot.afe, plot.ale,
             nrow=4)

## ---- message=FALSE, warning=FALSE, fig.width=4, fig.height=3, fig.align="center"----
# Define sample groups
    # Retrieve sample metadata
    df.pheno <- marvel.demo$SplicePheno
    
    # Group 1
    sample.ids.1 <- df.pheno[which(df.pheno$cell.type=="iPSC"), "sample.id"]
    
    # Group 2
    sample.ids.2 <- df.pheno[which(df.pheno$cell.type=="Endoderm"), "sample.id"]

    # Merge
    cell.group.list <- list("iPSC"=sample.ids.1,
                            "Endoderm"=sample.ids.2
                            )

# Assign modality dynamics
marvel.demo <- ModalityChange(MarvelObject=marvel.demo,
                       method=c("ad", "dts"),
                       psi.pval=c(0.10, 0.10)
                       )

marvel.demo$DE$Modality$Plot
head(marvel.demo$DE$Modality$Table)
marvel.demo$DE$Modality$Plot.Stats

## ---- message=FALSE, warning=FALSE, fig.width=8, fig.height=2, fig.align="center"----
# Example 1
tran_id <- "chr4:108620569:108620600|108620656:108620712:+@chr4:108621951:108622024"
  
marvel.demo <- PlotValues(MarvelObject=marvel.demo,
                          cell.group.list=cell.group.list,
                          feature=tran_id,
                          xlabels.size=5,
                          level="splicing",
                          min.cells=5
                          )

plot.1 <- marvel.demo$adhocPlot$PSI

# Example 2
tran_id <- "chr12:110502049:110502117:-@chr12:110499535:110499546:-@chr12:110496012:110496203"

marvel.demo <- PlotValues(MarvelObject=marvel.demo,
                          cell.group.list=cell.group.list,
                          feature=tran_id,
                          xlabels.size=5,
                          level="splicing",
                          min.cells=5
                          )

plot.2 <- marvel.demo$adhocPlot$PSI

# Example 3
tran_id <- "chr9:35685269:35685339:-@chr9:35685064:35685139:-@chr9:35684732:35684807:-@chr9:35684488:35684550"

marvel.demo <- PlotValues(MarvelObject=marvel.demo,
                          cell.group.list=cell.group.list,
                          feature=tran_id,
                          xlabels.size=5,
                          level="splicing",
                          min.cells=5
                          )

plot.3 <- marvel.demo$adhocPlot$PSI

# Example 4
tran_id <- "chr11:85981129:85981228:-@chr11:85978070:85978093:-@chr11:85976623:85976682"

marvel.demo <- PlotValues(MarvelObject=marvel.demo,
                          cell.group.list=cell.group.list,
                          feature=tran_id,
                          xlabels.size=5,
                          level="splicing",
                          min.cells=5
                          )

plot.4 <- marvel.demo$adhocPlot$PSI

# Arrange and view plots
# Read plots from right to left for each row
grid.arrange(plot.1, plot.2, 
             plot.3, plot.4,
             nrow=1)

## ---- message=FALSE, warning=FALSE, fig.width=8, fig.height=2, fig.align="center"----
# Example 1
tran_id <- "chr17:8383254:8382781|8383157:-@chr17:8382143:8382315"
  
marvel.demo <- PlotValues(MarvelObject=marvel.demo,
                          cell.group.list=cell.group.list,
                          feature=tran_id,
                          xlabels.size=5,
                          level="splicing",
                          min.cells=5
                          )

plot.1 <- marvel.demo$adhocPlot$PSI

# Example 2
tran_id <- "chr17:8383157:8383193|8382781:8383164:-@chr17:8382143:8382315"
  
marvel.demo <- PlotValues(MarvelObject=marvel.demo,
                          cell.group.list=cell.group.list,
                          feature=tran_id,
                          xlabels.size=5,
                          level="splicing",
                          min.cells=5
                          )

plot.2 <- marvel.demo$adhocPlot$PSI

# Example 3
tran_id <- "chr15:24962114:24962209:+@chr15:24967029:24967152:+@chr15:24967932:24968082"
  
marvel.demo <- PlotValues(MarvelObject=marvel.demo,
                          cell.group.list=cell.group.list,
                          feature=tran_id,
                          xlabels.size=5,
                          level="splicing",
                          min.cells=5
                          )

plot.3 <- marvel.demo$adhocPlot$PSI

# Example 4
tran_id <- "chr8:144792587:144792245|144792366:-@chr8:144791992:144792140"

marvel.demo <- PlotValues(MarvelObject=marvel.demo,
                          cell.group.list=cell.group.list,
                          feature=tran_id,
                          xlabels.size=5,
                          level="splicing",
                          min.cells=5
                          )

plot.4 <- marvel.demo$adhocPlot$PSI

# Arrange and view plots
# Read plots from right to left for each row
grid.arrange(plot.1, plot.2, 
             plot.3, plot.4,
             nrow=1)

## ---- message=FALSE, warning=FALSE, fig.width=8, fig.height=2, fig.align="center"----
# Example 1
tran_id <- "chr5:150449703:150449739|150449492:150449696:-@chr5:150447585:150447735"
  
marvel.demo <- PlotValues(MarvelObject=marvel.demo,
                          cell.group.list=cell.group.list,
                          feature=tran_id,
                          xlabels.size=5,
                          level="splicing",
                          min.cells=5
                          )

plot.1 <- marvel.demo$adhocPlot$PSI

# Example 2
tran_id <- "chr12:56725340:56724962|56725263:-@chr12:56724452:56724523"

marvel.demo <- PlotValues(MarvelObject=marvel.demo,
                          cell.group.list=cell.group.list,
                          feature=tran_id,
                          xlabels.size=5,
                          level="splicing",
                          min.cells=5
                          )

plot.2 <- marvel.demo$adhocPlot$PSI

# Example 3
tran_id <- "chr10:78037194:78037304:+@chr10:78037439:78037441:+@chr10:78040204:78040225"

marvel.demo <- PlotValues(MarvelObject=marvel.demo,
                          cell.group.list=cell.group.list,
                          feature=tran_id,
                          xlabels.size=5,
                          level="splicing",
                          min.cells=5
                          )

plot.3 <- marvel.demo$adhocPlot$PSI

# Example 4
tran_id <- "chr10:78037194:78037304:+@chr10:78037439|78040204:78040225"

marvel.demo <- PlotValues(MarvelObject=marvel.demo,
                          cell.group.list=cell.group.list,
                          feature=tran_id,
                          xlabels.size=5,
                          level="splicing",
                          min.cells=5
                          )

plot.4 <- marvel.demo$adhocPlot$PSI

# Arrange and view plots
# Read plots from right to left for each row
grid.arrange(plot.1, plot.2, 
             plot.3, plot.4,
             nrow=1)

## ---- message=FALSE, warning=FALSE, fig.width=4, fig.height=3, fig.align="center"----
marvel.demo <- IsoSwitch(MarvelObject=marvel.demo,
                         method=c("ad", "dts"),
                         psi.pval=c(0.10, 0.10),
                         psi.delta=0,
                         gene.pval=0.10,
                         gene.log2fc=0.5
                         )

marvel.demo$DE$Cor$Plot
head(marvel.demo$DE$Cor$Table)
marvel.demo$DE$Cor$Plot.Stats

## ---- message=FALSE, warning=FALSE, fig.width=5, fig.height=5, fig.align="center"----
# Define cell groups
    # Retrieve sample metadata
    df.pheno <- marvel.demo$SplicePheno
    
    # Group 1
    sample.ids.1 <- df.pheno[which(df.pheno$cell.type=="iPSC"), "sample.id"]
    
    # Group 2
    sample.ids.2 <- df.pheno[which(df.pheno$cell.type=="Endoderm"), "sample.id"]

    # Merge
    cell.group.list <- list("iPSC"=sample.ids.1,
                            "Endoderm"=sample.ids.2
                            )
    
# Example 1
  # Gene
  df.feature <- marvel.demo$GeneFeature
  gene_id <- df.feature[which(df.feature$gene_short_name=="CMC2"), "gene_id"]

  marvel.demo <- PlotValues(MarvelObject=marvel.demo,
                            cell.group.list=cell.group.list,
                            feature=gene_id,
                            maintitle="gene_short_name",
                            xlabels.size=7,
                            level="gene"
                            )
  
  plot.1_gene <- marvel.demo$adhocPlot$Exp

  # Splicing
  tran_id <- "chr16:80981806:80981877:-@chr16:80980808:80980879|80976003:80976179"
    
  marvel.demo <- PlotValues(MarvelObject=marvel.demo,
                            cell.group.list=cell.group.list,
                            feature=tran_id,
                            xlabels.size=7,
                            level="splicing",
                            min.cells=5
                            )
  
  plot.1_splicing <- marvel.demo$adhocPlot$PSI
  
# Example 2
  # Gene
  df.feature <- marvel.demo$GeneFeature
  gene_id <- df.feature[which(df.feature$gene_short_name=="HNRNPC"), "gene_id"]

  marvel.demo <- PlotValues(MarvelObject=marvel.demo,
                            cell.group.list=cell.group.list,
                            feature=gene_id,
                            maintitle="gene_short_name",
                            xlabels.size=7,
                            level="gene"
                            )
  
  plot.2_gene <- marvel.demo$adhocPlot$Exp

  # Splicing
  tran_id <- "chr14:21231072:21230958|21230997:-@chr14:21230319:21230366"
    
  marvel.demo <- PlotValues(MarvelObject=marvel.demo,
                            cell.group.list=cell.group.list,
                            feature=tran_id,
                            xlabels.size=7,
                            level="splicing",
                            min.cells=5
                           )
  
  plot.2_splicing <- marvel.demo$adhocPlot$PSI
  
# Arrange and view plots
# Read plots from right to left for each row
grid.arrange(plot.1_gene, plot.1_splicing, 
             plot.2_gene, plot.2_splicing,
             nrow=2)

## ---- message=FALSE, warning=FALSE, fig.width=5, fig.height=5, fig.align="center"----
# Example 1
  # Gene
  df.feature <- marvel.demo$GeneFeature
  gene_id <- df.feature[which(df.feature$gene_short_name=="APOO"), "gene_id"]

  marvel.demo <- PlotValues(MarvelObject=marvel.demo,
                            cell.group.list=cell.group.list,
                            feature=gene_id,
                            maintitle="gene_short_name",
                            xlabels.size=7,
                            level="gene"
                            )
  
  plot.1_gene <- marvel.demo$adhocPlot$Exp

  # Splicing
  tran_id <- "chrX:23840313:23840377:-@chrX:23833353:23833612|23833367:23833510"
    
  marvel.demo <- PlotValues(MarvelObject=marvel.demo,
                            cell.group.list=cell.group.list,
                            feature=tran_id,
                            xlabels.size=7,
                            level="splicing",
                            min.cells=5
                            )
  
  plot.1_splicing <- marvel.demo$adhocPlot$PSI
  
# Example 2
  # Gene
  df.feature <- marvel.demo$GeneFeature
  gene_id <- df.feature[which(df.feature$gene_short_name=="BUB3"), "gene_id"]

  marvel.demo <- PlotValues(MarvelObject=marvel.demo,
                            cell.group.list=cell.group.list,
                            feature=gene_id,
                            maintitle="gene_short_name",
                            xlabels.size=7,
                            level="gene"
                            )
  
  plot.2_gene <- marvel.demo$adhocPlot$Exp

  # Splicing
  tran_id <- "chr10:123162612:123162828:+@chr10:123163820:123170467|123165047:123165365"
    
  marvel.demo <- PlotValues(MarvelObject=marvel.demo,
                            cell.group.list=cell.group.list,
                            feature=tran_id,
                            xlabels.size=7,
                            level="splicing",
                            min.cells=5
                            )
  
  plot.2_splicing <- marvel.demo$adhocPlot$PSI
  
# Arrange and view plots
# Read plots from right to left for each row
grid.arrange(plot.1_gene, plot.1_splicing, 
             plot.2_gene, plot.2_splicing,
             nrow=2)

## ---- message=FALSE, warning=FALSE, fig.width=5, fig.height=5, fig.align="center"----
# Example 1
  # Gene
  df.feature <- marvel.demo$GeneFeature
  gene_id <- df.feature[which(df.feature$gene_short_name=="AC004086.1"), "gene_id"]

  marvel.demo <- PlotValues(MarvelObject=marvel.demo,
                            cell.group.list=cell.group.list,
                            feature=gene_id,
                            maintitle="gene_short_name",
                            xlabels.size=7,
                            level="gene"
                            )
    
  plot.1_gene <- marvel.demo$adhocPlot$Exp

  # Splicing
  tran_id <- "chr12:112409641:112409411|112409587:-@chr12:112408420:112408656"
    
  marvel.demo <- PlotValues(MarvelObject=marvel.demo,
                            cell.group.list=cell.group.list,
                            feature=tran_id,
                            xlabels.size=7,
                            level="splicing",
                            min.cells=5
                            )
      
  plot.1_splicing <- marvel.demo$adhocPlot$PSI
  
# Example 2
  # Gene
  df.feature <- marvel.demo$GeneFeature
  gene_id <- df.feature[which(df.feature$gene_short_name=="ACP1"), "gene_id"]

  marvel.demo <- PlotValues(MarvelObject=marvel.demo,
                            cell.group.list=cell.group.list,
                            feature=gene_id,
                            maintitle="gene_short_name",
                            xlabels.size=7,
                            level="gene"
                            )
  
  plot.2_gene <- marvel.demo$adhocPlot$Exp

  # Splicing
  tran_id <- "chr2:271866:271939:+@chr2:272037:272150:+@chr2:272192:272305:+@chr2:275140:275201"
    
  marvel.demo <- PlotValues(MarvelObject=marvel.demo,
                            cell.group.list=cell.group.list,
                            feature=tran_id,
                            xlabels.size=7,
                            level="splicing",
                            min.cells=5
                            )
  
  plot.2_splicing <- marvel.demo$adhocPlot$PSI
  
# Arrange and view plots
# Read plots from right to left for each row
grid.arrange(plot.1_gene, plot.1_splicing, 
             plot.2_gene, plot.2_splicing,
             nrow=2)

## ---- eval = FALSE------------------------------------------------------------
#  marvel.demo <- BioPathways(MarvelObject=marvel.demo,
#                             method=c("ad", "dts"),
#                             pval=0.10,
#                             species="human"
#                             )
#  
#  head(marvel.demo$DE$BioPathways$Table)

## ---- message=FALSE, warning=FALSE, fig.width=6, fig.height=4, fig.align="center"----
# Plot top pathways
df <- marvel.demo$DE$BioPathways$Table
go.terms <- df$Description[c(1:10)]
              
marvel.demo <- BioPathways.Plot(MarvelObject=marvel.demo,
                                go.terms=go.terms,
                                y.label.size=10
                                )

marvel.demo$DE$BioPathways$Plot

