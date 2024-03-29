% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Script_DROPLET_05_DE_5_IsoSwitch.R
\name{IsoSwitch.10x}
\alias{IsoSwitch.10x}
\title{Classify gene-splicing relationship}
\usage{
IsoSwitch.10x(
  MarvelObject,
  pval.sj = 0.05,
  log2fc.sj = NULL,
  delta.sj = 5,
  min.gene.norm = 0,
  pval.adj.gene = 0.05,
  log2fc.gene = 0.5
)
}
\arguments{
\item{MarvelObject}{Marvel object. S3 object generated from \code{CompareValues.Genes.10x} function.}

\item{pval.sj}{Numeric value. p-value from differential splicing analysis, below which, the splice junction is considered differentially spliced. Default is \code{0.05}.}

\item{log2fc.sj}{Numeric value. Absolute log2 fold change from differential splicing analysis, above which, the splice junction is considered differentially spliced. This option should be \code{NULL} if \code{delta.sj} has been specified.}

\item{delta.sj}{Numeric value. Absolute difference in average PSI values between the two cell groups, above which, the splice junction is considered differentially spliced. This option should be \code{NULL} if \code{log2fc.sj} has been specified.}

\item{min.gene.norm}{Numeric value. The average normalised gene expression across the two cell groups above which the splice junction is considered differentially spliced. Default is \code{0}.}

\item{pval.adj.gene}{Numeric value. Adjusted p-value from differential gene expression analysis, below which, the gene is considered differentially expressed. Default is \code{0.05}.}

\item{log2fc.gene}{Numeric value. Absolute log2 fold change from differential gene expression analysis, above which, the gene is considered differentially expressed. This option should be \code{NULL} if \code{delta.sj} has been specified.}
}
\value{
An object of class S3 containing new slots \code{MarvelObject$SJ.Gene.Cor$Data}, \code{MarvelObject$SJ.Gene.Cor$Proportion$Plot}, and \code{MarvelObject$SJ.Gene.Cor$Proportion$Table}.
}
\description{
Classify gene-splicing relative changes to each other from cell group 1 to group 2. Classifications are coordinated, opposing, isoform-switching, and complex. In coordinated relationship, both gene and splicing changes in the same direction from cell group 1 to group 2. In opposing relationship, gene changes in the opposite direction relative to splicing from cell group 1 to group 2. In isoform-switching, there is differential splice junction usage without differential expression of the corresponding gene between cell group 1 and group 2. Complex relationship involves genes with both coordinated and opposing relationships with splicing. Only differentially spliced junctions are included for analysis here.
}
\examples{

marvel.demo.10x <- readRDS(system.file("extdata/data",
                               "marvel.demo.10x.rds",
                               package="MARVEL")
                               )

marvel.demo.10x <- readRDS(system.file("extdata/data",
                               "marvel.demo.10x.rds",
                               package="MARVEL")
                               )

marvel.demo.10x <- IsoSwitch.10x(
                        MarvelObject=marvel.demo.10x,
                        pval.sj=0.05,
                        delta.sj=5,
                        min.gene.norm=1.0,
                        pval.adj.gene=0.05,
                        log2fc.gene=0.5
                        )

# Check outputs
marvel.demo.10x$SJ.Gene.Cor$Proportion$Plot
marvel.demo.10x$SJ.Gene.Cor$Proportion$Table
cols <- c("coord.intron", "gene_short_name", "cor.complete")
head(marvel.demo.10x$SJ.Gene.Cor$Data[,cols])
}
