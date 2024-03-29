% Generated by roxygen2: do not edit by hand
% Please edit documentation in
%   R/Script_DROPLET_09_ADHOC_GENE_1_TabulateExpression_Gene.R
\name{adhocGene.TabulateExpression.Gene.10x}
\alias{adhocGene.TabulateExpression.Gene.10x}
\title{Dotplot of gene expression values for a specified gene}
\usage{
adhocGene.TabulateExpression.Gene.10x(
  MarvelObject,
  cell.group.list,
  gene_short_name,
  log2.transform = TRUE,
  min.pct.cells = 10,
  downsample = FALSE,
  seed = 1
)
}
\arguments{
\item{MarvelObject}{Marvel object. S3 object generated from \code{CheckAlignment.10x} function.}

\item{cell.group.list}{List of character strings. Each element of the list is a vector of cell IDs corresponding to a cell group.}

\item{gene_short_name}{Character string. Gene names whose expression will be plotted.}

\item{log2.transform}{Logical value. If set to \code{TRUE} (default), normalised gene expression values will be off-set by 1 and then log2-transformed prior to plotting.}

\item{min.pct.cells}{Numeric value. Percentage of cell expressing the gene in a cell group, below which, the value be re-coded as missing and appear will be omitted from the plot. A gene is considered to be expressed in a given cell if it has non-zero normalised count.}

\item{downsample}{Logical value. If set to \code{TRUE}, the number of cells in each cell group will be down-sampled so that all cell groups will have the same number of cells. The number of cells to down-sample will be based on the smallest cell group. Default is \code{FALSE}.}

\item{seed}{Numeric value. Random number generator to be fixed for down-sampling.}
}
\value{
An object of class S3 with new slots \code{MarvelObject$adhocGene$Expression$Gene$Table}, \code{MarvelObject$adhocGene$Expression$Gene$Plot}, \code{MarvelObject$adhocGene$cell.group.list}, and \code{MarvelObject$adhocGene$gene_short_name}.
}
\description{
Creates a dotplot of average expression value of a specified gene across different cell groups.
}
\examples{

marvel.demo.10x <- readRDS(system.file("extdata/data",
                               "marvel.demo.10x.rds",
                               package="MARVEL")
                               )

# Define cell groups
    # Retrieve sample metadata
    sample.metadata <- marvel.demo.10x$sample.metadata

    # iPSC
    index <- which(sample.metadata$cell.type=="iPSC")
    cell.ids.1 <- sample.metadata[index, "cell.id"]
    length(cell.ids.1)

    # Cardio day 10
    index <- which(sample.metadata$cell.type=="Cardio day 10")
    cell.ids.2 <- sample.metadata[index, "cell.id"]
    length(cell.ids.2)

    # Save into list
    cell.group.list <- list("iPSC"=cell.ids.1,
                            "Cardio d10"=cell.ids.2
                            )

# Gene expression profiling
marvel.demo.10x <- adhocGene.TabulateExpression.Gene.10x(
                        MarvelObject=marvel.demo.10x,
                        cell.group.list=cell.group.list,
                        gene_short_name="TPM2",
                        min.pct.cells=10,
                        downsample=TRUE
                        )

# Check output
marvel.demo.10x$adhocGene$Expression$Gene$Plot
marvel.demo.10x$adhocGene$Expression$Gene$Table
}
