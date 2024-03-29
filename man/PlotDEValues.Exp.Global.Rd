% Generated by roxygen2: do not edit by hand
% Please edit documentation in
%   R/Script_PLATE_09_DE_2_4_PlotDEValues_Gene_Mean_Global.R
\name{PlotDEValues.Exp.Global}
\alias{PlotDEValues.Exp.Global}
\title{Plot global differential gene expression analysis results}
\usage{
PlotDEValues.Exp.Global(
  MarvelObject,
  pval = 0.1,
  log2fc = 0.5,
  point.size = 1,
  anno = FALSE,
  anno.gene_short_name = NULL,
  label.size = 2.5,
  y.upper.offset = 5,
  xlabel.size = 8
)
}
\arguments{
\item{MarvelObject}{Marvel object. S3 object generated from \code{CompareValues} function.}

\item{pval}{Numeric value. Adjusted p-value below which the genes are considered as statistically significant and will consequently be color-annotated on the plot.}

\item{log2fc}{Numeric value. The positive (and negative) value specified above (and below) which the genes are considered to be statistically significant and will consequently be color-annotated on the plot.}

\item{point.size}{Numeric value. The point size for the data points. Default value is \code{1}.}

\item{anno}{Logical value. If set to \code{TRUE}, the specific gene names will be annotated on the plot as defined in \code{anno.gene_short_name} option.}

\item{anno.gene_short_name}{Vector of character strings. When \code{anno} set to \code{TRUE}, the gene names to be annotated on the plot.}

\item{label.size}{Numeric value. Only applicable if \code{anno} set to TRUE. Size of the gene name labels.}

\item{y.upper.offset}{Numeric value. The value in -log10(p-value) to increase the upper limit of the y-axis. To be used when \code{anno} set to TRUE so that gene labels will not be truncated at the upper limit of the y-axis.}

\item{xlabel.size}{Numeric value. Font size of the xtick labels. Default is \code{8}.}
}
\value{
An object of class S3 with new slots \code{MarvelObject$DE$Exp.Global$Table}, \code{MarvelObject$DE$Exp.Global$Summary}, and \code{MarvelObject$DE$Exp.Global$Plot}
}
\description{
Volcano plot of differential splicing analysis results based on all expressed genes between 2 groups of cells. x-axis represents the log2 fold change in gene expression. y-axis represents the adjusted p-values.
}
\examples{
marvel.demo <- readRDS(system.file("extdata/data", "marvel.demo.rds", package="MARVEL"))

marvel.demo <- PlotDEValues.Exp.Global(MarvelObject=marvel.demo,
                                       pval=0.10,
                                       log2fc=0.5
                                       )

# Check output
head(marvel.demo$DE$Exp.Global$Table)
marvel.demo$DE$Exp.Global$Plot
marvel.demo$DE$Exp.Global$Summary
}
