% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Script_PLATE_11_ADHOC_PLOT_1_PlotValues_PSI.R
\name{PlotValues.PSI}
\alias{PlotValues.PSI}
\title{Plot percent spliced-in (PSI) values}
\usage{
PlotValues.PSI(
  MarvelObject,
  cell.group.list,
  feature,
  maintitle = "gene_short_name",
  xlabels.size = 8,
  max.cells.jitter = 10000,
  max.cells.jitter.seed = 1,
  min.cells = 25,
  sigma.sq = 0.001,
  bimodal.adjust = TRUE,
  seed = 1,
  modality.column = "modality.bimodal.adj",
  scale.y.log = FALSE,
  cell.group.colors = NULL,
  point.alpha = 0.2,
  point.size = 0.001,
  plot.type = "violin",
  fix.y.axis.interval = TRUE
)
}
\arguments{
\item{MarvelObject}{Marvel object. S3 object generated from \code{TransformExpValues} function.}

\item{cell.group.list}{List of character strings. Each element of the list is a vector of cell IDs corresponding to a cell group. The name of the element will be the cell group label.}

\item{feature}{Character string. Coordinates of splicing event to plot.}

\item{maintitle}{Character string. Column to use as plot main title as per \code{MarvelObject$ValidatedSpliceFeature}. Default is \code{"gene_short_name"} column.}

\item{xlabels.size}{Numeric value. Size of x-axis labels as per \code{ggplot2} function. Default is 8.}

\item{max.cells.jitter}{Numeric value. Maximum number of cells for jitter points. Cells are randomly downsampled to show on jitter plot. Useful when there are large number of cells so that individual jitter points do not overcrowd the violin plot.}

\item{max.cells.jitter.seed}{Numeric value. Cells downsampled are reproducible.}

\item{min.cells}{Numeric value. The minimum no. of cells expressing the splicing event to be included for analysis. Please refer to \code{AssignModality} function help page for more details.}

\item{sigma.sq}{Numeric value. The variance threshold below which the included/excluded modality will be defined as primary sub-modality, and above which it will be defined as dispersed sub-modality. Please refer to \code{AssignModality} function help page for more details. Default is 0.001.}

\item{bimodal.adjust}{Logical. When set to \code{TRUE}, MARVEL will identify false bimodal modalities and reassign them as included/excluded modality. Please refer to \code{AssignModality} function help page for more details.}

\item{seed}{Numeric value. Ensure the \code{fitdist} function returns the same values for alpha and beta paramters each time this function is executed using the same random number generator. Please refer to \code{AssignModality} function help page for more details.}

\item{modality.column}{Character string. Can take the value \code{"modality"}, \code{"modality.var"} or \code{"modality.bimodal.adj"}. Please refer to \code{AssignModality} function help page for more details. Default is \code{"modality.bimodal.adj"}.}

\item{scale.y.log}{Logical value. Only applicable when \code{level} set to \code{"splicing"}. If set to \code{TRUE}, the y-axis of will log10-scaled. Useful when most PSI values are extremely small (< 0.02) or big (> 0.98). Default is \code{FALSE}.}

\item{cell.group.colors}{Character string. Vector of colors for the cell groups specified for PCA analysis using \code{cell.type.columns}, \code{cell.type.variable}, and \code{cell.type.labels}. If not specified, default \code{ggplot2} colors will be used.}

\item{point.alpha}{Numeric value. Transparency of the data points. Takes any values between 0-1. Default value is \code{0.2}.}

\item{point.size}{Numeric value. Size of the data points. Default value is \code{0.001}.}

\item{plot.type}{Character string. Indicate to present PSI values using \code{"violin"} (default) or \code{"boxplot"}. The former is recommended for single-cell data while the latter is recommend for bulk data.}

\item{fix.y.axis.interval}{Logical. If set to \code{TRUE} (default), y-axis values will be fixed from 0-100. Only applicable to when \code{"plot.type"} option set to \code{"boxplot"}.}
}
\value{
An object of class S3 with new slot \code{MarvelObject$adhocPlot$PSI}.
}
\description{
Violin plot of percent spliced-in (PSI) values across different groups of cells.
}
\examples{
marvel.demo <- readRDS(system.file("extdata/data", "marvel.demo.rds", package="MARVEL"))

# Define cell groups to plot
df.pheno <- marvel.demo$SplicePheno
cell.group.g1 <- df.pheno[which(df.pheno$cell.type=="iPSC"), "sample.id"]
cell.group.g2 <- df.pheno[which(df.pheno$cell.type=="Endoderm"), "sample.id"]
cell.group.list <- list(cell.group.g1, cell.group.g2)
names(cell.group.list) <- c("iPSC", "Endoderm")

# Plot
marvel.demo <- PlotValues.PSI(MarvelObject=marvel.demo,
                              cell.group.list=cell.group.list,
                              feature="chr17:8383254:8382781|8383157:-@chr17:8382143:8382315",
                              min.cells=5,
                              xlabels.size=5
                              )

# Check output
marvel.demo$adhocPlot$PSI
}
