% Generated by roxygen2: do not edit by hand
% Please edit documentation in
%   R/Script_DROPLET_08_ADHOC_PLOT_BOXPLOT_2_PlotValues_PseudoBulk_Gene.R
\name{PlotValues.Gene.Pseudobulk.10x}
\alias{PlotValues.Gene.Pseudobulk.10x}
\title{Boxplot of pseudo-bulk PSI values vs cell groups}
\usage{
PlotValues.Gene.Pseudobulk.10x(
  MarvelObject,
  cell.group.list,
  gene_short_name,
  log2.transform = TRUE,
  cell.group.colors = NULL,
  xlabels.size = 10,
  min.n.cells.total = 10,
  method = "t.test",
  p.adjust.method = "fdr"
)
}
\arguments{
\item{MarvelObject}{Marvel object. S3 object generated from \code{CheckAlignment.10x} function.}

\item{cell.group.list}{List of character strings. Each element in the list (main list) is a list (sub-list). Each sub-list in turn is a list of pseudobulk samples represented by cell IDs. Each element in the main list is a cell group and the label of each element of the main list represents the cell group label.}

\item{gene_short_name}{Character string. Name of gene to plot.}

\item{log2.transform}{Logical value. If set to \code{TRUE} (default), normalised gene expression values will be off-set by 1 and then log2-transformed prior to plotting.}

\item{cell.group.colors}{Vector of character strings. Colors of cell groups and should be same length as \code{cell.group.list}. Default \code{ggplot2} colors are used.}

\item{xlabels.size}{Numeric value. Font size of x-tick labels. Default is \code{10}.}

\item{min.n.cells.total}{Numeric value. Minimum number of cells in a pseudobulk, below which, the pseudobulk will be omitted from plotting. Default is \code{10}.}

\item{method}{Character string. Statistical test for all possible pair-wise comparisons. Options are \code{"t.test"} (default) or \code{"wilcox"}.}

\item{p.adjust.method}{Character string. Method for multiple testing adjustment as per \code{method} option of \code{p.adjust} function. Default is \code{"fdr"}.}
}
\value{
An object of class S3 with new slots \code{MarvelObject$adhocPlot$Boxplot$Pseudobulk$Gene$Plot}, \code{MarvelObject$adhocPlot$Boxplot$Pseudobulk$Gene$Stats}, and \code{MarvelObject$adhocPlot$Boxplot$Pseudobulk$Gene$Data}.
}
\description{
Boxplot of pseudo-bulk PSI values on the y-axis against cell groups on the x-axis.
}
