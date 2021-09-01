#' @title Plot Gene Expression Values
#'
#' @description
#' \code{PlotValues.Exp} plots gene expression values across different groups of cells.
#'
#' @details
#' This function plots gene expression values across different groups of cells. Boxplot is used for gene expression values.
#'
#' @param MarvelObject S3 object generated from \code{CreateMarvelObject} or \code{ComputePSI} function.
#' @param cell.type.columns List. To indicate which columns in the \code{MarvelObject$GenePheno} slot to refer to when filtering samples (cells) for analysis. Each item in the list defines the columns for each group to plot.
#' @param cell.type.variables List. To indicate which specific variables of the corresponding columns to keep the samples (cells). Each item in the list should be the same length as the corresponding item in the \code{cell.type.columns} list.
#' @param cell.type.labels Character string. To indicate the cell group labels on the x-axis. Should be same length as the number of items in \code{cell.type.columns} and \code{cell.type.variables} lists.
#' @param feature Character string. \code{gene_id} for plotting. Should match \code{gene_id} column of \code{MarvelObject$GeneFeature} slot.
#' @param maintitle Character string. Column to use as plot main title as per \code{MarvelObject$GeneFeature}. Default is \code{"gene_short_name"} column.
#' @param xlabels.size Numeric value. Size of x-axis labels as per \code{ggplot2} function. Default is 8.
#' @export
#' @return An object of class S3 containing all the original slots as inputted by the user in addition to one new slot named \code{MarvelObject$adhocPlot$Exp}.
#' @author Sean Wen <sean.wenwx@gmail.com>
#' @importFrom plyr join
#' @import stats
#' @import methods
#' @import ggplot2
#' @examples
#' # Load input
#' marvel <- readRDS(system.file("extdata/Data", "MarvelObject.rds", package="MARVEL"))
#'
#' # Define event to plot
#' gene_id <- marvel$DE$PSI$Table$gene_id[1]
#'
#' # Run example
#' marvel <- PlotValues.Exp(MarvelObject=marvel,
#'                          cell.type.columns=list(c("cell.type"), c("cell.type")),
#'                          cell.type.variables=list(list("iPSC"), list("Endoderm")),
#'                          cell.type.labels=c("iPSC", "Endoderm"),
#'                          feature=gene_id,
#'                          xlabels.size=9
#'                          )
#'
#' # Check output
#' marvel$adhocPlot$Exp

PlotValues.Exp <- function(MarvelObject, cell.type.columns, cell.type.variables, cell.type.labels, feature, maintitle="gene_short_name", xlabels.size=8) {
    
    # Define arguments
    df <- MarvelObject$Exp
    df.pheno <- MarvelObject$GenePheno
    df.feature <- MarvelObject$GeneFeature
    cell.type.columns <- cell.type.columns
    cell.type.variables <- cell.type.variables
    cell.type.labels <- cell.type.labels
    feature <- feature
    maintitle <- maintitle
    xlabels.size <- xlabels.size
    
    # Example arguments
    #df <- marvel$Exp
    #df.pheno <- marvel$GenePheno
    #df.feature <- marvel$GeneFeature
    #cell.type.columns <- list(c("sample.type", "cell.type"), c("sample.type", "cell.type"))
    #cell.type.variables <- list(list("Single Cell", "iPSC"), list("Single Cell", "Endoderm"))
    #cell.type.labels <- c("Group 1", "Group 2")
    #feature <- "ENSG00000128739.22"
    #maintitle <- "gene_short_name"
    
    # Create row names for matrix
    row.names(df) <- df$gene_id
    df$gene_id <- NULL
  
    # Subset relevant feature
    df.feature.small <- df.feature[which(df.feature$gene_id ==  feature), ]
    df.small <- df[feature, , drop=FALSE]
    df.small <- as.data.frame(t(df.small))
    names(df.small) <- "exp"
    df.small$sample.id <- row.names(df.small)
    row.names(df.small) <- NULL
    
    # Retrieve sample ids for each group
    .list.2 <- list()
    
    for(j in 1:length(cell.type.columns)) {
        
        cell.type.column <- cell.type.columns[[j]]
        cell.type.variable <- cell.type.variables[[j]]
        
        .list <- list()
        
        for(i in 1:length(cell.type.column)) {
        
            .list[[i]] <- df.pheno[which(df.pheno[[cell.type.column[i]]] %in% cell.type.variable[[i]]), "sample.id"]
        
        }
        
        sample.ids <- Reduce(intersect, .list)
        
        .list.2[[j]] <- data.frame("sample.id"=sample.ids, cell.type.label=cell.type.labels[j], stringsAsFactors=FALSE)
        
    }
    
    md <- do.call(rbind.data.frame, .list.2)
    
    # Set factor levels
    md$cell.type.label <- factor(md$cell.type.label, levels=cell.type.labels)
    
    # Annotate group labels
    df.small <- join(df.small, md, by="sample.id", type="left")
    
    # Remove un-defined samples
    df.small <- df.small[!is.na(df.small$cell.type.label), ]
    
    # Compute statistics
        # n cells per cell type
        n.cells <- tapply(df.small$exp, df.small$cell.type.label, function(x) {sum(x > 0)})
        n.cells <- data.frame("cell.type.label"=names(n.cells), "freq"=n.cells, stringsAsFactors=FALSE)
        row.names(n.cells) <- NULL
        
        n.cells$label <- paste(n.cells$cell.type.label, "\n", "(n=", n.cells$freq, ")", sep="")
        
        # Average
        ave <- tapply(df.small$exp, df.small$cell.type.label, function(x) {mean(x, na.rm=TRUE)})
        ave <- data.frame("cell.type.label"=names(ave), "average"=ave, stringsAsFactors=FALSE)
        row.names(ave) <- NULL

    # Boxplot
        # Definition
        data <- df.small
        x <- data$cell.type.label
        y <- data$exp
        z <- data$cell.type.label
        maintitle <- df.feature.small[, maintitle]
        ytitle <- "Normalized expression"
        xtitle <- ""
        xlabels <- n.cells$label
        #fivenum(y) ; ymin <- 0 ; ymax <- 1 ; yinterval <- 0.25

        # Plot
        plot <- ggplot() +
            geom_boxplot(data, mapping=aes(x=x, y=y, fill=z), outlier.size=0.1) +
            #geom_jitter(data, mapping=aes(x=x, y=y), position=position_jitter(width=0.1, height=0), size=0.001) +
            stat_summary(data, mapping=aes(x=x, y=y), geom="point", fun="mean", fill="red", col="black", size=2, shape=23) +
            scale_x_discrete(labels=xlabels) +
            #scale_y_continuous(breaks=seq(ymin, ymax, by=yinterval), limits=c(ymin, ymax)) +
            labs(title=maintitle, x=xtitle, y=ytitle) +
            theme(panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.background = element_blank(),
                panel.border=element_blank(),
                plot.title=element_text(hjust = 0.5, size=12),
                plot.subtitle=element_text(hjust = 0.5, size=12),
                axis.line.y.left = element_line(color="black"),
                axis.line.x = element_line(color="black"),
                axis.title=element_text(size=12),
                axis.text.x=element_text(size=xlabels.size, colour="black"),
                axis.text.y=element_text(size=10, colour="black"),
                legend.position="none"
                )

    # Save into new slow
    MarvelObject$adhocPlot$Exp <- plot
    
    return(MarvelObject)
    
}
