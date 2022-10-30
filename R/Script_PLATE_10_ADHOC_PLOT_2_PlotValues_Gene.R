#' @title Plot gene expression values
#'
#' @description Boxplot of gene expression values across different groups of cells.
#'
#' @param MarvelObject Marvel object. S3 object generated from \code{TransformExpValues} function.
#' @param cell.group.list List of character strings. Each element of the list is a vector of cell IDs corresponding to a cell group. The name of the element will be the cell group label.
#' @param feature Character string. \code{gene_id} for plotting. Should match \code{gene_id} column of \code{MarvelObject$GeneFeature} slot.
#' @param maintitle Character string. Column to use as plot main title as per \code{MarvelObject$GeneFeature}. Default is \code{"gene_short_name"} column.
#' @param xlabels.size Numeric value. Size of x-axis labels as per \code{ggplot2} function. Default is 8.
#' @param cell.group.colors Character string. Vector of colors for the cell groups specified for PCA analysis using \code{cell.type.columns}, \code{cell.type.variable}, and \code{cell.type.labels}. If not specified, default \code{ggplot2} colors will be used.
#' @param point.alpha Numeric value. Transparency of the data points. Takes any values between 0-1. Default value is \code{0.2}.
#'
#' @return An object of class S3 with new slot \code{MarvelObject$adhocPlot$Exp}.
#'
#' @importFrom plyr join
#' @import methods
#' @import ggplot2
#' @import scales
#'
#' @export
#'
#' @examples
#' marvel.demo <- readRDS(system.file("extdata/data", "marvel.demo.rds", package="MARVEL"))
#'
# Define cell groups to plot
#' df.pheno <- marvel.demo$SplicePheno
#' cell.group.g1 <- df.pheno[which(df.pheno$cell.type=="iPSC"), "sample.id"]
#' cell.group.g2 <- df.pheno[which(df.pheno$cell.type=="Endoderm"), "sample.id"]
#' cell.group.list <- list(cell.group.g1, cell.group.g2)
#' names(cell.group.list) <- c("iPSC", "Endoderm")
#'
#' # Plot
#' marvel.demo <- PlotValues.Exp(MarvelObject=marvel.demo,
#'                               cell.group.list=cell.group.list,
#'                               feature="ENSG00000161970.15",
#'                               xlabels.size=8
#'                               )
#'
#' # Check output
#' marvel.demo$adhocPlot$Exp

PlotValues.Exp <- function(MarvelObject, cell.group.list, feature, maintitle="gene_short_name", xlabels.size=8, cell.group.colors=NULL, point.alpha=0.2) {
    
    # Define arguments
    df <- MarvelObject$Exp
    df.pheno <- MarvelObject$SplicePheno
    df.feature <- MarvelObject$GeneFeature
    cell.group.list <- cell.group.list
    feature <- feature
    maintitle <- maintitle
    xlabels.size <- xlabels.size
    cell.group.colors <- cell.group.colors
    point.alpha <- point.alpha
    
    # Example arguments
    #MarvelObject <- marvel
    #df <- MarvelObject$Exp
    #df.pheno <- MarvelObject$SplicePheno
    #df.feature <- MarvelObject$GeneFeature
    #cell.group.list <- cell.group.list
    #feature <- gene_id
    #maintitle <- "gene_short_name"
    #xlabels.size <- 8
    #cell.group.colors <- NULL
    #point.alpha <- 0.2
    
    ###############################################################
    
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
    .list <- list()
    
    for(i in 1:length(cell.group.list)) {
        
        
        .list[[i]] <- data.frame("sample.id"=cell.group.list[[i]],
                                 "cell.type.label"=names(cell.group.list)[[i]],
                                 stringsAsFactors=FALSE
                                 )
        
    }
    
    md <- do.call(rbind.data.frame, .list)
    
    # Set factor levels
    md$cell.type.label <- factor(md$cell.type.label, levels=names(cell.group.list))
    
    # Annotate group labels
    df.small <- join(df.small, md, by="sample.id", type="left")
    
    # Remove un-defined samples
    df.small <- df.small[!is.na(df.small$cell.type.label), ]
    
    # Compute statistics
        # n cells per cell type
        freq.expr <- tapply(df.small$exp, df.small$cell.type.label, function(x) {sum(x > 0)})
        freq.expr <- data.frame("cell.type.label"=names(freq.expr), "freq.expr"=as.numeric(freq.expr), stringsAsFactors=FALSE)
        freq.total <- as.data.frame(table(df.small$cell.type.label), stringsAsFactors=FALSE)
        freq.total <- data.frame("cell.type.label"=freq.total[,1], "freq.total"=freq.total[,2], stringsAsFactors=FALSE)
        
        n.cells <- join(freq.expr, freq.total, by="cell.type.label", type="left")
        
        n.cells$pct.expr <- round(n.cells$freq.expr / n.cells$freq.total * 100, digits=0)
        
        n.cells$label <- ifelse(n.cells$pct.expr > 1,
                            paste(n.cells$cell.type.label, "\n", "", n.cells$freq.expr, "/", n.cells$freq.total, " (", n.cells$pct.expr, "%)", sep=""),
                            paste(n.cells$cell.type.label, "\n", "", n.cells$freq.expr, "/", n.cells$freq.total, " (<1%)", sep="")
                            )
        
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
        
        # Color scheme
        if(is.null(cell.group.colors[1])) {
        
            gg_color_hue <- function(n) {
              hues = seq(15, 375, length = n + 1)
              grDevices::hcl(h = hues, l = 65, c = 100)[1:n]
            }
            
            n = length(levels(z))
            cols = gg_color_hue(n)
        
        } else {
            
            cols <- cell.group.colors
            
        }

        # Plot
        plot <- ggplot() +
            geom_boxplot(data, mapping=aes(x=x, y=y, fill=z), outlier.size=0.1) +
            geom_jitter(data, mapping=aes(x=x, y=y), position=position_jitter(width=0.1, height=0), size=0.001, alpha=point.alpha) +
            stat_summary(data, mapping=aes(x=x, y=y), geom="point", fun="mean", fill="red", col="black", size=2, shape=23) +
            scale_fill_manual(values=cols) +
            scale_x_discrete(labels=xlabels) +
            scale_y_continuous(breaks= pretty_breaks()) +
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
