#' @title Plot Gene Expression Values
#'
#' @description
#' \code{PlotValues.Exp} plots gene expression values across different groups of cells.
#'
#' @details
#' This function plots gene expression values across different groups of cells. Boxplot is used for gene expression values.
#'
#' @param MarvelObject S3 object generated from \code{CreateMarvelObject} or \code{ComputePSI} function.
#' @param cell.types Character string. To indicate which groups of cells that will be used for plotting. Group names should match those in \code{cell.type} column of \code{$GenePheno} slot.
#' @param feature Character string. \code{gene_id} for plotting. Should match \code{gene_id} column of \code{$GeneFeature} slot.
#' @param maintitle Character string. Column to use as plot main title as per \code{GeneFeature}.
#' @export
#' @return An object of class S3 containing all the original slots as inputted by the user in addition to one new slot named \code{$adhocPlot$Gene}.
#' @author Sean Wen <sean.wenwx@gmail.com>
#' @importFrom plyr join
#' @import stats
#' @import methods
#' @import ggplot2
#' @examples
#' marvel <- readRDS(system.file("extdata/Data", "MarvelObject.rds", package="MARVEL"))
#'
#' feature <- marvel$GeneFeature$gene_id[1]
#'
#' marvel <- PlotValues.Exp(MarvelObject=marvel,
#'                          cell.types=c("iPSC", "Endoderm"),
#'                          feature=feature,
#'                          maintitle="gene_short_name"
#'                          )
#'
#' marvel$adhocPlot$Gene


PlotValues.Exp <- function(MarvelObject, cell.types, feature, maintitle=NULL) {
    
    # Define arguments
    df <- MarvelObject$Exp
    df.pheno <- MarvelObject$GenePheno
    df.feature <- MarvelObject$GeneFeature
    cell.types <- cell.types
    feature <- feature
    maintitle <- maintitle
    
    #df <- marvel$Exp
    #df.pheno <- marvel$GenePheno
    #df.feature <- marvel$GeneFeature
    #cell.types <- c("iPSC", "Endoderm")
    #feature <- "ENSG00000128739.22"
    #maintitle <- "gene_short_name"
    
    # Create row names for matrix
    row.names(df) <- df$gene_id
    df$gene_id <- NULL
    
    # Subset overlapping samples in matrix and pheno file
    df <- df[, which(names(df) %in% df.pheno$sample.id)]

    # Subset relevant feature
    df.feature.small <- df.feature[which(df.feature$gene_id ==  feature), ]
    df.small <- df[feature, ]
    df.small <- as.data.frame(t(df.small))
    names(df.small) <- "exp"
    df.small$sample.id <- row.names(df.small)
    row.names(df.small) <- NULL
    
    # Subset relevant cell types
    df.small <- join(df.small, df.pheno[,c("sample.id", "cell.type")], by="sample.id", type="left")
    df.small <- df.small[which(df.small$cell.type %in% cell.types), ]
    df.small$cell.type <- factor(df.small$cell.type, levels=cell.types)
    
    # Compute statistics
        # n cells per cell type
        n.cells <- tapply(df.small$exp, df.small$cell.type, function(x) {sum(x > 0)})
        n.cells <- data.frame("cell.type"=names(n.cells), "freq"=n.cells, stringsAsFactors=FALSE)
        row.names(n.cells) <- NULL
        
        n.cells$label <- paste(n.cells$cell.type, "\n", "(n=", n.cells$freq, ")", sep="")
        
        # Average
        ave <- tapply(df.small$exp, df.small$cell.type, function(x) {mean(x, na.rm=TRUE)})
        ave <- data.frame("cell.type"=names(ave), "average"=ave, stringsAsFactors=FALSE)
        row.names(ave) <- NULL
        
    # Keep complete entries
    df.small <- df.small[which(!is.na(df.small$cell.type)), ]
    df.small <- df.small[which(!is.na(df.small$exp)), ]
    
    # Violin plot
        # Definition
        data <- df.small
        x <- data$cell.type
        y <- data$exp
        z <- data$cell.type
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
                plot.title=element_text(hjust = 0.5, size=15),
                plot.subtitle=element_text(hjust = 0.5, size=15),
                axis.line.y.left = element_line(color="black"),
                axis.line.x = element_line(color="black"),
                axis.title=element_text(size=15),
                axis.text=element_text(size=15),
                axis.text.x=element_text(size=11, colour="black"),
                axis.text.y=element_text(size=15, colour="black"),
                legend.position="none",
                legend.title=element_text(size=12),
                legend.text=element_text(size=12)
                )

    # Save into new slow
    MarvelObject$adhocPlot$Gene <- plot
    
    return(MarvelObject)
    
}
