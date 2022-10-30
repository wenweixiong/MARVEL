#' @title Plot differential splice junction analysis results for a specified gene
#'
#' @description Scatterplot of results from differential gene and splice junction analysis. x-axis represents the gene expression log2 fold change between the different pairs of cell groups. y-axis represents the PSI differences or log2 fold change between the different pairs of cell groups.
#'
#' @param MarvelObject Marvel object. S3 object generated from \code{adhocGene.DE.Gene.10x} and \code{adhocGene.DE.PSI.10x} functions.
#' @param coord.intron Character string. Coordinates of splice junction whose differential splice junction results will be plotted.
#' @param log2fc.gene Numeric value. Absolute log2 fold change, above which, the gene is considered differentially expressed.
#' @param delta.sj Numeric value. Absolute differences in average PSI values between the two cell groups, above which, the splice junction is considered differentially spliced.
#' @param label.size Numeric value. The font size of the group comparison labels on the plot will be adjusted to the size specified here. Default is \code{2}.
#' @param point.size Numeric value. Size of data points. Default is \code{2}.
#' @param xmin Numeric value. Minimum x-axis value.
#' @param xmax Numeric value. Maximum x-axis value.
#' @param ymin Numeric value. Minimum y-axis value.
#' @param ymax Numeric value. Maximum y-axis value.
#'
#' @return An object of class S3 with a new slots \code{MarvelObject$adhocGene$DE$VolcanoPlot$Plot} and \code{MarvelObject$adhocGene$DE$VolcanoPlot$Table}.
#'
#' @importFrom plyr join
#' @import ggplot2
#' @import scales
#' @importFrom grDevices hcl
#' @import Matrix
#'
#' @export
#'
#' @examples
#'
#' marvel.demo.10x <- readRDS(system.file("extdata/data",
#'                                "marvel.demo.10x.rds",
#'                                package="MARVEL")
#'                                )
#'
#' # Define SJ to plot
#' coord.intron <- marvel.demo.10x$adhocGene$DE$PSI$Data$coord.intron[1]
#'
#' # Plot SJ vs gene
#' marvel.demo.10x <- adhocGene.PlotDEValues.10x(
#'                         MarvelObject=marvel.demo.10x,
#'                         coord.intron=coord.intron,
#'                         log2fc.gene=0.5,
#'                         delta.sj=5,
#'                         label.size=2,
#'                         point.size=2,
#'                         xmin=-2.0,
#'                         xmax=2.0,
#'                         ymin=-25,
#'                         ymax=25
#'                         )
#'
#' # Check output
#' marvel.demo.10x$adhocGene$DE$VolcanoPlot$Plot

adhocGene.PlotDEValues.10x <- function(MarvelObject, coord.intron, log2fc.gene=0.5, delta.sj=5, label.size=2, point.size=2, xmin=NULL, xmax=NULL, ymin=NULL, ymax=NULL) {
        
    # Define arguments
    MarvelObject <- MarvelObject
    results.gene <- MarvelObject$adhocGene$DE$Gene$Data
    results.sj <- MarvelObject$adhocGene$DE$PSI$Data
    coord.intron <- coord.intron
    log2fc.gene <- log2fc.gene
    delta.sj <- delta.sj
    label.size <- label.size
    point.size <- point.size
    xmin <- xmin
    xmax <- xmax
    ymin <- ymin
    ymax <- ymax
    
    # Example arguments
    #MarvelObject <- marvel
    #results.gene <- MarvelObject$adhocGene$DE$Gene$Data
    #results.sj <- MarvelObject$adhocGene$DE$PSI$Data
    #coord.intron <- MarvelObject$adhocGene$DE$PSI$Table$coord.intron[1]
    #log2fc.gene <- 0.5
    #delta.sj <- 1
    #label.size <- 2.5
    #point.size <- 2
    #xmin <- -4.5
    #xmax <- 4.5
    #ymin <- -20
    #ymax <- 20
    
    #################################################
    
    # Subset define SJ
    results.sj <- results.sj[which(results.sj$coord.intron==coord.intron), ]
    
    # Merge DE gene, SJ
    results <- join(results.gene, results.sj[,c("group.pair", "delta")], by="group.pair", type="left")
    
    # Indicate sig genes/sj
    results$change <- NA
    
    index <- which(results$log2fc > log2fc.gene & results$delta > delta.sj)
    results$change[index] <- "Coordinated"
    
    index <- which(results$log2fc < (log2fc.gene * -1) & results$delta < (delta.sj * -1))
    results$change[index] <- "Coordinated"
    
    index <- which(results$log2fc > log2fc.gene & results$delta < (delta.sj * -1))
    results$change[index] <- "Opposing"
  
    index <- which(results$log2fc < (log2fc.gene * -1) & results$delta > delta.sj)
    results$change[index] <- "Opposing"
    
    index <- which(is.na(results$change) & abs(results$delta) > delta.sj)
    results$change[index] <- "Iso-Swicth"
    
    index <- which(is.na(results$change))
    results$change[index] <- "Gene/SJ n.s."
      
    # Set factor levels
    levels <- intersect(c("Coordinated", "Opposing", "Iso-Swicth", "Gene/SJ n.s."), unique(results$change))
    results$change <- factor(results$change, levels=levels)
    
    # Color scheme
    gg_color_hue <- function(n) {
        hues = seq(15, 375, length = n + 1)
        hcl(h = hues, l = 65, c = 100)[1:n]
    }
    
    cols = gg_color_hue(4)[c(1:3)]
    
    color <- data.frame("change"=c("Coordinated", "Opposing", "Iso-Swicth", "Gene/SJ n.s."),
                        "change.color"=c(cols, "gray"),
                        stringsAsFactors=FALSE
                        )
                        
    color <- color[which(color$change %in% levels), "change.color"]
    
    # Volcano plot
        # Definition
        data <- results
        x <- data$log2fc
        y <- data$delta
        z <- data$change
        maintitle <- coord.intron
        xtitle <- "log2FC expression"
        ytitle <- "Delta PSI"
        label <- data$group.pair
                        
        if(is.null(xmin)) {
            
            # Plot
            plot <- ggplot() +
                geom_point(data, mapping=aes(x=x, y=y, fill=z), size=point.size, alpha=0.75, pch=21, stroke=0.1) +
                geom_hline(yintercept=c(delta.sj * -1, delta.sj), linetype="dashed", color="black", size=0.1) +
                geom_vline(xintercept=c(log2fc.gene * -1, log2fc.gene), linetype="dashed", color="black", size=0.1) +
                ggrepel::geom_text_repel(data, mapping=aes(x=x, y=y, label=label), max.overlaps = Inf, box.padding = 0.5, size=label.size, max.time = 1, max.iter = 1e5, segment.alpha=0.5, segment.size=0.1, min.segment.length = 0.01) +
                scale_fill_manual(values=color) +
                #scale_x_continuous(limits=c(xmin, xmax)) +
                #scale_y_continuous(limits=c(ymin, ymax)) +
                labs(title=maintitle, x=xtitle, y=ytitle) +
                theme(panel.grid.major=element_blank(),
                    panel.grid.minor=element_blank(),
                    panel.background=element_blank(),
                    plot.title = element_text(size=12, hjust=0.5),
                    axis.line=element_line(colour = "black"),
                    axis.title=element_text(size=12),
                    axis.text.x=element_text(size=10, colour="black"),
                    axis.text.y=element_text(size=10, colour="black"),
                    legend.title=element_blank(),
                    legend.text=element_text(size=8),
                    legend.key=element_blank()
                    )  +
                guides(fill = guide_legend(override.aes=list(size=2, alpha=0.75, stroke=0.1), ncol=1))

        } else {
            
            # Plot
            plot <- ggplot() +
                geom_point(data, mapping=aes(x=x, y=y, fill=z), size=point.size, alpha=0.75, pch=21, stroke=0.1) +
                geom_hline(yintercept=c(delta.sj * -1, delta.sj), linetype="dashed", color="black", size=0.1) +
                geom_vline(xintercept=c(log2fc.gene * -1, log2fc.gene), linetype="dashed", color="black", size=0.1) +
                ggrepel::geom_text_repel(data, mapping=aes(x=x, y=y, label=label), max.overlaps = Inf, box.padding = 0.5, size=label.size, max.time = 1, max.iter = 1e5, segment.alpha=0.5, segment.size=0.1, min.segment.length = 0.01) +
                scale_fill_manual(values=color) +
                scale_x_continuous(limits=c(xmin, xmax)) +
                scale_y_continuous(limits=c(ymin, ymax)) +
                labs(title=maintitle, x=xtitle, y=ytitle) +
                theme(panel.grid.major=element_blank(),
                    panel.grid.minor=element_blank(),
                    panel.background=element_blank(),
                    plot.title = element_text(size=12, hjust=0.5),
                    axis.line=element_line(colour = "black"),
                    axis.title=element_text(size=12),
                    axis.text.x=element_text(size=10, colour="black"),
                    axis.text.y=element_text(size=10, colour="black"),
                    legend.title=element_blank(),
                    legend.text=element_text(size=8),
                    legend.key=element_blank()
                    )  +
                guides(fill = guide_legend(override.aes=list(size=2, alpha=0.75, stroke=0.1), ncol=1))
                        
        }

    # Save into new slot
    MarvelObject$adhocGene$DE$VolcanoPlot$Plot <- plot
    MarvelObject$adhocGene$DE$VolcanoPlot$Table <- results
    
    # Return final object
    return(MarvelObject)
            
}


