#' @title Classify gene-splicing relationship
#'
#' @description Classify gene-splicing relative changes to each other from cell group 1 to group 2. Classifications are coordinated, opposing, isoform-switching, and complex. In coordinated relationship, both gene and splicing changes in the same direction from cell group 1 to group 2. In opposing relationship, gene changes in the opposite direction relative to splicing from cell group 1 to group 2. In isoform-switching, there is differential splice junction usage without differential expression of the corresponding gene between cell group 1 and group 2. Complex relationship involves genes with both coordinated and opposing relationships with splicing. Only differentially spliced junctions are included for analysis here.
#'
#' @param MarvelObject Marvel object. S3 object generated from \code{CompareValues.Genes.10x} function.
#' @param pval.sj Numeric value. p-value from differential splicing analysis, below which, the splice junction is considered differentially spliced. Default is \code{0.05}.
#' @param log2fc.sj Numeric value. Absolute log2 fold change from differential splicing analysis, above which, the splice junction is considered differentially spliced. This option should be \code{NULL} if \code{delta.sj} has been specified.
#' @param delta.sj Numeric value. Absolute difference in average PSI values between the two cell groups, above which, the splice junction is considered differentially spliced. This option should be \code{NULL} if \code{log2fc.sj} has been specified.
#' @param min.gene.norm Numeric value. The average normalised gene expression across the two cell groups above which the splice junction is considered differentially spliced. Default is \code{0}.
#' @param pval.adj.gene Numeric value. Adjusted p-value from differential gene expression analysis, below which, the gene is considered differentially expressed. Default is \code{0.05}.
#' @param log2fc.gene Numeric value. Absolute log2 fold change from differential gene expression analysis, above which, the gene is considered differentially expressed. This option should be \code{NULL} if \code{delta.sj} has been specified.
#'
#' @return An object of class S3 containing new slots \code{MarvelObject$SJ.Gene.Cor$Data}, \code{MarvelObject$SJ.Gene.Cor$Proportion$Plot}, and \code{MarvelObject$SJ.Gene.Cor$Proportion$Table}.
#'
#' @importFrom plyr join
#' @import ggplot2
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
#' marvel.demo.10x <- readRDS(system.file("extdata/data",
#'                                "marvel.demo.10x.rds",
#'                                package="MARVEL")
#'                                )
#'
#' marvel.demo.10x <- IsoSwitch.10x(
#'                         MarvelObject=marvel.demo.10x,
#'                         pval.sj=0.05,
#'                         delta.sj=5,
#'                         min.gene.norm=1.0,
#'                         pval.adj.gene=0.05,
#'                         log2fc.gene=0.5
#'                         )
#'
#' # Check outputs
#' marvel.demo.10x$SJ.Gene.Cor$Proportion$Plot
#' marvel.demo.10x$SJ.Gene.Cor$Proportion$Table
#' cols <- c("coord.intron", "gene_short_name", "cor.complete")
#' head(marvel.demo.10x$SJ.Gene.Cor$Data[,cols])

IsoSwitch.10x <- function(MarvelObject, pval.sj=0.05, log2fc.sj=NULL, delta.sj=5, min.gene.norm=0, pval.adj.gene=0.05, log2fc.gene=0.5) {

    # Define arguments
    MarvelObject <- MarvelObject
    df <- MarvelObject$DE$SJ$Table
    pval.sj <- pval.sj
    log2fc.sj <- log2fc.sj
    delta.sj <- delta.sj
    pval.adj.gene <- pval.adj.gene
    log2fc.gene <- log2fc.gene
    min.gene.norm <- min.gene.norm
    
    # Example arguments
    #MarvelObject <- marvel
    #df <- MarvelObject$DE$SJ$Table
    #pval.sj <- 0.05
    #log2fc.sj <- NULL
    #delta.sj <- 5
    #pval.adj.gene <- 0.05
    #log2fc.gene <- 0.5
    #min.gene.norm <- 1
    
    ######################### CLASSIFY RELATIONSHIP ###########################
    
    # Indicate sig events and direction: SJ
    if(!is.null(log2fc.sj)) {
        
        df$sig <- NA
        df$sig[which(df$pval < pval.sj & df$log2fc > log2fc.sj & df$mean.expr.gene.norm.g1.g2 > min.gene.norm)] <- "up"
        df$sig[which(df$pval < pval.sj & df$log2fc < (log2fc.sj*-1) & df$mean.expr.gene.norm.g1.g2 > min.gene.norm)] <- "down"
        df$sig[is.na(df$sig)] <- "n.s."
        df$sig <- factor(df$sig, levels=c("up", "down", "n.s."))
        table(df$sig)
        
    } else if(!is.null(delta.sj)){
        
        df$sig <- NA
        df$sig[which(df$pval < pval.sj & df$delta > delta.sj & df$mean.expr.gene.norm.g1.g2 > min.gene.norm)] <- "up"
        df$sig[which(df$pval < pval.sj & df$delta < (delta.sj*-1) & df$mean.expr.gene.norm.g1.g2 > min.gene.norm)] <- "down"
        df$sig[is.na(df$sig)] <- "n.s."
        df$sig <- factor(df$sig, levels=c("up", "down", "n.s."))
        table(df$sig)
        
    }
    
    names(df)[which(names(df)=="sig")] <- "sig.sj"
    
    # Subset sig SJ
    df <- df[which(df$sig.sj %in% c("up", "down")), ]
    
    # Indicate sig events and direction: Gene
    df$sig <- NA
    df$sig[which(df$pval.adj.gene.norm < pval.adj.gene & df$log2fc.gene.norm > log2fc.gene)] <- "up"
    df$sig[which(df$pval.adj.gene.norm < pval.adj.gene & df$log2fc.gene.norm < (log2fc.gene*-1))] <- "down"
    df$sig[is.na(df$sig)] <- "n.s."
    df$sig <- factor(df$sig, levels=c("up", "down", "n.s."))
    table(df$sig)
    
    names(df)[which(names(df)=="sig")] <- "sig.gene"
    
    # Classify SJ-gene relationship
    table(df$sig.sj, df$sig.gene)
    
    df$cor <- NA
    
    df$cor[which(df$sig.sj=="up" & df$sig.gene=="up")] <- "Coordinated"
    df$cor[which(df$sig.sj=="down" & df$sig.gene=="down")] <- "Coordinated"
    df$cor[which(df$sig.sj=="down" & df$sig.gene=="up")] <- "Opposing"
    df$cor[which(df$sig.sj=="up" & df$sig.gene=="down")] <- "Opposing"
    df$cor[which(df$sig.sj=="up" & df$sig.gene=="n.s.")] <- "Iso-Switch"
    df$cor[which(df$sig.sj=="down" & df$sig.gene=="n.s.")] <- "Iso-Switch"
    
    sum(is.na(df$cor))
    
    # Find mixed relationship
    gene_short_names <- unique(df$gene_short_name)

    corr <- NULL
    
    for(i in 1:length(gene_short_names)) {
        
        # Subset gene
        df.small <- df[which(df$gene_short_name == gene_short_names[i]), ]
        
        # Stratify cor
        test.unique <- unique(df.small$cor)
        
        if(length(test.unique) != 1) {
            
            corr[i] <- "Complex"
            
        } else {
        
            corr[i] <- test.unique
        
        }
        
    }
    
    # Update SJ-gene cor column
    results <- data.frame("gene_short_name"=gene_short_names, "cor.complete"=corr, stringsAsFactors=FALSE)
    df <- join(df, results, by="gene_short_name", type="left")
    
    ############################ PLOT PROPORTION ###########################
    
    # Set factor levels
    levels <- intersect(c("Coordinated", "Opposing", "Iso-Switch", "Complex"), unique(df$cor.complete))

    # Tabulate freq (for genes, not splicing)
    #. <- as.data.frame(table(df$cor.complete), stringsAsFactors=FALSE)
    df.small <- unique(df[, c("gene_short_name", "cor.complete")])
    . <- as.data.frame(table(df.small$cor.complete), stringsAsFactors=FALSE)
    names(.) <- c("sj.gene.cor", "freq")
    .$pct <- .$freq / sum(.$freq) * 100
    
    # Set factor levels
    .$sj.gene.cor <- factor(.$sj.gene.cor, levels=levels)
    . <- .[order(.$sj.gene.cor), ]
    
    # Compute statistics for plot
    .$fraction <- .$freq / sum(.$freq)
    .$ymax <- cumsum(.$fraction)
    .$ymin = c(0, .$ymax[-length(.$ymax)])
    
    # Definitions
    data <- .
    xmax <- nrow(data) + 1
    xmin <- nrow(data)
    ymax <- data$ymax
    ymin <- data$ymin
    z <- data$sj.gene.cor
    maintitle <- ""
    xtitle <- ""
    ytitle <- ""
    legendtitle <- "Gene-SJ Relationship"
    
    # Plot
    plot <- ggplot() +
        geom_rect(data=data, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=z), color="black") +
        coord_polar(theta="y") +
        xlim(c(2, 4)) +
        #scale_fill_manual(values=colors) +
        labs(title=maintitle, x=xtitle, y=ytitle, fill=legendtitle) +
        theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            panel.border=element_blank(),
            plot.title=element_text(hjust = 0.5, size=15),
            plot.subtitle=element_text(hjust = 0.5, size=15),
            axis.line = element_blank(),
            axis.ticks=element_blank(),
            axis.text=element_blank(),
            legend.title=element_text(size=9),
            legend.text=element_text(size=9)
            )
            
    #########################################################################
    
    # Save into new slot
    MarvelObject$SJ.Gene.Cor$Data <- df
    MarvelObject$SJ.Gene.Cor$Proportion$Plot <- plot
    MarvelObject$SJ.Gene.Cor$Proportion$Table <- data[,c("sj.gene.cor", "freq", "pct")]
    
    # Return final object
    return(MarvelObject)
        
}
