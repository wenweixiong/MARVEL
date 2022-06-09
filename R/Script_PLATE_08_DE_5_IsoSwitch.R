#' @title Classify gene-splicing relationship
#'
#' @description Classify gene-splicing relative changes to each other from cell group 1 to group 2. Classifications are coordinated, opposing, isoform-switching, and complex. In coordinated relationship, both gene and splicing changes in the same direction from cell group 1 to group 2. In opposing relationship, gene changes in the opposite direction relative to splicing from cell group 1 to group 2. In isoform-switching, there is differential splice junction usage without differential expression of the corresponding gene between cell group 1 and group 2. Complex relationship involves genes with both coordinated and opposing relationships with splicing. Only differentially spliced junctions are included for analysis here.
#'
#' @param MarvelObject Marvel object. S3 object generated from \code{CompareValues.Genes.10x} function.
#' @param method Character string. The statistical method used for differential splicing analysis.
#' @param psi.pval Numeric value. Adjusted p-value below which the splicing event is considered differentially spliced and included for isoform switching analysis. To be used in conjunction with \code{psi.delta}.
#' @param psi.delta Numeric value. The absolute mininum difference in PSI values between the two cell groups above which the splicing event is considered differentially spliced nd included for isoform switching analysis. To be used in conjunction with \code{psi.pval}. Specify \code{0} (default) to switch this threshold off.
#' @param gene.pval Numeric value. Adjusted p-value below which the gene is considered differentially expressed. Default value is \code{0.1}.
#' @param gene.log2fc Numeric value. The absolute log2 fold change in mean gene expression values between the two cell groups above which the gene is considered differentially expressed. To be used in conjunction with \code{gene.pval}. Specify \code{0} to switch this threshold off. Default value is \code{0.5}.
#'
#' @return An object of class S3 containing with new slots \code{MarvelObject$DE$Cor$Table}, \code{MarvelObject$DE$Cor$Plot}, and \code{MarvelObject$DE$Cor$Plot.Stats}.
#'
#' @importFrom plyr join
#' @import stats
#' @import methods
#' @import ggplot2
#'
#' @export

IsoSwitch <- function(MarvelObject, method, psi.pval=0.1, psi.delta=0, gene.pval=0.1, gene.log2fc=0.5) {

    # Define arguments
    method <- method
    psi.pval <- psi.pval
    psi.delta <- psi.delta
    de.exp <- MarvelObject$DE$Exp$Table
    gene.pval <- gene.pval
    gene.log2fc <- gene.log2fc
    
    # Example arguments
    #MarvelObject <- marvel
    #method <- c("ad", "dts")
    #psi.pval <- c(0.10, 0.10)
    #de.exp <- MarvelObject$DE$Exp$Table
    #gene.pval <- 0.10
    #gene.log2fc <- 0.5
    #psi.delta <- 0
    
    # Tabulate sig events
    .list <- list()
    
    for(j in 1:length(method)) {
    
        # Subset relevent splicing DE results
        de.psi <- MarvelObject$DE$PSI$Table[[method[j]]]

        # Indicate psi direction
        de.psi$direction.psi <- NA
        index <- which(de.psi$mean.diff > psi.delta & de.psi$p.val.adj < psi.pval[j] & de.psi$outlier==FALSE)
        de.psi$direction.psi[index] <- "Up"
        index <- which(de.psi$mean.diff < (psi.delta * -1) & de.psi$p.val.adj < psi.pval[j] & de.psi$outlier==FALSE)
        de.psi$direction.psi[index] <- "Down"
        de.psi$direction.psi[is.na(de.psi$direction.psi)] <- "No change"
        
        # Subset sig events
        de.psi <- de.psi[which(de.psi$direction.psi %in% c("Up", "Down")), ]
        
        # Subset relevant columns
        de.psi <- de.psi[,c("tran_id", "event_type", "gene_id", "gene_short_name", "gene_type", "direction.psi")]
        
        # Save into list
        .list[[j]] <- de.psi
        
    }
    
    de.psi <- do.call(rbind.data.frame, .list)
    de.psi <- unique(de.psi)
        
    # Annotate gene DE results
    de.exp <- MarvelObject$DE$Exp$Table
    de.exp <- de.exp[, c("gene_id", "log2fc", "p.val.adj")]
    names(de.exp)[-1] <- paste(names(de.exp)[-1], ".gene", sep="")
    de <- join(de.psi, de.exp, by="gene_id", type="left")
        
    # Indicate gene direction
    de$direction.gene <- NA
    de$direction.gene[which(de$log2fc.gene > gene.log2fc & de$p.val.adj.gene < gene.pval)] <- "Up"
    de$direction.gene[which(de$log2fc.gene < (gene.log2fc * -1) & de$p.val.adj < gene.pval)] <- "Down"
    de$direction.gene[is.na(de$direction.gene)] <- "No change"
    
    # Tabulate gene id freq
    freq <- as.data.frame(table(de$gene_id), stringsAsFactors=FALSE)
    names(freq) <- c("gene_id", "freq")
        
    # psi-gene cor:single entries
        # Subset single gene ids
        gene_ids <- freq[which(freq$freq == 1), "gene_id"]
        de.small <- de[which(de$gene_id %in% gene_ids), ]
        
        # Stratify cor
        de.small$cor <- NA
        de.small$cor[which(de.small$direction.psi=="Up" & de.small$direction.gene=="Up")] <- "Coordinated"
        de.small$cor[which(de.small$direction.psi=="Down" & de.small$direction.gene=="Down")] <- "Coordinated"
        de.small$cor[which(de.small$direction.psi=="Up" & de.small$direction.gene=="Down")] <- "Opposing"
        de.small$cor[which(de.small$direction.psi=="Down" & de.small$direction.gene=="Up")] <- "Opposing"
        de.small$cor[which(de.small$direction.psi=="Up" & de.small$direction.gene=="No change")] <- "Iso-Switch"
        de.small$cor[which(de.small$direction.psi=="Down" & de.small$direction.gene=="No change")] <- "Iso-Switch"
        
        # Save as new object
        de.single <- de.small
        
    # psi-gene cor: duplicate entries
        # Subset duplicate gene ids
        gene_ids <- freq[which(freq$freq != 1), "gene_id"]
                        
        gene_ids <- unique(gene_ids)
        de.small <- de[which(de$gene_id %in% gene_ids), ]
        
        # Stratify cor
        de.small$cor <- NA
        de.small$cor[which(de.small$direction.psi=="Up" & de.small$direction.gene=="Up")] <- "Coordinated"
        de.small$cor[which(de.small$direction.psi=="Down" & de.small$direction.gene=="Down")] <- "Coordinated"
        de.small$cor[which(de.small$direction.psi=="Up" & de.small$direction.gene=="Down")] <- "Opposing"
        de.small$cor[which(de.small$direction.psi=="Down" & de.small$direction.gene=="Up")] <- "Opposing"
        de.small$cor[which(de.small$direction.psi=="Up" & de.small$direction.gene=="No change")] <- "Iso-Switch"
        de.small$cor[which(de.small$direction.psi=="Down" & de.small$direction.gene=="No change")] <- "Iso-Switch"
        
        # Stratify cor: Find mixed correlation
        .list <- list()
        
        for(i in 1:length(gene_ids)) {
            
            de.small. <- de.small[which(de.small$gene_id==gene_ids[i]), ]
                          
            cor.type <- unique(de.small.$cor)
            
            if(length(cor.type)==1) {
                
                de.small.$cor <- cor.type
                .list[[i]] <- de.small.
                
            } else if(length(cor.type) >= 2) {
                
                de.small.$cor <- "Complex"
                .list[[i]] <- de.small.
                
            }
            
        }
        
        # Save as new object
        de.multi <- do.call(rbind.data.frame, .list)
    
        # Merge
        de <- rbind.data.frame(de.single, de.multi)
    
    # Save gene data only
    cols <- c("gene_id", "gene_short_name", "gene_type", "cor")
    results <- de[, cols]
    results <- unique(results)

    # Doughnut plot
        # Tabulate freq (for gene, not splicing)
        #. <- as.data.frame(table(de$psi.gene.cor), stringsAsFactors=FALSE)
        #de.small <- unique(de[,c("gene_short_name", "psi.gene.cor")])
        . <- as.data.frame(table(results$cor), stringsAsFactors=FALSE)
        names(.) <- c("cor", "freq")
        .$pct <- .$freq / sum(.$freq) * 100
        
        # Set factor levels
        levels <- intersect(c("Coordinated", "Opposing", "Iso-Switch", "Complex"), unique(.$cor))
        .$cor <- factor(.$cor, levels=c("Coordinated", "Opposing", "Iso-Switch", "Complex"))
        . <- .[order(.$cor), ]
        
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
        z <- data$cor
        maintitle <- ""
        xtitle <- ""
        ytitle <- ""
        legendtitle <- "Gene-Splicing Relationship"
        
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
        
    # Save to new slots
    MarvelObject$DE$Cor$Table <- results
    MarvelObject$DE$Cor$Plot <- plot
    MarvelObject$DE$Cor$Plot.Stats <- .[,c("cor", "freq", "pct")]
  
    return(MarvelObject)
        
}
