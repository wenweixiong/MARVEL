#' @title Detect Isoform Switch
#'
#' @description
#' \code{IsoSwitch} compare changes in PSI values with corresponding gene expression values to detect genes that have undergone isoform switching.
#'
#' @details
#' This function compare changes in PSI values with corresponding gene expression values to detect genes that have undergone isoform switching. Isoform switch occurs when there is significant change in PSI values between 2 groups of cells in the absence of any significant change in gene expression values. This function also detect PSI values that change in the same or opposite direction with gene expression values.
#'
#' @param MarvelObject S3 object generated from \code{CompareValues} function.
#' @param psi.de.sig Numeric value. Adjusted p-value below which the splicing event is considered differentially spliced and included for isoform switching analysis.
#' @param method Character string. Statistical test for differential gene expression analysis. Can take values \code{"ks"}, \code{"wilcox"} or \code{"t.test"}. Please refer to \code{CompareValues} function for more details.
#' @param method.adjust Character strings. Adjust for multiple testing as per \code{p.adjust} function.
#' @param gene.de.sig Numeric value. Adjusted p-value below which the gene is considered differentially expressed.
#'
#' @export
#'
#' @return An object of class S3 containing all the original slots as inputted by the user in addition to three new slots named \code{$DE$Cor*}. \code{$DE$Cor} Original data frame generated from \code{CompareValues} function with an additional column to indicate if isoform switching has taken place. \code{$DE$CorProp} Tabulated proportion for each PSI-gene expression relationship. \code{$DE$CorPlot} Doughnut plot representing the values in \code{$DE$CorProp}.
#'
#' @author Sean Wen <sean.wenwx@gmail.com>
#'
#' @importFrom plyr join
#' @import stats
#' @import methods
#' @import ggplot2
#'
#' @examples
#' # Load input
#' marvel <- readRDS(system.file("extdata/Data", "MarvelObject.rds", package="MARVEL"))
#'
#' # Run example
#' marvel <- IsoSwitch(MarvelObject=marvel,
#'                     psi.de.sig=0.10,
#'                     method="t.test",
#'                     method.adjust="fdr",
#'                     gene.de.sig=0.10
#'                     )
#'
#' # Check output
#' marvel$DE$Cor$Table
#' marvel$DE$Cor$Plot
#' marvel$DE$Cor$Plot.Stats

IsoSwitch <- function(MarvelObject, psi.de.sig, method, method.adjust, gene.de.sig) {

    # Define arguments
    de <- MarvelObject$DE$PSI$Table
    df <- MarvelObject$Exp
    df.pheno <- MarvelObject$GenePheno
    df.feature <- MarvelObject$GeneFeature
    cell.type.columns.1 <- MarvelObject$DE$PSI$cell.type.columns.1
    cell.type.variables.1 <- MarvelObject$DE$PSI$cell.type.variables.1
    cell.type.columns.2 <- MarvelObject$DE$PSI$cell.type.columns.2
    cell.type.variables.2 <- MarvelObject$DE$PSI$cell.type.variables.2
    method <- method
    method.adjust <- method.adjust
    gene.de.sig <- gene.de.sig
    psi.de.sig <- psi.de.sig
    
    # Example arguments
    #de <- marvel$DE$PSI$Table
    #df <- marvel$Exp
    #df.pheno <- marvel$GenePheno
    #df.feature <- marvel$GeneFeature
    #cell.type.columns.1 <- marvel$DE$PSI$cell.type.columns.1
    #cell.type.variables.1 <- marvel$DE$PSI$cell.type.variables.1
    #cell.type.columns.2 <- marvel$DE$PSI$cell.type.columns.2
    #cell.type.variables.2 <- marvel$DE$PSI$cell.type.variables.2
    #method <- "wilcox"
    #method.adjust <- "fdr"
    #gene.de.sig <- 0.10
    #psi.de.sig <- 0.10
         
    # Create row names for matrix
    row.names(df) <- df$gene_id
    df$gene_id <- NULL
        
    # Subset sig psi genes
    gene_ids <- de[which(de$p.val.adj < psi.de.sig), "gene_id"]
    gene_ids <- unique(gene_ids)
    df.feature <- df.feature[which(df.feature$gene_id %in% gene_ids), ]
    df <- df[df.feature$gene_id, ]
    
    # Subset sample IDs (for DE later): Group 1
    .list <- list()
    
    for(i in 1:length(cell.type.columns.1)) {
        
        .list[[i]] <- df.pheno[which(df.pheno[[cell.type.columns.1[i]]] %in% cell.type.variables.1[[i]]), "sample.id"]
        
    }
    
    sample.ids.1 <- Reduce(intersect, .list)
    
    # Subset sample IDs (for DE later): Group 1
    .list <- list()
    
    for(i in 1:length(cell.type.columns.2)) {
        
        .list[[i]] <- df.pheno[which(df.pheno[[cell.type.columns.2[i]]] %in% cell.type.variables.2[[i]]), "sample.id"]
        
    }
    
    sample.ids.2 <- Reduce(intersect, .list)

    ##############################################################################
    
    # Statistical test
    gene_ids <- df.feature$gene_id
    
    n.cells.x <- NULL
    n.cells.y <- NULL
    mean.x <- NULL
    mean.y <- NULL
    mean.fc <- NULL
    p.val <- NULL

    #pb <- txtProgressBar(1, length(gene_ids), style=3)

    for(i in 1:length(gene_ids)) {

        # Subset event
        . <- df[gene_ids[i], ]
        . <- as.data.frame(t(.))
        . <- na.omit(.)
        names(.) <- "exp"
        .$sample.id <- row.names(.)
        row.names(.) <- NULL
        
        # Retrieve values
        x <- .[which(.$sample.id %in% sample.ids.1), "exp"]
        y <- .[which(.$sample.id %in% sample.ids.2), "exp"]
        
        # Compute statistics
        n.cells.x[i] <- length(which(x > 0))
        n.cells.y[i] <- length(which(y > 0))
        mean.x[i] <- mean(x)
        mean.y[i] <- mean(y)
        mean.fc[i] <- mean(y) - mean(x)
        
        # Statistical test
        if(method=="ks") {
        
            p.val[i] <- ks.test(x, y)$p.value
        
        } else if(method=="wilcox") {
        
            p.val[i] <- wilcox.test(x, y)$p.value
        
         
        } else {
        
            p.val[i] <- t.test(x, y)$p.value
        
        }
        
        # Track progress
        #setTxtProgressBar(pb, i)
        
    }
        
    # Save into data frame
    results <- data.frame("gene_id"=gene_ids,
                          "n.cells.g1"=n.cells.x, "n.cells.g2"=n.cells.y,
                          "mean.g1"=mean.x, "mean.g2"=mean.y,
                          "mean.fc"=mean.fc,
                          "p.val"=p.val,
                          stringsAsFactors=FALSE)
    
    # Reorder by p-value
    results <- results[order(results$p.val), ]
    
    # Adjust for multiple testing
    results$p.val.adj <- p.adjust(results$p.val, method=method.adjust, n=length(results$p.val))
    
    # Annotate with feature metadata
    results <- join(results, df.feature, by="gene_id", type="left")
    cols.1 <- names(df.feature)
    cols.2 <- setdiff(names(results), names(df.feature))
    results <- results[,c(cols.1, cols.2)]
    
    ##############################################################################
    
    # Indicate psi direction
    de <- de[which(de$p.val.adj < psi.de.sig), ]
    de$psi.direction <- ifelse(de$mean.diff > 0, "Up", "Down")
    
    # Indicate gene direction
    results$gene.direction <- ifelse(results$mean.fc > 0, "Up", "Down")
    results$gene.direction[which(results$p.val.adj > gene.de.sig)] <- "No change"
    
    # Merge
    de <- join(de, results[,c("gene_id", "gene.direction")], by="gene_id", type="left")
    
    # Tabulate gene id freq
    freq <- as.data.frame(table(de$gene_id), stringsAsFactors=FALSE)
    names(freq) <- c("gene_id", "freq")
    
    # psi-gene cor:single entries
        # Subset single gene ids
        gene_ids <- freq[which(freq$freq == 1), "gene_id"]
        de.small <- de[which(de$gene_id %in% gene_ids), ]
        
        # Stratify cor
        de.small$cor <- NA
        de.small$cor[which(de.small$psi.direction=="Up" & de.small$gene.direction=="Up")] <- "Positive"
        de.small$cor[which(de.small$psi.direction=="Down" & de.small$gene.direction=="Down")] <- "Positive"
        de.small$cor[which(de.small$psi.direction=="Up" & de.small$gene.direction=="Down")] <- "Negative"
        de.small$cor[which(de.small$psi.direction=="Down" & de.small$gene.direction=="Up")] <- "Negative"
        de.small$cor[which(de.small$psi.direction=="Up" & de.small$gene.direction=="No change")] <- "Iso-Switch"
        de.small$cor[which(de.small$psi.direction=="Down" & de.small$gene.direction=="No change")] <- "Iso-Switch"
        
        # Save as new object
        de.single <- de.small
    
    # psi-gene cor:duplicate entries
        # Subset single gene ids
        gene_ids <- freq[which(freq$freq != 1), "gene_id"]
        
        if(length(gene_ids) != 0) {
            
            gene_ids <- unique(gene_ids)
            de.small <- de[which(de$gene_id %in% gene_ids), ]
            
            # Stratify cor
            de.small$cor <- NA
            de.small$cor[which(de.small$psi.direction=="Up" & de.small$gene.direction=="Up")] <- "Positive"
            de.small$cor[which(de.small$psi.direction=="Down" & de.small$gene.direction=="Down")] <- "Positive"
            de.small$cor[which(de.small$psi.direction=="Up" & de.small$gene.direction=="Down")] <- "Negative"
            de.small$cor[which(de.small$psi.direction=="Down" & de.small$gene.direction=="Up")] <- "Negative"
            de.small$cor[which(de.small$psi.direction=="Up" & de.small$gene.direction=="No change")] <- "Iso-Switch"
            de.small$cor[which(de.small$psi.direction=="Down" & de.small$gene.direction=="No change")] <- "Iso-Switch"
            
            # Stratify cor: 2nd-pass
            cor.revised <- NULL
            
            for(i in 1:length(gene_ids)) {
                
                # Subset gene
                . <- de.small[which(de.small$gene_id == gene_ids[i]), ]
                
                # Stratify cor
                test.unique <- unique(.$cor)
                
                if(length(test.unique) != 1) {
                    
                    cor.revised[i] <- paste(test.unique, collapse=" + ")
                    cor.revised[i] <- "Mixed"
                    
                } else {
                
                    cor.revised[i] <- test.unique
                
                }
                
            }
            
            . <- data.frame("gene_id"=gene_ids, "cor"=cor.revised, stringsAsFactors=FALSE)
            
            de.small$cor <- NULL
            de.small <- join(de.small, ., by="gene_id", type="left")
            
            # Save as new object
            de.multi <- de.small
        
            # Merge
            de.small <- rbind.data.frame(de.single, de.multi)
        
        } else {
            
            de.small <- de.single
        }
    
    # Annotate original data frame
    de <- join(de, de.small[,c("tran_id", "cor")], by="tran_id", type="left")
    names(de)[which(names(de)=="cor")] <- "psi.gene.cor"
    
    # Remove intermediate columns
    de$psi.direction <- NULL
    de$gene.direction <- NULL
    
    # Doughnut plot
        # Tabulate freq
        . <- as.data.frame(table(de$psi.gene.cor), stringsAsFactors=FALSE)
        names(.) <- c("psi.gene.cor", "freq")
        .$pct <- .$freq / sum(.$freq) * 100
        
        # Set factor levels
        .$psi.gene.cor <- factor(.$psi.gene.cor, levels=c("Positive", "Negative", "Iso-Switch", "Mixed"))
        . <- .[order(.$psi.gene.cor), ]
        
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
        z <- data$psi.gene.cor
        maintitle <- ""
        xtitle <- ""
        ytitle <- ""
        legendtitle <- "Gene-Splicing Correlation"
        
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
    MarvelObject$DE$Cor$Table <- de
    MarvelObject$DE$Cor$Plot <- plot
    MarvelObject$DE$Cor$Plot.Stats <- .[,c("psi.gene.cor", "freq", "pct")]
  
    return(MarvelObject)
        
}
