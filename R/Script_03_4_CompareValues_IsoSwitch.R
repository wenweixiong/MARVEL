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
#' @param cell.types Character string. To indicate which 2 groups of cells that will be used for differential splicing analysis. Group names should match those in \code{cell.type} column of \code{$SplicePheno} slot.
#' @param method Character string. Statistical test for differential gene expression analysis. Can take values \code{"ks"}, \code{"wilcox"} or \code{"t.test"}. Please refer to \code{CompareValues} function for more details.
#' @param method.adjust Character strings. Adjust for multiple testing as per \code{p.adjust} function.
#' @param gene.de.sig Numeric value. Adjusted p-value below which the gene is considered differentially expressed.
#' @export
#' @return An object of class S3 containing all the original slots as inputted by the user in addition to three new slots named \code{$DE$Cor*}. \code{$DE$Cor} Original data frame generated from \code{CompareValues} function with an additional column to indicate if isoform switching has taken place. \code{$DE$CorProp} Tabulated proportion for each PSI-gene expression relationship. \code{$DE$CorPlot} Doughnut plot representing the values in \code{$DE$CorProp}.
#' @author Sean Wen <sean.wenwx@gmail.com>
#' @importFrom plyr join
#' @import stats
#' @import methods
#' @import ggplot2
#' @examples
#' marvel <- readRDS(system.file("extdata/Data", "MarvelObject.rds", package="MARVEL"))
#'
#' marvel <- IsoSwitch(MarvelObject=marvel,
#'                     psi.de.sig=0.05,
#'                     cell.types=c("iPSC", "Endoderm"),
#'                     method="t.test",
#'                     method.adjust="fdr",
#'                     gene.de.sig=0.05
#'                     )
#'
#' marvel$DE$Cor[1:5, ]
#' marvel$DE$CorPlot
#' marvel$DE$CorProp

IsoSwitch <- function(MarvelObject, psi.de.sig, cell.types, method, method.adjust, gene.de.sig) {

    # Define arguments
    de <- MarvelObject$DE$PSI
    df <- MarvelObject$Exp
    df.pheno <- MarvelObject$GenePheno
    df.feature <- MarvelObject$GeneFeature
    cell.types <- cell.types
    method <- method
    method.adjust <- method.adjust
    gene.de.sig <- gene.de.sig
    psi.de.sig <- psi.de.sig
    
    #de <- marvel$DE$PSI
    #df <- marvel$Exp
    #df.pheno <- marvel$GenePheno
    #df.feature <- marvel$GeneFeature
    #cell.types <- c("iPSC", "Endoderm")
    #method <- "t.test"
    #method.adjust <- "fdr"
    #gene.de.sig <- 0.05
    #psi.de.sig <- 0.05
         
    # Create row names for matrix
    row.names(df) <- df$gene_id
    df$gene_id <- NULL
    
    # Subset overlapping samples in matrix and pheno file
    df <- df[, which(names(df) %in% df.pheno$sample.id)]
    
    # Subset sig psi genes
    gene_ids <- de[which(de$p.val.adj < psi.de.sig), "gene_id"]
    gene_ids <- unique(gene_ids)
    df.feature <- df.feature[which(df.feature$gene_id %in% gene_ids), ]
    df <- df[df.feature$gene_id, ]
    
    # Check if matrix column and rows align with metadata
        # Column
        index.check <- which(unique((names(df)==df.pheno$sample.id))==FALSE)
        
        if(length(index.check)==0) {
            
            print("Checking... Matrix column (sample) names match sample metadata")
            
        } else {
            
            print("Checking... Matrix column (sample) names DO NOT match sample metadata")
            
        }
        
        # Row
        index.check <- which(unique((row.names(df)==df.feature$tran.id))==FALSE)
        
        if(length(index.check)==0) {
            
            print("Checking... Matrix row (feature) names match feature metadata")
            
        } else {
            
            print("Checking... Matrix row (feature) names DO NOT match feature metadata")
            
        }
    
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
        
        # Subset relevant cell types
        . <- join(., df.pheno[,c("sample.id", "cell.type")], by="sample.id", type="left")
        . <- .[which(.$cell.type %in% c(cell.types[1], cell.types[2])), ]
        
        # Retrieve values
        x <- .[which(.$cell.type==cell.types[1]), "exp"]
        y <- .[which(.$cell.type==cell.types[2]), "exp"]
        
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
        
            p.val[i] <- wilcox.test(exp ~ cell.type, .)$p.value
        
         
        } else {
        
            p.val[i] <- t.test(exp ~ cell.type, .)$p.value
        
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
        de.small$cor[which(de.small$psi.direction=="Up" & de.small$gene.direction=="Up")] <- "Same direction"
        de.small$cor[which(de.small$psi.direction=="Down" & de.small$gene.direction=="Down")] <- "Same direction"
        de.small$cor[which(de.small$psi.direction=="Up" & de.small$gene.direction=="Down")] <- "Opposite direction"
        de.small$cor[which(de.small$psi.direction=="Down" & de.small$gene.direction=="Up")] <- "Same direction"
        de.small$cor[which(de.small$psi.direction=="Up" & de.small$gene.direction=="No change")] <- "Iso-switch"
        de.small$cor[which(de.small$psi.direction=="Down" & de.small$gene.direction=="No change")] <- "Iso-switch"
        
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
            de.small$cor[which(de.small$psi.direction=="Up" & de.small$gene.direction=="Up")] <- "Same direction"
            de.small$cor[which(de.small$psi.direction=="Down" & de.small$gene.direction=="Down")] <- "Same direction"
            de.small$cor[which(de.small$psi.direction=="Up" & de.small$gene.direction=="Down")] <- "Opposite direction"
            de.small$cor[which(de.small$psi.direction=="Down" & de.small$gene.direction=="Up")] <- "Same direction"
            de.small$cor[which(de.small$psi.direction=="Up" & de.small$gene.direction=="No change")] <- "Iso-switch"
            de.small$cor[which(de.small$psi.direction=="Down" & de.small$gene.direction=="No change")] <- "Iso-switch"
            
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
        .$psi.gene.cor <- factor(.$psi.gene.cor, levels=c("Same direction", "Opposite direction", "Iso-switch", "Mixed"))
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
        legendtitle <- "Gene-PSI Correlation"
        
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
    MarvelObject$DE$Cor <- de
    MarvelObject$DE$CorProp <- .[,c("psi.gene.cor", "freq", "pct")]
    MarvelObject$DE$CorPlot <- plot
  
    return(MarvelObject)
        
}
