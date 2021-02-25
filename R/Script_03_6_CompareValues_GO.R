#' @title Gene Ontology Analysis
#'
#' @description
#' \code{BioPathways} performs gene ontology analysis on genes that are differentially spliced.
#'
#' @details
#' This function performs gene ontology analysis on genes that are differentially spliced to identify significantly regulated biological pathways.
#'
#' @param MarvelObject S3 object generated from \code{CompareValues} function.
#' @param psi.de.sig Numeric value. Adjusted p-value below which the splicing event is considered differentially spliced and included for isoform switching analysis.
#' @param method.adjust Character string. Adjust p-values for multiple testing. Options available as per \code{p.adjust} function.
#' @param min.genes Numeric value. Number of differentially spliced genes required for analysis.
#' @param p.val.adj.return Numeric value. Return pathways with adjusted p-value below this value.
#' @param plot.top.n Numeric value. Plot these top most significant pathways.
#' @export
#' @return An object of class S3 containing all the original slots as inputted by the user in addition to two new slots. \code{$DE$BioPathways} contains the pathways significantly enriched among differentially spliced genes as specified in \code{p.val.adj.return}. \code{$DE$BioPathwaysPlot} contains the plot of top most significant pathways ranked by adjusted p-values as specified in \code{plot.top.n}.
#' @author Sean Wen <sean.wenwx@gmail.com>
#' @importFrom plyr join
#' @import stats
#' @import methods
#' @import ggplot2
#' @import org.Hs.eg.db
#' @import GO.db
#' @import GOstats
#' @import stringr
#' @importFrom AnnotationDbi select
#' @examples
#' marvel <- readRDS(system.file("extdata/Data", "MarvelObject.rds", package="MARVEL"))
#'
#' marvel <- BioPathways(MarvelObject=marvel,
#'                       psi.de.sig=0.05,
#'                       method.adjust="fdr",
#'                       min.genes=50,
#'                       p.val.adj.return=0.05,
#'                       plot.top.n=10
#'                       )
#'
#' marvel$DE$BioPathways
#' marvel$DE$BioPathwaysPlot

BioPathways <- function(MarvelObject, psi.de.sig, method.adjust, min.genes, p.val.adj.return, plot.top.n) {
    
    # Define arguments
    de <- MarvelObject$DE$PSI
    psi.de.sig <- psi.de.sig
    method.adjust <- method.adjust
    min.genes <- min.genes
    p.val.adj.return <- p.val.adj.return
    plot.top.n <- plot.top.n
    
    #de <- marvel$DE$PSI
    #psi.de.sig <- 0.05
    #method.adjust <- "fdr"
    #p.val.adj.return <- 0.05
    #plot.top.n <- 10
    
    # Define gene universe
    genes.universe <- unique(de$gene_short_name)

    # Subset significant genes
    genes.sig <- unique(de[which(de$p.val.adj < 0.05), "gene_short_name"])
    
    if(length(genes.sig) > min.genes) {
    
    ######################################################
    
    # Retrieve Entrez Gene ID for all genes
        # Retrieve
        ID_all <- select(org.Hs.eg.db, keys=genes.universe, columns=c("ENTREZID", "SYMBOL"), keytype="SYMBOL")

        # Remove non-matches
        ID_all <- ID_all[which(!is.na(ID_all$ENTREZID)), "ENTREZID"]

    # Retrieve Entrez Gene ID for sig genes
        # Retrieve Entrez Gene ID for significant genes
        ID <- select(org.Hs.eg.db, keys=genes.sig, columns=c("ENTREZID", "SYMBOL"), keytype="SYMBOL")

        # Remove non-matches
        ID <- ID[which(!is.na(ID$ENTREZID)), "ENTREZID"]

    # Test of over-representation
        # Set up parameters for analysis
        length(ID) ; length(ID_all)
        #params <- new('GOHyperGParams', geneIds=ID, universeGeneIds=ID_all, annotation='org.Hs.eg.db', ontology="BP", pvalueCutoff=0.05, conditional=TRUE, testDirection="over")
        params <- new('GOHyperGParams', geneIds=ID, universeGeneIds=ID_all, annotation='org.Hs.eg.db', ontology="BP", pvalueCutoff=1.00, conditional=TRUE, testDirection="over")

        # Analyse
        go <- hyperGTest(params)

        # Generate result table
        go.table <- summary(go)

        # Filter significant terms after adjustment
        go.table$p.val.adj <- p.adjust(go.table$Pvalue, method=method.adjust)
        
        # Indicate no. of hits
        #go.table$Total_genes_sig <- length(unique(genes.sig))
        
        # Calculate precentage of hits
        #go.table$Count_pct <- round((go.table$Count/go.table$Total_genes_sig)*100, 2)
        
        # Reorder columns
        #go.table <- go.table[, c(1, 3:5, 10:11, 6:7, 2, 8:9)]
        
        # Subset pathways to annotate gene hits
        go.table <- go.table[which(go.table$p.val.adj < p.val.adj.return), ]

        if(nrow(go.table) != 0) {
            
            # Annotate with gene hits
                # Retrieve BP genes
                genes <- select(org.Hs.eg.db, keys=go.table$GOBPID, columns=c("GOALL", "SYMBOL"), keytype="GOALL")
                genes <- genes[which(genes$ONTOLOGYALL=="BP"), ]
                
                overlap <- NULL
                
                #pb <- txtProgressBar(1, length(go.table$GOBPID), style=3)
                
                # Retrieve overlap
                for(i in 1:length(go.table$GOBPID)) {
                
                    # Subset all genes belonging to GO term
                    genes.all <- unique(genes[which(genes$GOALL==go.table$GOBPID[i]), "SYMBOL"])
                    
                    # Find overlap
                    overlap[i] <- paste(intersect(genes.all, genes.sig), collapse="|")
                    
                    # Track progress
                    #setTxtProgressBar(pb, i)
                    
                }
                
                # Annotate
                go.table$genes.hit <- overlap
                
        # Save into new slot
        MarvelObject$DE$BioPathways <- go.table

        #####################################################################################
        
        # Subset top genes to plot
            if(nrow(go.table) > plot.top.n) {
                
                    go.table.small <- go.table[c(1:plot.top.n), ]
                
                } else {
                    
                    go.table.small <- go.table
                    
                }
            
            # transform adjusted p-val
            go.table.small$p.val.adj <- -log10(go.table.small$p.val.adj)
            
            # Subset relevant columns
            go.table.small <- go.table.small[,c("Term", "OddsRatio", "p.val.adj")]
            
            # Wrap label
            labels <- data.frame("x"=as.character(go.table.small$Term), stringsAsFactors=FALSE)
            labels$y <- ifelse(nchar(labels$x) < 25, 10, 20)
            labels$newx <- str_wrap(labels$x, width=25)
            go.table.small$Term <- labels$newx
            go.table.small$Term <- factor(go.table.small$Term, levels=rev(go.table.small$Term))
            
            # Dotplot
                # Definitions
                data <- go.table.small
                x <- data$Term
                y <- data$p.val.adj
                z <- data$OddsRatio
                maintitle <- ""
                xtitle <- ""
                ytitle <- "-log10(p-val)"
                #fivenum(y) ; ymin <- 3.5 ; ymax <- 5.5 ; yinterval <- 0.5
                legendtitle <- "Odds Ratio"
                
                # Plot
                plot <- ggplot() +
                    geom_point(data, mapping=aes(x=x, y=y, size=z)) +
                    #scale_y_continuous(breaks=seq(ymin, ymax, by=yinterval), limits=c(ymin, ymax)) +
                    labs(title=maintitle, x=xtitle, y=ytitle, size=legendtitle) +
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
                        axis.text.x=element_text(size=10, colour="black"),
                        axis.text.y=element_text(size=8, colour="black"),
                        legend.title=element_text(size=12),
                        legend.text=element_text(size=12)
                        ) +
                    coord_flip()
            
                # Save into new slot
                MarvelObject$DE$BioPathwaysPlot <- plot
            
        } else {
        
            MarvelObject$DE$BioPathways <- "No statistically significant pathways identified at user-specific adjusted p-value threshold"
            MarvelObject$DE$BioPathwaysPlot <- "No statistically significant pathways identified at user-specific adjusted p-value threshold"
    
        }
        
        
    } else {
        
        MarvelObject$DE$BioPathways <- "Number of differentially spliced genes less than user-specified threshold in min.genes argument"
        MarvelObject$DE$BioPathwaysPlot <- "Number of differentially spliced genes less than user-specified threshold in min.genes argument"
        
    }
        
    return(MarvelObject)
    
}
