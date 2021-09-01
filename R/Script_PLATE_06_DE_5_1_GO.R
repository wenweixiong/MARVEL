#' @title Gene Ontology Analysis
#'
#' @description
#' \code{BioPathways} performs gene ontology analysis on genes that are differentially spliced.
#'
#' @details
#' This function performs gene ontology analysis on genes that are differentially spliced to identify significantly regulated gene sets (biological pathways).
#'
#' @param MarvelObject S3 object generated from \code{CompareValues} function.
#' @param de.p.val.adj Numeric value. Adjusted p-value below which the splicing events are considered differentially spliced and their corresponding genes are included for gene ontology analysis. If this argument is specified, then \code{n.top} must not be specified.
#' @param n.top Numeric value. Alternative to \code{de.p.val.adj}. Indicate the top n splicing events with the smallest adjusted p-values are differentially spliced and their corresponding genes are included for gene ontology analysis. If this argument is specified, then \code{de.p.val.adj} must not be specified.
#' @param min.gene.set.size Numeric value. Only evaluate enrichment of gene sets if the gene set is equal or bigger than this value. Default value is 10.
#' @param method.adjust Character string. Adjust p-values for multiple testing. Options available as per \code{p.adjust} function.
#' @param remove.ribo Logical. If set to \code{TRUE}, all ribosomal genes (genes with prefix RPS and RPL) meeting \code{de.p.val.adj} or \code{n.top} criteria will be removed prior to gene ontology analysis. This can prevent ribosomal genes from obscuring enriched gene sets that are not related to ribosomal genes.
#' @param annotate Logical. If set to \code{TRUE}, the genes differentially spliced that are found in each gene set are returned.
#'
#' @export
#'
#' @return An object of class S3 containing all the original slots as inputted by the user in addition to a new slot named \code{MarvelObject$DE$BioPathways} containing the gene sets significantly enriched among differentially spliced genes.
#'
#' @author Sean Wen <sean.wenwx@gmail.com>
#'
#' @importFrom plyr join
#' @import stats
#' @import methods
#' @import org.Hs.eg.db
#' @import GO.db
#' @import GOstats
#' @importFrom AnnotationDbi select
#'
#' @examples
#' # Load input
#' marvel <- readRDS(system.file("extdata/Data", "MarvelObject.rds", package="MARVEL"))
#'
#' # Run example
#' marvel <- BioPathways(MarvelObject=marvel,
#'                       n.top=5,
#'                       method.adjust="fdr"
#'                       )

BioPathways <- function(MarvelObject, de.p.val.adj=NULL, n.top=NULL, min.gene.set.size=10, method.adjust, remove.ribo=FALSE, annotate=FALSE) {
    
    # Define arguments
    de <- MarvelObject$DE$PSI$Table
    de.p.val.adj <- de.p.val.adj
    n.top <- n.top
    min.gene.set.size <- min.gene.set.size
    method.adjust <- method.adjust
    annotate <- annotate
    
    # Example arguments
    #de <- results.ad
    #de.p.val.adj <- 0.05
    #n.sig.top <- 10000
    #min.gene.set.size <- 10
    #method.adjust <- "fdr"
    #annotate <- TRUE
    
    # Define gene universe
    genes.universe <- unique(de$gene_short_name)

    # Subset significant genes
    if(!is.null(de.p.val.adj)) {
        
        genes.sig <- unique(de[which(de$p.val.adj < de.p.val.adj), "gene_short_name"])
        
    } else if(!is.null(n.top)) {
        
        genes.sig <- unique(de$gene_short_name[c(1:n.top)])
        
    }
    
    # Check if sufficient genes for GO analysis
    if(length(genes.sig) < 10) {
        
        return(print("Need to have at least 10 significant genes for GO analysis!"))
        
    }
 
    # Remove ribo genes
    if(remove.ribo==TRUE) {
        
        index.ribo.1 <- grep("^RPL", genes.sig)
        index.ribo.2 <- grep("^RPS", genes.sig)
        index.ribo <- unique(c(index.ribo.1, index.ribo.2))
        genes.sig <- genes.sig[-index.ribo]
        
    }
    
    # Report progress
    print(paste(length(genes.universe), " unique genes identified", sep=""))
    print(paste(length(genes.sig), " unique differentially spliced genes identified", sep=""))

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
        #params <- new('GOHyperGParams', geneIds=ID, universeGeneIds=ID_all, annotation='org.Hs.eg.db', ontology="BP", pvalueCutoff=0.05, conditional=TRUE, testDirection="over")
        params <- new('GOHyperGParams', geneIds=ID, universeGeneIds=ID_all, annotation='org.Hs.eg.db', ontology="BP", pvalueCutoff=1.00, conditional=FALSE, testDirection="over")

        # Analyse
        go <- hyperGTest(params)

        # Generate result table
        go.table <- summary(go)

        # Subset gene set size above minimum threshold
        go.table <- go.table[which(go.table$Size >= min.gene.set.size), ]
        
        # Adjust for multiple testing
        go.table$p.val.adj <- p.adjust(go.table$Pvalue, method=method.adjust)
        
        # Indicate no. of hits
        #go.table$Total_genes_sig <- length(unique(genes.sig))
        
        # Calculate precentage of hits
        #go.table$Count_pct <- round((go.table$Count/go.table$Total_genes_sig)*100, 2)
        
        # Reorder columns
        #go.table <- go.table[, c(1, 3:5, 10:11, 6:7, 2, 8:9)]
        
        # Subset pathways to annotate gene hits
        #go.table <- go.table[which(go.table$p.val.adj < p.val.adj.return), ]

    # Annotate with gene hits
    if(annotate==TRUE) {
        
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
            
    }
            
    # Save into new slot
    MarvelObject$DE$BioPathways$Table <- go.table
        
    return(MarvelObject)
    
}
