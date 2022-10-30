#' @title Pathway enrichment analysis
#'
#' @description Performs pathway enrichment analysis on differentially spliced genes or user-specified custom set of genes.
#'
#' @param MarvelObject Marvel object. S3 object generated from \code{CompareValues} function.
#' @param method Character string. The statistical method used for differential splicing analysis.
#' @param pval Numeric value. Alternative to \code{n.top} and \code{custom.genes}, i.e. choose one of these three options. Adjusted p-value below which the splicing events are considered differentially spliced and their corresponding genes are included for gene ontology analysis. If this argument is specified, then \code{n.top} must not be specified.
#' @param delta Numeric value. The absolute difference between the means PSI values of cell group 1 and 2, above which, the splicing event is considered differentially spliced and their corresponding genes are included for gene ontology analysis.
#' @param n.top Numeric value. Alternative to \code{pval} to \code{custom.genes}, i.e. choose one of these three options.. Indicate the top n splicing events with the smallest adjusted p-values are differentially spliced and their corresponding genes are included for gene ontology analysis. If this argument is specified, then \code{pval} must not be specified.
#' @param method.adjust Character string. Adjust p-values for multiple testing. Options available as per \code{p.adjust} function.
#' @param custom.genes Character strings. Alternative to \code{pval} and \code{n.top}, i.e. choose one of these three options.. Vector of gene names to be assessed for enrichment of biological pathways.
#' @param species Character strings. Takes the value \code{"human"} or \code{"mouse"}, which corresponds to human and mouse genes, respectively. Default value is \code{"human"}. This will enable \code{MARVEL} to retrieve the relevant database for GO analysis.
#'
#' @return An object of class S3 with new slot \code{MarvelObject$DE$BioPathways$Table}.
#'
#' @importFrom plyr join
#' @importFrom stats p.adjust p.adjust.methods
#' @import methods
#'
#' @export
#'
#' @examples
#' marvel.demo <- readRDS(system.file("extdata/data", "marvel.demo.rds", package="MARVEL"))
#'
#' marvel.demo <- BioPathways(MarvelObject=marvel.demo,
#'                            method="ad",
#'                            custom.genes=c("RPL26", "SNRPN")
#'                            )

BioPathways <- function(MarvelObject, method=NULL, pval=NULL, delta=0, n.top=NULL, method.adjust="fdr", custom.genes=NULL, species="human") {
    
    # Define arguments
    method <- method
    pval <- pval
    delta <- delta
    n.top <- n.top
    method.adjust <- method.adjust
    custom.genes <- custom.genes
    species <- species
    
    # Example arguments
    #MarvelObject <- marvel
    #method <- c("ad", "dts")
    #pval <- 0.10
    #delta <- delta
    #n.top <- NULL
    #method.adjust <- "fdr"
    #custom.genes <- NULL
    #species <- "human"
    
    if(is.null(custom.genes[1])) {
            
        gene_short_names.list <- list()
        
        for(i in 1:length(method)) {
            
            # Retrieve DE result table
            de <- MarvelObject$DE$PSI$Table[[method[i]]]
            
            # Define sig genes
            if(!is.null(pval)) {
                
                index <- which(de$p.val.adj < pval & abs(de$mean.diff) > delta & de$outlier==FALSE)
                gene_short_names <- de[index, "gene_short_name"]
                
            } else if(!is.null(n.top)) {
                
                gene_short_names <- de[c(1:n.top), "gene_short_name"]
                
            }
            
            gene_short_names.list[[i]] <- unique(gene_short_names)
            
        }
        
        gene_short_names <- unique(unlist(gene_short_names.list))
        
    } else {
    
        gene_short_names <- custom.genes
        
    }
    
    message(paste(length(gene_short_names), " unique genes identified for GO analysis", sep=""))
    
    # Check if sufficient no. of genes for GO analysis
    if(length(gene_short_names) < 10) {
        
        message("Require mininum 10 genes for GO analysis")
        
        return(MarvelObject)
        
    }
    
    # Retrieve entrez IDs
    if(species=="human") {
        
        ID <- AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db, keys=gene_short_names, columns=c("ENTREZID", "SYMBOL"), keytype="SYMBOL")
        
    } else if(species=="mouse"){
        
        ID <- AnnotationDbi::select(org.Mm.eg.db::org.Mm.eg.db, keys=gene_short_names, columns=c("ENTREZID", "SYMBOL"), keytype="SYMBOL")
        
    }

    # GO analysis
    if(species=="human") {
        
        ego <- clusterProfiler::enrichGO(ID$ENTREZID, OrgDb = "org.Hs.eg.db", ont="BP", readable=TRUE, pAdjustMethod=method.adjust, pvalueCutoff=0.10)
        head(ego)
        
    } else if(species=="mouse"){
        
        ego <- clusterProfiler::enrichGO(ID$ENTREZID, OrgDb = "org.Mm.eg.db", ont="BP", readable=TRUE, pAdjustMethod=method.adjust, pvalueCutoff=0.10)
        head(ego)
        
    }
 
    # GO analysis (Remove redundant terms)
    if(nrow(ego) != 0) {

        ego2 <- clusterProfiler::simplify(ego, cutoff=0.7, by="p.adjust", select_fun=min)
        ego2 <- as.data.frame(ego2)
        
        # Compute gene ratio
        . <- strsplit(ego2$GeneRatio, split="/", fixed=TRUE)
        numerator <- as.numeric(sapply(., function(x) {x[1]}))
        denominator <- as.numeric(sapply(., function(x) {x[2]}))
        GeneRatio <- numerator/denominator
        
        # Compute background ratio
        . <- strsplit(ego2$BgRatio, split="/", fixed=TRUE)
        numerator <- as.numeric(sapply(., function(x) {x[1]}))
        denominator <- as.numeric(sapply(., function(x) {x[2]}))
        BgRatio <- numerator/denominator
        
        # Compute enrichment
        ego2$enrichment <- GeneRatio/BgRatio
        
        # Reorder columns
        cols.1 <- c("ID", "Description", "GeneRatio", "BgRatio", "enrichment")
        cols.2 <- names(ego2)[-which(names(ego2) %in% cols.1)]
        ego2 <- ego2[, c(cols.1, cols.2)]
        
        # Save into new slot
        MarvelObject$DE$BioPathways$Table <- ego2
        
    } else {

        message("No significant GO terms (pathways) found")
        
        # Save into new slot
        MarvelObject$DE$BioPathways$Table <- NULL

    }
            
    return(MarvelObject)
    
}
