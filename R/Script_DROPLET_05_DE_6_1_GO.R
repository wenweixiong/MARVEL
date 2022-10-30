#' @title Pathway enrichment analysis
#'
#' @description Performs pathway enrichment analysis on differentially spliced genes or user-specified custom set of genes.
#'
#' @param MarvelObject Marvel object. S3 object generated from \code{CompareValues.Genes.10x} function.
#' @param method Character string. The statistical method used for differential splicing analysis.
#' @param pval Numeric value. p-value, above which, the splice junction is considered differentially spliced. Default is \code{0.05}.
#' @param log2fc Numeric value. Absolute log2 fold change from differential splicing analysis, above which, the splice junction is considered differentially spliced. This option should be \code{NULL} if \code{delta} has been specified.
#' @param delta Numeric value. Absolute difference in average PSI values between the two cell groups, above which, the splice junction is considered differentially spliced. This option should be \code{NULL} if \code{log2fc} has been specified.
#' @param min.gene.norm Numeric value. The average normalised gene expression across the two cell groups above which the splice junction is considered differentially spliced. Default is \code{0}.
#' @param method.adjust Character string. Adjust p-values for multiple testing. Options available as per \code{p.adjust} function.
#' @param custom.genes Character strings. Alternative to \code{pval} and \code{delta}. Vector of gene names to be assessed for enrichment of biological pathways.
#' @param species Character strings. Takes the value \code{"human"} or \code{"mouse"}, which corresponds to human and mouse genes, respectively. Default value is \code{"human"}.
#' @param remove.ribo Logical value. If set to \code{TRUE}, ribosomal genes will be removed prior to GO analysis. This may prevent high-expressing ribosomal genes from overshadowing more biological relevant genes for GO analysis. Default value is \code{FALSE}.
#'
#' @return An object of class S3 containing new slot \code{MarvelObject$DE$BioPathways$Table}.
#'
#' @importFrom plyr join
#' @importFrom stats p.adjust p.adjust.methods
#' @import methods
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
#' marvel.demo.10x <- BioPathways.10x(
#'                         MarvelObject=marvel.demo.10x,
#'                         custom.genes=c("TPM2", "GNAS"),
#'                         species="human"
#'                         )

BioPathways.10x <- function(MarvelObject, pval=0.05, log2fc=NULL, delta=5, min.gene.norm=0, method.adjust="fdr", custom.genes=NULL, species="human", remove.ribo=FALSE) {
    
    # Define arguments
    MarvelObject <- MarvelObject
    df <- MarvelObject$DE$SJ$Table
    pval <- pval
    delta <- delta
    log2fc <- log2fc
    method.adjust <- method.adjust
    custom.genes <- custom.genes
    species <- species
    min.gene.norm <- min.gene.norm
    remove.ribo <- remove.ribo
    
    # Example arguments
    #MarvelObject <- marvel
    #df <- MarvelObject$DE$SJ$Table
    #pval <- 0.05
    #delta <- 10
    #log2fc <- NULL
    #method.adjust <- "fdr"
    #custom.genes <- gene_short_names
    #species <- "human"
    #min.gene.norm <- 1
    #remove.ribo <- TRUE
    
    #############################################################
    
    if(is.null(custom.genes[1])) {
        
        # Subset expressed genes
        df <- df[which(df$mean.expr.gene.norm.g1.g2 > min.gene.norm), ]
            
        # Indicate sig events and direction
        if(!is.null(log2fc)) {
            
            df$sig <- NA
            df$sig[which(df$pval < pval & df$log2fc > log2fc)] <- "up"
            df$sig[which(df$pval < pval & df$log2fc < (log2fc*-1))] <- "down"
            df$sig[is.na(df$sig)] <- "n.s."
            df$sig <- factor(df$sig, levels=c("up", "down", "n.s."))
            table(df$sig)
            
        } else if(!is.null(delta)){
            
            df$sig <- NA
            df$sig[which(df$pval < pval & df$delta > delta)] <- "up"
            df$sig[which(df$pval < pval & df$delta < (delta*-1))] <- "down"
            df$sig[is.na(df$sig)] <- "n.s."
            df$sig <- factor(df$sig, levels=c("up", "down", "n.s."))
            table(df$sig)
            
        }
        
        # Subset genes with DE SJ
        df <- df[which(df$sig %in% c("up", "down")), ]
        
        gene_short_names <- unique(df$gene_short_name)
        
    } else {
    
        gene_short_names <- custom.genes
        
    }
    
    # Remove ribo genes
    gene_short_names <- grep("^RP[S/L]", gene_short_names, value=TRUE, invert=TRUE)
    gene_short_names <- grep("^rp[s/l]", gene_short_names, value=TRUE, invert=TRUE)

    # Print progress
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

    } else {

        message("No significant GO terms (pathways) found")

    }
    
    # Compute pathway enrichment
    if(nrow(ego2) != 0) {

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
        
    }
    
    # Save into new slot
    if(nrow(ego) != 0) {
        
        MarvelObject$DE$BioPathways$Table <- ego2
        
    } else {
        
        MarvelObject$DE$BioPathways$Table <- NULL
        
    }
        
    return(MarvelObject)
    
}
