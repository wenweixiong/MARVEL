#' @title Create pseudobulk for downstream differential splicing analysis
#'
#' @description Create pseudobulk by merging gene counts and SJ counts by donor for each cell type. The results in each cell type from each donor having aggregated gene counts and SJ counts.
#'
#' @param MarvelObject Marvel object. S3 object generated from \code{CheckAlignment.10x} function.
#' @param donor.id.column Character strings. Column name of the sample metadata the corresponds to the donor ID for pseudo-bulking.
#' @param feature.column Character strings. Column name of the sample metadata the corresponds to the cell type for pseudo-bulking.
#' @param comparison.column Character strings. Column name of the sample metadata the corresponds to the cell groups for downstream pair-wise comparison.
#' @param min.pct.cells.genes Numeric value. Minimum percentage of cells in which the gene is expressed for that gene to be included for pseudobulk-ing. Expressed genes defined as genes with non-zero normalised UMI counts.
#'
#' @return An object of class S3 with a new slots \code{MarvelObject$DE$PseudoBulk$GeneCountMatrix}, \code{MarvelObject$DE$PseudoBulk$SJCountMatrix}, and \code{MarvelObject$DE$PseudoBulk$sample.metadata}.
#'
#' @importFrom plyr join rbind.fill
#' @import utils
#' @import Matrix
#'
#' @export

CreatePseudoBulk.SJ.10x <- function(MarvelObject, donor.id.column, feature.column, comparison.column, min.pct.cells.genes=10) {
        
    # Define arguments
    MarvelObject <- MarvelObject
    sample.metadata <- MarvelObject$sample.metadata
    sj.metadata <- MarvelObject$sj.metadata
    df.gene.norm <- MarvelObject$gene.norm.matrix
    df.gene.count <- MarvelObject$gene.count.matrix
    df.sj.count <- MarvelObject$sj.count.matrix
    donor.id.column <- donor.id.column
    feature.column <- feature.column
    comparison.column <- comparison.column
    min.pct.cells.genes <- min.pct.cells.genes
    
    # Example arguments
    #MarvelObject <- marvel
    #sample.metadata <- MarvelObject$sample.metadata
    #sj.metadata <- MarvelObject$sj.metadata
    #df.gene.norm <- MarvelObject$gene.norm.matrix
    #df.gene.count <- MarvelObject$gene.count.matrix
    #df.sj.count <- MarvelObject$sj.count.matrix
    #donor.id.column <- "individual"
    #feature.column <- "cluster"
    #comparison.column <- "diagnosis"
    #min.pct.cells.genes <- 5
    
    #################################################################
    
    # Define donor ids
    donor.ids <- unique(sample.metadata[[donor.id.column]])[c(1:2)]
    
    .list.gene.count.collapsed.2 <- list()
    .list.sj.count.collapsed.2 <- list()
    
    for(j in 1:length(donor.ids)) {
        
        # Define donor id
        donor.id <- donor.ids[j]
        
        # Subset
        sample.metadata.small <- sample.metadata[which(sample.metadata[[donor.id.column]]==donor.id), ]
        
        # Define cell types
        features <- unique(sample.metadata[[feature.column]])[c(1:2)]
        
        .list.gene.count.collapsed <- list()
        .list.sj.count.collapsed <- list()
        
        for(i in 1:length(features)) {
            
            # Define cell type
            feature <- features[i]
            
            # Subset
            sample.metadata.small. <- sample.metadata.small[which(sample.metadata.small[[feature.column]]==feature), ]
            
            # Retrieve cell ids
            cell.ids <- sample.metadata.small.$cell.id
            
            # Retrieve expressed genes
            df.gene.norm.small <- df.gene.norm[, cell.ids]
            n.cells <- apply(df.gene.norm.small, 1, function(x) {sum(x != 0)})
            pct.cells <- n.cells/ncol(df.gene.norm.small) * 100
            gene_short_names <- names(pct.cells)[which(pct.cells > min.pct.cells.genes)]
            
            # Collapse gene count matrix
            df.gene.count.small <- df.gene.count[gene_short_names, cell.ids]
            . <- Matrix::rowSums(df.gene.count.small)
            . <- as.data.frame(.)
            names(.) <- paste(donor.id, feature, sep="@")
            df.gene.count.small.collapsed <- .
            
            # Retrieve expressed SJ
            #df.sj.count.small <- df.sj.count[, cell.ids]
            #n.cells <- apply(df.sj.count.small, 1, function(x) {sum(x != 0)})
            #pct.cells <- n.cells/ncol(df.sj.count.small) * 100
            #coord.introns <- names(pct.cells)[which(pct.cells > min.pct.cells.sj)]
            
            # Collapse sj count matrix
            #df.sj.count.small <- df.sj.count[coord.introns, cell.ids]
            df.sj.count.small <- df.sj.count[, cell.ids]
            . <- Matrix::rowSums(df.sj.count.small)
            . <- as.data.frame(.)
            names(.) <- paste(donor.id, feature, sep="@")
            df.sj.count.small.collapsed <- .
            
            # Save into list
            .list.gene.count.collapsed[[i]] <- df.gene.count.small.collapsed
            .list.sj.count.collapsed[[i]] <- df.sj.count.small.collapsed

        }
        
        # Merge all cell types
            # Gene counts
            . <- lapply(.list.gene.count.collapsed, function(x) {as.data.frame(t(x))})
            sample.ids <- unlist(lapply(., function(x) {row.names(x)}))
            . <- do.call(plyr::rbind.fill, .)
            row.names(.) <- sample.ids
            . <- as.data.frame(t(.))
            .list.gene.count.collapsed.2[[j]] <- .
            
            # SJ counts
            . <- lapply(.list.sj.count.collapsed, function(x) {as.data.frame(t(x))})
            sample.ids <- unlist(lapply(., function(x) {row.names(x)}))
            . <- do.call(plyr::rbind.fill, .)
            row.names(.) <- sample.ids
            . <- as.data.frame(t(.))
            .list.sj.count.collapsed.2[[j]] <- .

        # Track progress
        message(paste("donor ", donor.id, " collapsed", sep=""))
            
    }
    
    # Merge all donors
        # Gene counts
        . <- lapply(.list.gene.count.collapsed.2, function(x) {as.data.frame(t(x))})
        sample.ids <- unlist(lapply(., function(x) {row.names(x)}))
        . <- do.call(plyr::rbind.fill, .)
        row.names(.) <- sample.ids
        . <- as.data.frame(t(.))
        df.gene.count.collapsed <- .
        
        # SJ counts
        . <- lapply(.list.sj.count.collapsed.2, function(x) {as.data.frame(t(x))})
        sample.ids <- unlist(lapply(., function(x) {row.names(x)}))
        . <- do.call(plyr::rbind.fill, .)
        row.names(.) <- sample.ids
        . <- as.data.frame(t(.))
        df.sj.count.collapsed <- .
    
    # Collapse metadata
        # Rename sample.metadata object
        df.pheno <- sample.metadata
        
        # Collapse
        df.pheno.collapsed <- unique(df.pheno[,c(donor.id.column, feature.column, comparison.column)])
        
        # Create unique donor.id-feature id
        df.pheno.collapsed$donor.id_feature <- paste(df.pheno.collapsed[,donor.id.column], df.pheno.collapsed[,feature.column], sep="@")
        
        # Subset overlapping donor.id-feature id
        overlap <- intersect(names(df.gene.count.collapsed), df.pheno.collapsed$donor.id_feature)
        df.pheno.collapsed <- df.pheno.collapsed[which(df.pheno.collapsed$donor.id_feature %in% overlap), ]
        
        # Align matrix w/ metadata
        df.gene.count.collapsed <- df.gene.count.collapsed[, df.pheno.collapsed$donor.id_feature]
        df.sj.count.collapsed <- df.sj.count.collapsed[, df.pheno.collapsed$donor.id_feature]
        
        # Check alignment: gene matrix vs metadata
        index.l <- table(df.pheno.collapsed$donor.id_feature==colnames(df.gene.count.collapsed))
        index.true <- length(which(names(index.l)==TRUE))
        index.false <- length(which(names(index.l)==FALSE))
         
        if(index.true==1 & index.false==0) {
        
            message("Collapsed metadata and collapsed gene count matrix column names MATCHED")
            
        } else {
            
            message("Collapsed metadata and collapsed gene count matrix column names NOT MATCHED")

        }
    
        # Check alignment: sj matrix vs metadata
        index.l <- table(df.pheno.collapsed$donor.id_feature==colnames(df.sj.count.collapsed))
        index.true <- length(which(names(index.l)==TRUE))
        index.false <- length(which(names(index.l)==FALSE))
         
        if(index.true==1 & index.false==0) {
        
            message("Collapsed metadata and collapsed SJ count matrix column names MATCHED")
            
        } else {
            
            message("Collapsed metadata and collapsed SJ count matrix column names NOT MATCHED")

        }
    
    ################################################################
     
    # Save into new slot
    MarvelObject$DE$PseudoBulk$GeneCountMatrix <- df.gene.count.collapsed
    MarvelObject$DE$PseudoBulk$SJCountMatrix <- df.sj.count.collapsed
    MarvelObject$DE$PseudoBulk$sample.metadata <- df.pheno.collapsed
    
    # Return final object
    return(MarvelObject)
            
}


