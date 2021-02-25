#' @title Principle Component Analysis on Gene Expression Values
#'
#' @description
#' \code{RunPCA.Exp} performs principle component analysis on gene expression values.
#'
#' @details
#' This function performs principle component analysis on gene expression values and visualise cells on a reducted dimension space, i.e. 2D scatterplot.
#'
#' @param MarvelObject S3 object generated from \code{CreateMarvelObject} or \code{ComputePSI} function.
#' @param cell.types Character string. To indicate which groups of cells that will be used for analysis. Group names should match those in \code{cell.type} column of \code{$GenePheno} slot.
#' @param n.cells Numeric value. The minimum no. of cells expressing the splicing event or gene for the event or gene, respectively, to be included for analysis.
#' @param features Character string. Vector of \code{gene_id}s for analysis. Should match \code{gene_id} column of \code{$GeneFeature} slot.
#' @param point.size Numeric value. Size of data points on reduced dimension space.
#' @export
#' @return An object of class S3 containing all the original slots as inputted by the user in addition to one new slot named \code{MarvelObject$PCA$Gene}. Contains both 2D scatterplot in \code{MarvelObject$PCA$Gene$Plot} and the corresponding x- and y-coordinates for each sample in \code{MarvelObject$PCA$Gene$Results}.
#' @author Sean Wen <sean.wenwx@gmail.com>
#' @importFrom plyr join
#' @import stats
#' @import methods
#' @import FactoMineR
#' @import factoextra
#' @import ggplot2
#' @examples
#'
#' marvel <- readRDS(system.file("extdata/Data", "MarvelObject.rds", package="MARVEL"))
#'
#' features <- marvel$GeneFeature$gene_id
#'
#' marvel <- RunPCA.Exp(MarvelObject=marvel,
#'                      cell.types="all",
#'                      n.cells=3,
#'                      features=features,
#'                      point.size=2.5
#'                      )
#'
#' marvel$PCA$Gene$Results
#' marvel$PCA$Gene$Plot

RunPCA.Exp <- function(MarvelObject, cell.types, n.cells, features, point.size) {

    # Define arguments
    df <- MarvelObject$Exp
    df.pheno <- MarvelObject$GenePheno
    df.feature <- MarvelObject$GeneFeature
    cell.types <- cell.types
    n.cells <- n.cells
    features <- features
    
    # Create row names for matrix
    row.names(df) <- df$gene_id
    df$gene_id <- NULL
    
    # Subset overlapping samples in matrix and pheno file
    df <- df[, which(names(df) %in% df.pheno$sample.id)]
    
    # Subset relevant cell types
    if(cell.types[1] != "all") {
        
        df.pheno <- df.pheno[which(df.pheno$cell.type %in% cell.types), ]
        df <- df[, which(names(df) %in% df.pheno$sample.id)]
        
    }
    
    # Subset genes with sufficient cells
    . <- apply(df, 1, function(x) {sum(x > 0)})
    index.keep <- which(. >= n.cells)
    df <- df[index.keep, ]
    df.feature <- df.feature[which(df.feature$gene_id %in% row.names(df)), ]
    
    # Subset features to reduce
    df.feature <- df.feature[which(df.feature$gene_id %in% features), ]
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
        index.check <- which(unique((row.names(df)==df.feature$gene_id))==FALSE)
        
        if(length(index.check)==0) {
            
            print("Checking... Matrix row (feature) names match feature metadata")
            
        } else {
            
            print("Checking... Matrix row (feature) names DO NOT match feature metadata")
            
        }
        
    # Reduce dimension
    res.pca <- PCA(as.data.frame(t(df)), scale.unit=TRUE, ncp=20, graph=FALSE)
        
    # Scatterplot
        # Definition
        data <- as.data.frame(res.pca$ind$coord)
        x <- data[,1]
        y <- data[,2]
        z <- df.pheno$cell.type
        maintitle <- paste(nrow(df), " genes", sep="")
        xtitle <- paste("PC1 (", round(get_eigenvalue(res.pca)[1,2], digits=1), "%)" ,sep="")
        ytitle <- paste("PC2 (", round(get_eigenvalue(res.pca)[2,2], digits=1), "%)" ,sep="")
        legendtitle <- "Group"
        
        # Plot
        plot <- ggplot() +
            geom_point(data, mapping=aes(x=x, y=y, fill=z), size=point.size, pch=21) +
            labs(title=maintitle, x=xtitle, y=ytitle, fill=legendtitle) +
            theme(panel.grid.major=element_blank(),
                panel.grid.minor=element_blank(),
                panel.background=element_blank(),
                plot.title = element_text(size=12, hjust=0.5),
                axis.line=element_line(colour = "black"),
                axis.title=element_text(size=12),
                axis.text=element_text(size=10, colour="black"),
                legend.title=element_text(size=8),
                legend.text=element_text(size=8)
                )
                          
    # Save to new slot
    MarvelObject$PCA$Gene$Results <- res.pca
    MarvelObject$PCA$Gene$Plot <- plot
    
    return(MarvelObject)
        
}
