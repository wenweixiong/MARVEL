#' @title Principle Component Analysis on Percent Spliced-in (PSI) Values
#'
#' @description
#' \code{RunPCA.PSI} performs principle component analysis on percent spliced-in (PSI) values.
#'
#' @details
#' This function performs principle component analysis on percent spliced-in (PSI) values and visualise cells on a reducted dimension space, i.e. 2D scatterplot.
#'
#' @param MarvelObject S3 object generated from \code{CreateMarvelObject} or \code{ComputePSI} function.
#' @param cell.types Character string. To indicate which groups of cells that will be used for analysis. Group names should match those in \code{cell.type} column of \code{$SplicePheno} slot.
#' @param n.cells Numeric value. The minimum no. of cells expressing the splicing event or gene for the event or gene, respectively, to be included for analysis.
#' @param features Character string. Vector of \code{tran_id}s for analysis. Should match \code{tran_id} column of \code{$ValidatedSpliceFeature} slot.
#' @param point.size Numeric value. Size of data points on reduced dimension space.
#' @param event.type Character string. Indicate which splicing event type to include for analysis. Can take value \code{"SE"}, \code{"MXE"}, \code{"RI"}, \code{"A5SS"}, or \code{"A3SS"} which represents skipped-exon (SE), mutually-exclusive exons (MXE), retained-intron (RI), alternative 5' splice site (A5SS), and alternative 3' splice site (A3SS), respectively.
#' @param seed Numeric value.  Ensures imputed values for NA PSIs are reproducible.
#' @export
#' @return An object of class S3 containing all the original slots as inputted by the user in addition to one new slot named \code{MarvelObject$PCA$PSI}. Contains both 2D scatterplot in \code{MarvelObject$PCA$PSI$Plot} and the corresponding x- and y-coordinates for each sample in \code{MarvelObject$PCA$PSI$Results}.
#' @author Sean Wen <sean.wenwx@gmail.com>
#' @importFrom plyr join
#' @import stats
#' @import methods
#' @import FactoMineR
#' @import factoextra
#' @import ggplot2
#' @examples
#' marvel <- readRDS(system.file("extdata/Data", "MarvelObject.rds", package="MARVEL"))
#'
#' features <- do.call(rbind.data.frame, marvel$SpliceFeatureValidated)
#' features <- features$tran_id
#'
#' marvel <- RunPCA.PSI(MarvelObject=marvel,
#'                      cell.types="all",
#'                      n.cells=3,
#'                      features=features,
#'                      point.size=2.5,
#'                      event.type="all",
#'                      seed=1
#'                      )
#'
#' marvel$PCA$PSI$Results
#' marvel$PCA$PSI$Plot


RunPCA.PSI <- function(MarvelObject, cell.types, n.cells, features, point.size, event.type, seed) {

    # Define arguments
    df <- do.call(rbind.data.frame, MarvelObject$PSI)
    df.pheno <- MarvelObject$SplicePheno
    df.feature <- do.call(rbind.data.frame, MarvelObject$SpliceFeatureValidated)
    cell.types <- cell.types
    n.cells <- n.cells
    features <- features
    seed <- seed
    event.type <- event.type
        
    # Create row names for matrix
    row.names(df) <- df$tran_id
    df$tran_id <- NULL
    
    # Subset overlapping samples in matrix and pheno file
    df <- df[, which(names(df) %in% df.pheno$sample.id)]
    
    # Subset relevant cell types
    if(cell.types[1] != "all") {
        
        df.pheno <- df.pheno[which(df.pheno$cell.type %in% cell.types), ]
        df <- df[, which(names(df) %in% df.pheno$sample.id)]
        
    }
    
    # Subset relevant event type
    if(event.type[1] != "all") {
        
        df.feature <- df.feature[which(df.feature$event_type %in% event.type), ]
        df <- df[df.feature$tran_id, ]
        
    }
    
    # Subset events with sufficient cells
    . <- apply(df, 1, function(x) {sum(!is.na(x))})
    index.keep <- which(. >= n.cells)
    df <- df[index.keep, ]
    df.feature <- df.feature[which(df.feature$tran_id %in% row.names(df)), ]
    
    # Subset features to reduce
    df.feature <- df.feature[which(df.feature$tran_id %in% features), ]
    df <- df[df.feature$tran_id, ]
                        
    # Check if matrix column and rows align with metadata
        # Column
        index.check <- which(unique((names(df)==df.pheno$sample.id))==FALSE)
        
        if(length(index.check)==0) {
            
            print("Checking... Matrix column (sample) names match sample metadata")
            
        } else {
            
            print("Checking... Matrix column (sample) names DO NOT match sample metadata")
            
        }
        
        # Row
        index.check <- which(unique((row.names(df)==df.feature$tran_id))==FALSE)
        
        if(length(index.check)==0) {
            
            print("Checking... Matrix row (feature) names match feature metadata")
            
        } else {
            
            print("Checking... Matrix row (feature) names DO NOT match feature metadata")
            
        }
    
    # Impute missing values
    set.seed(seed)
    df[is.na(df)] <- runif(n=sum(is.na(df)), min=0, max=1)

    # Reduce dimension
    res.pca <- PCA(as.data.frame(t(df)), scale.unit=TRUE, ncp=20, graph=FALSE)
        
    # Scatterplot
        # Definition
        data <- as.data.frame(res.pca$ind$coord)
        x <- data[,1]
        y <- data[,2]
        z <- df.pheno$cell.type
        maintitle <- paste(nrow(df), " splicing events", sep="")
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
    MarvelObject$PCA$PSI$Results <- res.pca
    MarvelObject$PCA$PSI$Plot <- plot
    
    return(MarvelObject)
        
}
