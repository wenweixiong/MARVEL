#' @title Principle component analysis for gene Data
#'
#' @description Performs principle component analysis using gene expression values.
#'
#' @param MarvelObject Marvel object. S3 object generated from \code{TransformExpValues} function.
#' @param cell.group.list List of character strings. Each element of the list is a vector of cell IDs corresponding to a cell group. The name of the element will be the cell group label.
#' @param min.cells Numeric value. The minimum no. of cells expressing the splicing event to be included for analysis.
#' @param features Character string. Vector of \code{gene_id} for analysis. Should match \code{gene_id} column of \code{MarvelObject$GeneFeature}.
#' @param point.size Numeric value. Size of data points on reduced dimension space.
#' @param point.alpha Numeric value. Transparency of the data points on reduced dimension space. Take any values between 0 to 1. The smaller the value, the more transparent the data points will be.
#' @param point.stroke Numeric value. The thickness of the outline of the data points. The larger the value, the thicker the outline of the data points.
#' @param retrieve.non.outliers Logical. If set to \code{TRUE}, this function will retrieve \code{sample.id} of non-outliers based on the intial PCA. Define the non-outliers based on the initial PCA coordinates. Use in conjunction with arguments \code{pc1.min}, \code{pc1.max}, \code{pc2.min}, and \code{pc2.max}.
#' @param pc1.min Numeric value. When \code{retrieve.non.outliers} set to \code{TRUE}. To indicate the PC1 value above which to retrieve the sample IDs. Use in conjunction with \code{pc1.max}, \code{pc2.min}, and \code{pc2.max}.
#' @param pc1.max Numeric value. When \code{retrieve.non.outliers} set to \code{TRUE}. To indicate the PC1 value below which to retrieve the sample IDs. Use in conjunction with \code{pc1.min}, \code{pc2.min}, and \code{pc2.max}.
#' @param pc2.min Numeric value. When \code{retrieve.non.outliers} set to \code{TRUE}. To indicate the PC2 value above which to retrieve the sample IDs. Use in conjunction with \code{pc1.max}, \code{pc1.max}, and \code{pc2.max}.
#' @param pc2.max Numeric value. When \code{retrieve.non.outliers} set to \code{TRUE}. To indicate the PC2 value below which to retrieve the sample IDs. Use in conjunction with \code{pc1.min}, \code{pc2.min}, and \code{pc2.min}.
#' @param remove.outliers Logical. If set to \code{TRUE}, re-run PCA by only including non-outliers. Use after running the function with \code{retrieve.non.outliers} set to \code{TRUE}.
#' @param cell.group.colors Character string. Vector of colors for the cell groups specified for PCA analysis using \code{cell.type.columns}, \code{cell.type.variable}, and \code{cell.type.labels}. If not specified, default \code{ggplot2} colors will be used.
#'
#' @return An object of class S3 containing with new slots \code{MarvelObject$PCA$Exp$Results}, \code{MarvelObject$PCA$Exp$Plot}, and \code{MarvelObject$PCA$Exp$Plot.Elbow}.
#'
#' @importFrom plyr join
#' @import stats
#' @import methods
#' @import FactoMineR
#' @import factoextra
#' @import ggplot2
#' @importFrom grDevices hcl
#'
#' @export


RunPCA.Exp <- function(MarvelObject, cell.group.list,
                       min.cells=25, features, point.size=0.5, point.alpha=0.75, point.stroke=0.1,
                       retrieve.non.outliers=FALSE, pc1.min=NULL, pc1.max=NULL, pc2.min=NULL, pc2.max=NULL,
                       remove.outliers=FALSE, cell.group.colors=NULL
                       ) {

    # Define arguments
    df <- MarvelObject$Exp
    df.pheno <- MarvelObject$SplicePheno
    df.feature <- MarvelObject$GeneFeature
    cell.group.list <- cell.group.list
    min.cells <- min.cells
    features <- features
    point.size <- point.size
    point.alpha <- point.alpha
    point.stroke <- point.stroke
    retrieve.non.outliers <- retrieve.non.outliers
    pc1.min <- pc1.min
    pc1.max <- pc1.max
    pc2.min <- pc2.min
    pc2.max <- pc2.max
    remove.outliers <- remove.outliers
    cell.group.colors <- cell.group.colors
    
    # Example arguments
    #MarvelObject <- marvel
    #df <- MarvelObject$Exp
    #df.pheno <-  MarvelObject$SplicePheno
    #df.feature <-  MarvelObject$GeneFeature
    #cell.group.list <- cell.group.list
    #min.cells <- 25
    #features <- gene_ids
    #point.size <- 0.5
    #cell.group.colors <- NULL
    #remove.outliers <- FALSE
    #retrieve.non.outliers <- FALSE
    #point.alpha <- 0.75
    #point.stroke <- 0.1
    #point.size <- 0.5
    
    ######################################################################
    
    # Create row names for matrix
    row.names(df) <- df$gene_id
    df$gene_id <- NULL
    
    # Retrieve sample ids for each group
    .list <- list()
    
    for(i in 1:length(cell.group.list)) {
        
        
        .list[[i]] <- data.frame("sample.id"=cell.group.list[[i]],
                                 "cell.type.label"=names(cell.group.list)[[i]],
                                 stringsAsFactors=FALSE
                                 )
        
    }
    
    md <- do.call(rbind.data.frame, .list)
       
    # Set factor levels
    md$cell.type.label <- factor(md$cell.type.label, levels=names(cell.group.list))
    
    # Subset relevant cells
    df.pheno <- df.pheno[which(df.pheno$sample.id %in% md$sample.id), ]
    df.pheno <- join(df.pheno, md, by="sample.id", type="left")
    df <- df[, which(names(df) %in% df.pheno$sample.id)]
    
    # Subset genes with sufficient cells
    . <- apply(df, 1, function(x) {sum(x > 0)})
    index.keep <- which(. >= min.cells)
    df <- df[index.keep, ]
    df.feature <- df.feature[which(df.feature$gene_id %in% row.names(df)), ]
    
    # Subset features to reduce
    df.feature <- df.feature[which(df.feature$gene_id %in% features), ]
    df <- df[df.feature$gene_id, ]
                    
    # Remove outliers (previously identified)
    if(remove.outliers==TRUE) {
        
        df.pheno <- df.pheno[which(df.pheno$sample.id %in% MarvelObject$PCA$Exp$sample.ids.non.outliers), ]
        df <- df[, df.pheno$sample.id]
        
    }
    
    # Reduce dimension
    res.pca <- PCA(as.data.frame(t(df)), scale.unit=TRUE, ncp=20, graph=FALSE)
    
    # Retrieve sample IDs of non-outliers (to be used for 2nd round of PCA)
    if(retrieve.non.outliers==TRUE) {
        
        # Retrieve coordinates
        . <- as.data.frame(res.pca$ind$coord)
        
        # Subset PC1
        . <- .[which(.$Dim.1 > pc1.min & .$Dim.1 < pc1.max), ]
        
        # Subset PC2
        . <- .[which(.$Dim.2 > pc2.min & .$Dim.2 < pc2.max), ]
        
        # Save non-outliers
        MarvelObject$PCA$Exp$sample.ids.non.outliers <- row.names(.)
        
    }
    
    # Scatterplot
        # Definition
        data <- as.data.frame(res.pca$ind$coord)
        x <- data[,1]
        y <- data[,2]
        z <- df.pheno$cell.type.label
        maintitle <- paste(nrow(df), " genes", sep="")
        xtitle <- paste("PC1 (", round(get_eigenvalue(res.pca)[1,2], digits=1), "%)" ,sep="")
        ytitle <- paste("PC2 (", round(get_eigenvalue(res.pca)[2,2], digits=1), "%)" ,sep="")
        legendtitle <- "Group"
        
        # Color scheme
        if(is.null(cell.group.colors[1])) {
        
            gg_color_hue <- function(n) {
              hues = seq(15, 375, length = n + 1)
              hcl(h = hues, l = 65, c = 100)[1:n]
            }
            
            n = length(levels(z))
            cols = gg_color_hue(n)
        
        } else {
            
            cols <- cell.group.colors
            
        }
        
        # Plot
        plot <- ggplot() +
            geom_point(data, mapping=aes(x=x, y=y, fill=z), size=point.size, pch=21, alpha=point.alpha, stroke=point.stroke) +
            scale_fill_manual(values=cols) +
            labs(title=maintitle, x=xtitle, y=ytitle, fill=legendtitle) +
            theme(panel.grid.major=element_blank(),
                panel.grid.minor=element_blank(),
                panel.background=element_blank(),
                plot.title = element_text(size=12, hjust=0.5),
                axis.line=element_line(colour = "black"),
                axis.title=element_text(size=12),
                axis.text.x=element_text(size=10, colour="black"),
                axis.text.y=element_text(size=10, colour="black"),
                legend.title=element_text(size=8),
                legend.text=element_text(size=8)
                )
    
    # Elbow plot
        # Definition
        if(nrow(get_eigenvalue(res.pca)) >= 20) {
            
            data.2 <- as.data.frame(get_eigenvalue(res.pca)[c(1:20),])
            
        } else {
            
            data.2 <- as.data.frame(get_eigenvalue(res.pca)[c(1:nrow(get_eigenvalue(res.pca))),])
        
        }
        
        x.2 <- c(1:nrow(data.2))
        y.2 <- data.2$variance.percent
        maintitle.2 <- ""
        xtitle.2 <- "PC"
        ytitle.2 <- "Variance Explained (%)"
        
        # Plot
        plot.2 <- ggplot() +
            geom_point(data.2, mapping=aes(x=x.2, y=y.2), size=1, pch=21) +
            labs(title=maintitle.2, x=xtitle.2, y=ytitle.2) +
            theme(panel.grid.major=element_blank(),
                panel.grid.minor=element_blank(),
                panel.background=element_blank(),
                plot.title = element_text(size=12, hjust=0.5),
                axis.line=element_line(colour = "black"),
                axis.title=element_text(size=12),
                axis.text=element_text(size=10, colour="black"),
                )
                          
    # Save to new slot
    MarvelObject$PCA$Exp$Results <- res.pca
    MarvelObject$PCA$Exp$Plot <- plot
    MarvelObject$PCA$Exp$Plot.Elbow <- plot.2
    
    return(MarvelObject)
        
}
