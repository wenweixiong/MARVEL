#' @title Principle component analysis for gene Data
#'
#' @description Performs principle component analysis using gene expression values.
#'
#' @param MarvelObject Marvel object. S3 object generated from \code{ComputePSI} function.
#' @param sample.ids Character strings. Specific cells to plot.
#' @param cell.group.column Character string. The name of the sample metadata column in which the variables will be used to label the cell groups on the PCA.
#' @param cell.group.order Character string. The order of the variables under the sample metadata column specified in \code{cell.group.column} to appear in the PCA cell group legend.
#' @param cell.group.colors Character string. Vector of colors for the cell groups specified for PCA analysis using \code{cell.type.columns} and \code{cell.group.order}. If not specified, default \code{ggplot2} colors will be used.
#' @param min.cells Numeric value. The minimum no. of cells expressing the gene to be included for analysis.
#' @param features Character string. Vector of \code{gene_id} for analysis. Should match \code{gene_id} column of \code{MarvelObject$GeneFeature}.
#' @param point.size Numeric value. Size of data points on reduced dimension space.
#' @param point.alpha Numeric value. Transparency of the data points on reduced dimension space. Take any values between 0 to 1. The smaller the value, the more transparent the data points will be.
#' @param point.stroke Numeric value. The thickness of the outline of the data points. The larger the value, the thicker the outline of the data points.
#' @param pcs Numeric vector. The principal components (PCs) to plot. Default is the first two PCs, i.e., \code{c(1,2)}. If a vector of 3 is specified, a 3D scatterplot is returned.
#' @param mode Character string. Specify \code{"pca"} for linear dimension reduction analysis or \code{"umap"} for non-linear dimension reduction analysis. Specify \code{"elbow.plot"} to return eigen values. Default is \code{"pca"}.
#' @param seed.umap Numeric value. Only applicable when \code{mode} set to \code{"umap"}. To sure reproducibility of analysis. Default value is \code{42}.
#' @param npc.umap Numeric value. Only applicable when \code{mode} set to \code{"umap"}. Incidate the number of PCs to include for UMAP. Default value is \code{30}.
#' @param remove.outliers Logical value. If set to \code{TRUE}, outliers will be removed. Outliers defined as data points beyond 1.5 times the interquartile range (IQR) from the 1st and 99th percentile. Default is \code{FALSE}.
#' @param npc.elbow.plot  Numeric value. Only applicable when \code{mode} set to \code{"elbow.plot"}. Incidate the number of PCs to for elbow plot. Default value is \code{50}.
#'
#' @return An object of class S3 containing with new slots \code{MarvelObject$PCA$Exp$Results}, \code{MarvelObject$PCA$Exp$Plot}, and \code{MarvelObject$PCA$Exp$Plot.Elbow}.
#'
#' @importFrom plyr join
#' @import methods
#' @import ggplot2
#' @importFrom grDevices hcl
#'
#' @export
#'
#' @examples
#' marvel.demo <- readRDS(system.file("extdata/data", "marvel.demo.rds", package="MARVEL"))
#'
#' # Define genes for analysis
#' gene_ids <- marvel.demo$Exp$gene_id
#'
#' # PCA
#' marvel.demo <- RunPCA.Exp(MarvelObject=marvel.demo,
#'                           sample.ids=marvel.demo$SplicePheno$sample.id,
#'                           cell.group.column="cell.type",
#'                           cell.group.order=c("iPSC", "Endoderm"),
#'                           min.cells=5,
#'                           features=gene_ids,
#'                           point.size=2
#'                           )
#'
#' # Check outputs
#' head(marvel.demo$PCA$Exp$Results$ind$coord)
#' marvel.demo$PCA$Exp$Plot

RunPCA.Exp <- function(MarvelObject, sample.ids=NULL, cell.group.column, cell.group.order=NULL, cell.group.colors=NULL,
                       features, min.cells=25,
                       point.size=0.5, point.alpha=0.75, point.stroke=0.1,
                       pcs=c(1,2),
                       mode="pca", seed.umap=42, npc.umap=30,
                       remove.outliers=FALSE, npc.elbow.plot=50
                       ) {

    # Define arguments
    MarvelObject <- MarvelObject
    df <- MarvelObject$Exp
    df.pheno <- MarvelObject$SplicePheno
    df.feature <- MarvelObject$GeneFeature
    sample.ids <- sample.ids
    cell.group.column <- cell.group.column
    cell.group.order <- cell.group.order
    cell.group.colors <- cell.group.colors
    features <- features
    min.cells <- min.cells
    point.size <- point.size
    point.alpha <- point.alpha
    point.stroke <- point.stroke
    mode <- mode
    seed.umap <- seed.umap
    npc.umap <- npc.umap
    
    # Example arguments
    #MarvelObject <- marvel
    #df <- MarvelObject$Exp
    #df.pheno <- MarvelObject$SplicePheno
    #df.feature <- MarvelObject$GeneFeature
    #sample.ids <- NULL
    #cell.group.column <- cell.group.column
    #cell.group.order <- cell.group.order
    #cell.group.colors <- c("red", "blue", "lightgreen", "purple")
    #features <- gene_ids
    #min.cells <- 25
    #point.size <- 1.5
    #point.alpha <- 0.75
    #point.stroke <- 0.1
    #pcs <- c(1,2)
    #mode <- "umap"
    #seed.umap <- 42
    #npc.umap <- 30
    
    ######################################################################
        
    # Create row names for matrix
    row.names(df) <- df$gene_id
    df$gene_id <- NULL
    
    # Rename cell group label/impute columns
    names(df.pheno)[which(names(df.pheno)==cell.group.column)] <- "pca.cell.group.label"
    
    # Subset relevant cells: overall
    if(!is.null(sample.ids[1])) {
        
        df.pheno <- df.pheno[which(df.pheno$sample.id %in% sample.ids), ]
        
    }
    
    # Subset relevant cells
        # Check if cell group order is defined
        if(is.null(cell.group.order[1])) {
            
            cell.group.order <- unique(df.pheno$pca.cell.group.label)
            
        }
    
        # Cell group
        index <- which(df.pheno$pca.cell.group.label %in% cell.group.order)
        df.pheno <- df.pheno[index, ]
        
        # Subset matrix
        df <- df[, df.pheno$sample.id]
        
    # Set factor levels
    levels <- intersect(cell.group.order, unique(df.pheno$pca.cell.group.label))
    df.pheno$pca.cell.group.label <- factor(df.pheno$pca.cell.group.label, levels=levels)

    # Subset features to reduce
    df.feature <- df.feature[which(df.feature$gene_id %in% features), ]
    df <- df[df.feature$gene_id, ]
 
    # Subset events with sufficient cells
    . <- apply(df, 1, function(x) {sum(x != 0)})
    index.keep <- which(. >= min.cells)
    df <- df[index.keep, ]
    df.feature <- df.feature[which(df.feature$gene_id %in% row.names(df)), ]
 
     # Define n PCs to return
     if(nrow(df.pheno) >= 50) {
         
         npc <- 50
         
     } else {
         
         npc <- nrow(df)
         
     }
     
    # Reduce dimension
    res.pca <- FactoMineR::PCA(as.data.frame(t(df)), scale.unit=TRUE, ncp=npc, graph=FALSE)
    
    ##############################################
                               
    # Return elbow plot
    if(mode=="elbow.plot") {
       
       # Retrieve eigenvalues
       results.eigen <- as.data.frame(factoextra::get_eigenvalue(res.pca))
       . <- data.frame("pc"=row.names(results.eigen))
       results.eigen <- cbind.data.frame(., results.eigen)
       results.eigen$pc <- gsub("Dim.", "", results.eigen$pc, fixed=TRUE)
       results.eigen$pc <- as.numeric(results.eigen$pc)
       
       # Subset dimensions
       index <- c(1:npc.elbow.plot)
       results.eigen <- results.eigen[index, ]
        
       # Dot plot
           # Definition
           data <- results.eigen
           x <- data$pc
           y <- data$variance.percent
           maintitle <- ""
           xtitle <- "PC"
           ytitle <- "Variance explained (%)"
           
           # Plot
           plot <- ggplot() +
               geom_point(data, mapping=aes(x=x, y=y), fill="black", size=1) +
               labs(title=maintitle, x=xtitle, y=ytitle) +
               theme(panel.grid.major=element_blank(),
                   panel.grid.minor=element_blank(),
                   panel.background=element_blank(),
                   plot.title = element_text(size=12, hjust=0.5),
                   axis.line=element_line(colour = "black"),
                   axis.title=element_text(size=12),
                   axis.text.x=element_text(size=10, colour="black"),
                   axis.text.y=element_text(size=10, colour="black")
                   )
       
       # Save into new slot
       MarvelObject$PCA$Exp$EigenValues <- results.eigen
       MarvelObject$PCA$Exp$ElbowPlot <- plot
       
       # Return MARVEL object
       return(MarvelObject)
       
    }
                               
    ##############################################
    
    # Define no. of columns for legends
    n.groups <- length(levels(df.pheno$pca.cell.group.label))
    ncol.legends <- ifelse(n.groups < 6, 1, 2)
    
    # Scatterplot: 2D or 3D
    if(mode=="pca") {
        
        if(length(pcs)==2){
            
            # Definition
            data <- as.data.frame(res.pca$ind$coord)
            x <- data[,pcs[1]]
            y <- data[,pcs[2]]
            z <- df.pheno$pca.cell.group.label
            maintitle <- paste(nrow(df), " genes", sep="")
            xtitle <- paste("PC", pcs[1], " (", round(factoextra::get_eigenvalue(res.pca)[pcs[1],2], digits=2), "%)" ,sep="")
            ytitle <- paste("PC", pcs[2], " (", round(factoextra::get_eigenvalue(res.pca)[pcs[2],2], digits=2), "%)" ,sep="")
            legendtitle <- "Cell group"
            
            # Remove outliers
            if(remove.outliers==TRUE){
                
                # Find outliers
                    # PC1
                    lower.limit <- quantile(x, 0.01) - (1.5*IQR(x))
                    upper.limit <- quantile(x, 0.99) + (1.5*IQR(x))
                    index.rm.x <- which(x < lower.limit | x >  upper.limit)
                
                    # PC2
                    lower.limit <- quantile(y, 0.01) - (1.5*IQR(y))
                    upper.limit <- quantile(y, 0.99) + (1.5*IQR(y))
                    index.rm.y <- which(y < lower.limit | y >  upper.limit)
                    
                    # Merge
                    index.rm <- unique(c(index.rm.x, index.rm.y))
                    
                # Remove outliers
                if(length(index.rm) != 0) {
                    
                    data <- data[-index.rm, ]
                    x <- x[-index.rm]
                    y <- y[-index.rm]
                    z <- z[-index.rm]
                
                }
                
                # Track progress
                print(paste(length(index.rm), " outliers removed", sep=""))
                    
            }
            
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
                    legend.text=element_text(size=8),
                    legend.key=element_rect(fill="white")
                    )  +
            guides(fill = guide_legend(override.aes=list(size=2, alpha=point.alpha, stroke=point.stroke), ncol=ncol.legends))

        } else if(length(pcs)==3){
        
            # Definition
            data <- as.data.frame(res.pca$ind$coord)
            x <- data[,pcs[1]]
            y <- data[,pcs[2]]
            z <- data[,pcs[3]]
            group <- df.pheno$pca.cell.group.label
            maintitle <- ""
            xtitle <- ""
            ytitle <- ""
            legendtitle <- "Cell group"
            
            # Remove outliers
            if(remove.outliers==TRUE){
                
                # Find outliers
                    # PC1
                    lower.limit <- quantile(x, 0.01) - (1.5*IQR(x))
                    upper.limit <- quantile(x, 0.99) + (1.5*IQR(x))
                    index.rm.x <- which(x < lower.limit | x >  upper.limit)
                
                    # PC2
                    lower.limit <- quantile(y, 0.01) - (1.5*IQR(y))
                    upper.limit <- quantile(y, 0.99) + (1.5*IQR(y))
                    index.rm.y <- which(y < lower.limit | y >  upper.limit)
                    
                    # PC3
                    lower.limit <- quantile(z, 0.01) - (1.5*IQR(z))
                    upper.limit <- quantile(z, 0.99) + (1.5*IQR(z))
                    index.rm.z <- which(z < lower.limit | z >  upper.limit)
                    
                    # Merge
                    index.rm <- unique(c(index.rm.x, index.rm.y, index.rm.z))
                    
                # Remove outliers
                if(length(index.rm) != 0) {
                    
                    data <- data[-index.rm, ]
                    x <- x[-index.rm]
                    y <- y[-index.rm]
                    z <- z[-index.rm]
                    group <- group[-index.rm]
                
                }
                
                # Track progress
                print(paste(length(index.rm), " outliers removed", sep=""))
                    
            }

            # Color scheme
            if(is.null(cell.group.colors[1])) {
            
                gg_color_hue <- function(n) {
                  hues = seq(15, 375, length = n + 1)
                  hcl(h = hues, l = 65, c = 100)[1:n]
                }
                
                n = length(levels(group))
                cols = gg_color_hue(n)
            
            } else {
                
                cols <- cell.group.colors
                
            }
            
            # Report variance explained
            print(paste("PC", pcs[1], " (", round(factoextra::get_eigenvalue(res.pca)[pcs[1],2], digits=2), "%)" ,sep=""))
            print(paste("PC", pcs[2], " (", round(factoextra::get_eigenvalue(res.pca)[pcs[2],2], digits=2), "%)" ,sep=""))
            print(paste("PC", pcs[3], " (", round(factoextra::get_eigenvalue(res.pca)[pcs[3],2], digits=2), "%)" ,sep=""))

            
            # Plot
            plot <- ggplot(data, aes(x=x, y=y, z=z, fill=group)) +
                theme_void() +
                axes_3D(color="grey75") +
                stat_3D(size=point.size, pch=21, alpha=point.alpha, stroke=point.stroke) +
                scale_fill_manual(values=cols)  +
                #theme(legend.title=element_text(size=8),
                      #legend.text=element_text(size=8),
                      #legend.key=element_rect(fill="white")
                      #) +
                guides(fill = guide_legend(override.aes=list(size=2, alpha=point.alpha, stroke=point.stroke), ncol=ncol.legends))

        }
        
    }
        
    # Retrieve eigenvalues
    results.eigen <- as.data.frame(factoextra::get_eigenvalue(res.pca))
    . <- data.frame("pc"=row.names(results.eigen))
    results.eigen <- cbind.data.frame(., results.eigen)
    results.eigen$pc <- gsub("Dim.", "", results.eigen$pc, fixed=TRUE)
    results.eigen$pc <- as.numeric(results.eigen$pc)
    
    # Non-linear dimension reduction
    if(mode=="umap") {
        
        # Subset first PCs
        data <- as.data.frame(res.pca$ind$coord)
        data.small <- data[,c(1:npc.umap)]
        
        # Reduce dimension
            # Non-linear
            set.seed(seed.umap)
            umap_out <- umap::umap(data.small)
            
        # Scatterplot: Annotate by donor ID
            # Definition
            data <- as.data.frame(umap_out$layout)
            x <- data[,1]
            y <- data[,2]
            z <- df.pheno$pca.cell.group.label
            maintitle <- ""
            xtitle <- "UMAP-1"
            ytitle <- "UMAP-2"
            legendtitle <- "Cell group"
            
            # Color scheme
            if(is.null(cell.group.colors[1])) {
            
                gg_color_hue <- function(n) {
                  hues = seq(15, 375, length = n + 1)
                  hcl(h = hues, l = 65, c = 100)[1:n]
                }
                
                n = length(levels(group))
                cols = gg_color_hue(n)
            
            } else {
                
                cols <- cell.group.colors
                
            }

            # Plot
            plot <- ggplot() +
                geom_point(data, mapping=aes(x=x, y=y, fill=z), size=1.5, pch=21, alpha=0.8, stroke=0.1) +
                scale_fill_manual(values=cols) +
                labs(title=maintitle, x=xtitle, y=ytitle, fill=legendtitle) +
                theme(panel.grid.major=element_blank(),
                    panel.grid.minor=element_blank(),
                    panel.background=element_blank(),
                    plot.title = element_text(size=12, hjust=0.5),
                    axis.line=element_line(colour = "black"),
                    axis.title=element_text(size=12),
                    axis.text=element_text(size=10, colour="black"),
                    legend.title=element_text(size=8),
                    legend.text=element_text(size=8),
                    legend.key=element_rect(fill="white")
                    )  +
            guides(fill = guide_legend(override.aes=list(size=2, alpha=0.8, stroke=0.1),  ncol=ncol.legends))
    
    }
    
    ######################################################################
    
    # Save to new slot
    MarvelObject$PCA$Exp$RawData <- df
    MarvelObject$PCA$Exp$Results <- res.pca
    MarvelObject$PCA$Exp$Plot.Data <- data
    MarvelObject$PCA$Exp$Plot <- plot
    MarvelObject$PCA$Exp$EigenValues <- results.eigen

    return(MarvelObject)
        
}
