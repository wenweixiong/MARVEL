#' @title Integrated splicing and gene non-linear dimension reduction analysis
#'
#' @description Non-linear dimension reduction analysis based on both splicing and gene expression data.
#'
#' @param MarvelObject Marvel object. S3 object generated from \code{RunPCA} function.
#' @param n.dim. Indicate the first number of principal components to use for analysis . Default value is \code{20}, i.e., the first 20 PCs.
#' @param seed Numeric value. To sure reproducibility of analysis. Default value is \code{42}.
#' @param sample.ids Character strings. Specific cells to plot.
#' @param cell.group.column Character string. The name of the sample metadata column in which the variables will be used to label the cell groups on the PCA.
#' @param cell.group.order Character string. The order of the variables under the sample metadata column specified in \code{cell.group.column} to appear in the PCA cell group legend.
#' @param cell.group.colors Character string. Vector of colors for the cell groups specified for PCA analysis using \code{cell.type.columns} and \code{cell.group.order}. If not specified, default \code{ggplot2} colors will be used.
#' @param point.size Numeric value. Size of data points on reduced dimension space.
#' @param point.alpha Numeric value. Transparency of the data points on reduced dimension space. Take any values between 0 to 1. The smaller the value, the more transparent the data points will be.
#' @param point.stroke Numeric value. The thickness of the outline of the data points. The larger the value, the thicker the outline of the data points.
#'
#' @return An object of class S3 containing with new slots \code{ MarvelObject$PCA$PSI$ElbowPlot} or \code{MarvelObject$PCA$Exp$ElbowPlot} when \code{level} option set to \code{"splicing"} or \code{"gene"}, respectively.
#'
#' @importFrom plyr join
#' @import methods
#' @import ggplot2
#'
#' @export
#'

RunPCA.PSI.Exp <- function(MarvelObject, n.dim=20, seed=42,
                           sample.ids=NULL,
                           cell.group.column, cell.group.order, cell.group.colors=NULL,
                           point.size=0.5, point.alpha=0.75, point.stroke=0.1
                           ) {

    # Define arguments
    n.dim <- n.dim
    seed <- seed
    sample.ids <- sample.ids
    cell.group.column <- cell.group.column
    cell.group.order <- cell.group.order
    cell.group.colors <- cell.group.colors
    point.size <- point.size
    point.alpha <- point.alpha
    point.stroke <- point.stroke
    
    # Define options
    #MarvelObject <- marvel
    #n.dim <- 10
    #seed <- 42
    #sample.ids <- NULL
    #cell.group.column <- "cell.type"
    #cell.group.order <- c("0-hrs", "24-hrs", "48-hrs", "72-hrs")
    #cell.group.colors <- NULL
    #point.size <- 2
    #point.alpha <- 0.8
    #point.stroke <- 0.1
    
    ##########################################

    # Create feature matrix
        # Retrieve splicing PCs
        df.psi <- as.data.frame(MarvelObject$PCA$PSI$Results$ind$coord)
        df.psi <- df.psi[,c(1:n.dim)]

        # Retrieve gene PCA raw data
        df.exp <- as.data.frame(MarvelObject$PCA$Exp$Results$ind$coord)
        df.exp <- df.exp[,c(1:n.dim)]

        # Subset overlapping samples
        overlap <- intersect(row.names(df.psi), row.names(df.exp))
        df.psi <- df.psi[overlap, ]
        df.exp <- df.exp[overlap, ]

        # Merge
        names(df.psi) <- paste(names(df.psi), "_splicing", sep="")
        names(df.exp) <- paste(names(df.exp), "_gene", sep="")
        df <- cbind.data.frame(df.psi, df.exp)

    # Match cell metadata
        # Retreve sample metadata
        df.pheno <- marvel$SplicePheno

        # Subset overlapping samples
        index <- which(df.pheno$sample.id %in% row.names(df))
        df.pheno <- df.pheno[index, ]

        # Match sample order
        df.pheno$sample.id <- factor(df.pheno$sample.id, levels=row.names(df))
        df.pheno <- df.pheno[order(df.pheno$sample.id), ]
        df.pheno$sample.id <- as.character(df.pheno$sample.id)

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
        df <- df[df.pheno$sample.id, ]
        
    # Set factor levels
    levels <- intersect(cell.group.order, unique(df.pheno$pca.cell.group.label))
    df.pheno$pca.cell.group.label <- factor(df.pheno$pca.cell.group.label, levels=levels)

    # Reduce dimension
        # Non-linear
        set.seed(seed)
        umap_out <- umap::umap(df)

    # Scatterplot
        # Definition
        data <- as.data.frame(umap_out$layout)
        x <- data[,1]
        y <- data[,2]
        z <- df.pheno$pca.cell.group.label
        maintitle <- ""
        xtitle <- "UMAP-1"
        ytitle <- "UMAP-2"
        legendtitle <- "Group"

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
            geom_point(data, mapping=aes(x=x, y=y, fill=z), size=point.size, pch=21, alpha=point.alpha, stroke=point.stroke) +
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
                legend.text=element_text(size=8)
                )  +
        guides(fill = guide_legend(override.aes=list(size=2, alpha=0.8, stroke=0.1), ncol=1))

    ######################################################################
    
    # Save to new slot
    MarvelObject$PCA$Integrated$Plot <- plot
    MarvelObject$PCA$Integrated$Results <- data

    return(MarvelObject)
        
}
