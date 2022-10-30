#' @title Principle component analysis for splicing data
#'
#' @description Performs principle component analysis using PSI values.
#'
#' @param MarvelObject Marvel object. S3 object generated from \code{TransformExpValues} function.
#' @param cell.group.column Character string. The name of the sample metadata column in which the variables will be used to label the cell groups on the PCA.
#' @param cell.group.order Character string. The order of the variables under the sample metadata column specified in \code{cell.group.column} to appear in the PCA cell group legend.
#' @param cell.group.colors Character string. Vector of colors for the cell groups specified for PCA analysis using \code{cell.type.columns} and \code{cell.group.order}. If not specified, default \code{ggplot2} colors will be used.
#' @param sample.ids Character strings. Specific cells to plot.
#' @param min.cells Numeric value. The minimum no. of cells expressing the splicing event to be included for analysis.
#' @param features Character string. Vector of \code{tran_id} for analysis. Should match \code{tran_id} column of \code{MarvelObject$ValidatedSpliceFeature}.
#' @param point.size Numeric value. Size of data points on reduced dimension space.
#' @param point.alpha Numeric value. Transparency of the data points on reduced dimension space. Take any values between 0 to 1. The smaller the value, the more transparent the data points will be.
#' @param point.stroke Numeric value. The thickness of the outline of the data points. The larger the value, the thicker the outline of the data points.
#' @param seed Numeric value. Ensures imputed values for NA PSIs are reproducible.
#' @param method.impute Character string. Indicate the method for imputing missing PSI values (low coverage). \code{"random"} method randomly assigns any values between 0-1. \code{"population.mean"} method uses the mean PSI value for each cell population. Default option is \code{"population.mean"}.
#' @param cell.group.column.impute Character string. Only applicable when \code{method.impute} set to \code{"population.mean"}. The name of the sample metadata column in which the variables will be used to impute missing values.
#'
#' @return An object of class S3 containing with new slots \code{MarvelObject$PCA$PSI$Results} and  \code{MarvelObject$PCA$PSI$Plot}
#'
#' @importFrom plyr join
#' @importFrom stats rnorm runif sd
#' @import methods
#' @import ggplot2
#' @importFrom grDevices hcl
#'
#' @export
#'
#' @examples
#' marvel.demo <- readRDS(system.file("extdata/data", "marvel.demo.rds", package="MARVEL"))
#'
#' # Define splicing events for analysis
#' df <- do.call(rbind.data.frame, marvel.demo$PSI)
#' tran_ids <- df$tran_id
#'
#' # PCA
#' marvel.demo <- RunPCA.PSI(MarvelObject=marvel.demo,
#'                           sample.ids=marvel.demo$SplicePheno$sample.id,
#'                           cell.group.column="cell.type",
#'                           cell.group.order=c("iPSC", "Endoderm"),
#'                           cell.group.colors=NULL,
#'                           min.cells=5,
#'                           features=tran_ids,
#'                           point.size=2
#'                           )
#'
#' # Check outputs
#' head(marvel.demo$PCA$PSI$Results$ind$coord)
#' marvel.demo$PCA$PSI$Plot

RunPCA.PSI <- function(MarvelObject, sample.ids=NULL, cell.group.column, cell.group.order, cell.group.colors=NULL,
                       features, min.cells=25,
                       point.size=0.5, point.alpha=0.75, point.stroke=0.1,
                       seed=1, method.impute="random", cell.group.column.impute=NULL
                       ) {

    # Define arguments
    df <- do.call(rbind.data.frame, MarvelObject$PSI)
    df.pheno <- MarvelObject$SplicePheno
    df.feature <- do.call(rbind.data.frame, MarvelObject$SpliceFeatureValidated)
    sample.ids <- sample.ids
    cell.group.column <- cell.group.column
    cell.group.order <- cell.group.order
    cell.group.colors <- cell.group.colors
    features <- features
    min.cells <- min.cells
    method.impute <- method.impute
    cell.group.column.impute <- cell.group.column.impute
    seed <- seed
    point.size <- point.size
    point.alpha <- point.alpha
    point.stroke <- point.stroke
    
    # Example arguments
    #MarvelObject <- marvel
    #df <- do.call(rbind.data.frame, MarvelObject$PSI)
    #df.pheno <- MarvelObject$SplicePheno
    #df.feature <- do.call(rbind.data.frame, MarvelObject$SpliceFeatureValidated)
    #sample.ids <- sample.ids
    #cell.group.column <- "genotype"
    #cell.group.order <- cell.group.order
    #cell.group.colors <- cell.group.colors
    #features <- tran_ids
    #min.cells <- 20
    #method.impute <- "random"
    #cell.group.column.impute <- "genotype.impute"
    #seed <- 1
    #point.size <- 2.5
    #point.alpha <- 0.75
    #point.stroke <- 0.1
    
    ######################################################################
    
    # Create row names for matrix
    row.names(df) <- df$tran_id
    df$tran_id <- NULL
    
    # Rename cell group label/impute columns
    #names(df.pheno)[which(names(df.pheno)==cell.group.column)] <- "pca.cell.group.label"
    df.pheno$pca.cell.group.label <- df.pheno[[cell.group.column]]
    
    if(!is.null(cell.group.column.impute)) {
      
        #names(df.pheno)[which(names(df.pheno)==cell.group.column.impute)] <- "pca.cell.group.impute"
        df.pheno$pca.cell.group.impute <- df.pheno[[cell.group.column.impute]]
        
    } else {
        
        #df.pheno$pca.cell.group.impute <- df.pheno$pca.cell.group.label
        df.pheno$pca.cell.group.impute <- df.pheno$pca.cell.group.label
        
    }
    
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
    df.feature <- df.feature[which(df.feature$tran_id %in% features), ]
    df <- df[df.feature$tran_id, ]

    # Subset relevant event type
    #df.feature <- df.feature[which(df.feature$event_type %in% event.type), ]
    #df <- df[df.feature$tran_id, ]
    
    # Subset events with sufficient cells
    . <- apply(df, 1, function(x) {sum(!is.na(x))})
    index.keep <- which(. >= min.cells)
    df <- df[index.keep, ]
    df.feature <- df.feature[which(df.feature$tran_id %in% row.names(df)), ]
        
    # Impute missing values
    if(method.impute=="random") {
        
        set.seed(seed)
        df[is.na(df)] <- runif(n=sum(is.na(df)), min=0, max=1)
        
    } else if(method.impute=="population.mean"){
        
        # Define cell groups
        groups <- unique(df.pheno$pca.cell.group.impute)
        
        # Subset events with sufficient cells
        .list <- list()
        
        for(i in 1:length(groups)) {
            
            # Define cell groups
            group <- groups[i]
            
            # Subset cell group
            sample.ids <- df.pheno[which(df.pheno$pca.cell.group.impute==group), "sample.id"]
            df.small <- df[, sample.ids]
            
            # Retrieve expressed events
            . <- apply(df.small, 1, function(x) {sum(!is.na(x))})
            tran_ids <- names(.)[which(. >= min.cells)]
            #message(length(tran_ids))
            
            # Save into list
            .list[[i]] <- tran_ids
            
        }
        
        tran_ids <- Reduce(intersect, .list)

        df.feature <- df.feature[which(df.feature$tran_id %in% tran_ids), ]
        df <- df[df.feature$tran_id ,]
        
        # Impute values for each cell group
        for(i in 1:length(groups)) {
            
            # Define cell groups
            group <- groups[i]
            
            # Subset cell group
            sample.ids <- df.pheno[which(df.pheno$pca.cell.group.impute==group), "sample.id"]
            df.small <- df[, sample.ids]
            
            # Impute
            set.seed(seed)
            
            df.small <- apply(df.small, 1, function(x) {
                            
                            # Example
                            #x <- as.numeric(df.small[4, ])
                            
                            # Find mean, std
                            ave <- mean(x[!is.na(x)], na.rm=TRUE)
                            std.dev <- sd(x[!is.na(x)], na.rm=TRUE)
                            
                            # Impute
                            values <- rnorm(sum(is.na(x)), mean=ave, sd=std.dev)
                            
                            # Jitter values
                            #values <- jitter(values)
                            
                            # Re-code missing values
                            x[is.na(x)] <- values
                            
                            # Return values
                            return(x)
                            
                        })
            
            # Check alignment
            df.small <- as.data.frame(t(df.small))
            #message(table(names(df.small)==sample.ids))
            
            # Save into list
            .list[[i]] <- df.small
            
        }
        
        df <- do.call(cbind.data.frame, .list)
        
        # Match matrix column to phenoData
        df <- df[, df.pheno$sample.id]
    
    }

    ##############################################
    
    # Reduce dimension
    res.pca <- FactoMineR::PCA(as.data.frame(t(df)), scale.unit=TRUE, ncp=20, graph=FALSE)
    
    # Scatterplot
        # Definition
        data <- as.data.frame(res.pca$ind$coord)
        x <- data[,1]
        y <- data[,2]
        z <- df.pheno$pca.cell.group.label
        maintitle <- paste(nrow(df), " splicing events", sep="")
        xtitle <- paste("PC1 (", round(factoextra::get_eigenvalue(res.pca)[1,2], digits=1), "%)" ,sep="")
        ytitle <- paste("PC2 (", round(factoextra::get_eigenvalue(res.pca)[2,2], digits=1), "%)" ,sep="")
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
                axis.text=element_text(size=10, colour="black"),
                legend.title=element_text(size=8),
                legend.text=element_text(size=8)
                )  +
        guides(fill = guide_legend(override.aes=list(size=2, alpha=point.alpha, stroke=point.stroke), ncol=1))
     
    ##############################################
     
    # Save to new slot
    MarvelObject$PCA$PSI$Results <- res.pca
    MarvelObject$PCA$PSI$Plot <- plot
    
    return(MarvelObject)
        
}
