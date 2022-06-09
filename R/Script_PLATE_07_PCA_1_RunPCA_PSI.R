#' @title Principle component analysis for splicing data
#'
#' @description Performs principle component analysis using PSI values.
#'
#' @param MarvelObject Marvel object. S3 object generated from \code{TransformExpValues} function.
#' @param cell.group.list List of character strings. Each element of the list is a vector of cell IDs corresponding to a cell group. The name of the element will be the cell group label.
#' @param min.cells Numeric value. The minimum no. of cells expressing the splicing event to be included for analysis.
#' @param features Character string. Vector of \code{tran_id} for analysis. Should match \code{tran_id} column of \code{MarvelObject$ValidatedSpliceFeature}.
#' @param point.size Numeric value. Size of data points on reduced dimension space.
#' @param point.alpha Numeric value. Transparency of the data points on reduced dimension space. Take any values between 0 to 1. The smaller the value, the more transparent the data points will be.
#' @param point.stroke Numeric value. The thickness of the outline of the data points. The larger the value, the thicker the outline of the data points.
#' @param event.type Character string. Indicate which splicing event type to include for analysis. Can take value \code{"SE"}, \code{"MXE"}, \code{"RI"}, \code{"A5SS"}, or \code{"A3SS"} which represents skipped-exon (SE), mutually-exclusive exons (MXE), retained-intron (RI), alternative 5' splice site (A5SS), and alternative 3' splice site (A3SS), respectively.
#' @param seed Numeric value. Ensures imputed values for NA PSIs are reproducible.
#' @param method.impute Character string. Indicate the method for imputing missing PSI values (low coverage). \code{"random"} method randomly assigns any values between 0-1. \code{"population.mean"} method uses the mean PSI value for each cell population. Default option is \code{"population.mean"}.
#' @param retrieve.non.outliers Logical. If set to \code{TRUE}, this function will retrieve \code{sample.id} of non-outliers based on the intial PCA. Define the non-outliers based on the initial PCA coordinates. Use in conjunction with arguments \code{pc1.min}, \code{pc1.max}, \code{pc2.min}, and \code{pc2.max}.
#' @param pc1.min Numeric value. When \code{retrieve.non.outliers} set to \code{TRUE}. To indicate the PC1 value above which to retrieve the sample IDs. Use in conjunction with \code{pc1.max}, \code{pc2.min}, and \code{pc2.max}.
#' @param pc1.max Numeric value. When \code{retrieve.non.outliers} set to \code{TRUE}. To indicate the PC1 value below which to retrieve the sample IDs. Use in conjunction with \code{pc1.min}, \code{pc2.min}, and \code{pc2.max}.
#' @param pc2.min Numeric value. When \code{retrieve.non.outliers} set to \code{TRUE}. To indicate the PC2 value above which to retrieve the sample IDs. Use in conjunction with \code{pc1.max}, \code{pc1.max}, and \code{pc2.max}.
#' @param pc2.max Numeric value. When \code{retrieve.non.outliers} set to \code{TRUE}. To indicate the PC2 value below which to retrieve the sample IDs. Use in conjunction with \code{pc1.min}, \code{pc2.min}, and \code{pc2.min}.
#' @param remove.outliers Logical. If set to \code{TRUE}, re-run PCA by only including non-outliers. Use after running the function with \code{retrieve.non.outliers} set to \code{TRUE}.
#' @param cell.group.colors Character string. Vector of colors for the cell groups specified for PCA analysis using \code{cell.type.columns}, \code{cell.type.variable}, and \code{cell.type.labels}. If not specified, default \code{ggplot2} colors will be used.
#'
#' @return An object of class S3 containing with new slots \code{MarvelObject$PCA$PSI$Results}, \code{MarvelObject$PCA$PSI$Plot}, and \code{MarvelObject$PCA$PSI$Plot.Elbow}.
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


RunPCA.PSI <- function(MarvelObject, cell.group.list, min.cells=25, features,
                       point.size=0.5, point.alpha=0.75, point.stroke=0.1,
                       event.type, seed=1, method.impute="population.mean",
                       retrieve.non.outliers=FALSE, pc1.min=NULL, pc1.max=NULL, pc2.min=NULL, pc2.max=NULL,
                       remove.outliers=FALSE, cell.group.colors=NULL
                       ) {

    # Define arguments
    df <- do.call(rbind.data.frame, MarvelObject$PSI)
    df.pheno <- MarvelObject$SplicePheno
    df.feature <- do.call(rbind.data.frame, MarvelObject$SpliceFeatureValidated)
    cell.group.list <- cell.group.list
    min.cells <- min.cells
    features <- features
    point.size <- point.size
    point.alpha <- point.alpha
    point.stroke <- point.stroke
    seed <- seed
    method.impute <- method.impute
    event.type <- event.type
    retrieve.non.outliers <- retrieve.non.outliers
    pc1.min <- pc1.min
    pc1.max <- pc1.max
    pc2.min <- pc2.min
    pc2.max <- pc2.max
    remove.outliers <- remove.outliers
    cell.group.colors <- cell.group.colors
    
    # Example arguments
    #df <- do.call(rbind.data.frame, marvel$PSI)
    #df.pheno <- marvel$SplicePheno
    #df.feature <- do.call(rbind.data.frame, marvel$SpliceFeatureValidated)
    #cell.group.list <- cell.group.list
    #min.cells <- 25
    #features <- tran_ids
    #method.impute <- "random"
    #seed <- 1
    #event.type <- c("SE", "MXE", "RI", "A5SS", "A3SS", "ALE", "AFE")
    #point.size <- 2.5
    #point.alpha <- 0.75
    #point.stroke <- 0.1
    #retrieve.non.outliers <- FALSE
    #pc1.min <- NULL
    #pc1.max <- NULL
    #pc2.min <- NULL
    #pc2.max <- NULL
    #cell.group.colors <- NULL
    #remove.outliers <- FALSE
    
    ######################################################################
    
    # Create row names for matrix
    row.names(df) <- df$tran_id
    df$tran_id <- NULL
    
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

    # Subset relevant event type
    df.feature <- df.feature[which(df.feature$event_type %in% event.type), ]
    df <- df[df.feature$tran_id, ]
    
    # Subset events with sufficient cells
    . <- apply(df, 1, function(x) {sum(!is.na(x))})
    index.keep <- which(. >= min.cells)
    df <- df[index.keep, ]
    df.feature <- df.feature[which(df.feature$tran_id %in% row.names(df)), ]
    
    # Subset features to reduce
    df.feature <- df.feature[which(df.feature$tran_id %in% features), ]
    df <- df[df.feature$tran_id, ]
    
    # Remove outliers (previously identified)
    if(remove.outliers==TRUE) {
        
        df.pheno <- df.pheno[which(df.pheno$sample.id %in% MarvelObject$PCA$PSI$sample.ids.non.outliers), ]
        df <- df[, df.pheno$sample.id]
        
    }
    
    # Impute missing values
    if(method.impute=="random") {
        
        set.seed(seed)
        df[is.na(df)] <- runif(n=sum(is.na(df)), min=0, max=1)
        
    } else if(method.impute=="population.mean"){
        
        # Define cell groups
        cell.type.labels <- unique(md$cell.type.label)
        
        .list <- list()
        
        for(i in 1:length(cell.type.labels)) {
            
            # Define cell group
            cell.type.label <- cell.type.labels[i]
            
            # Subset cell group
            sample.ids <- md[which(md$cell.type.label==cell.type.label), "sample.id"]
            df.small <- df[, sample.ids]
            
            # Impute
            df.small <- apply(df.small, 1, function(x) {
                
                            #x <- as.numeric(df.small[4, ])
                            x[is.na(x)] <- mean(x, na.rm=TRUE)
                            return(x)
                            
                        })
            
            # Check alignment
            df.small <- as.data.frame(t(df.small))
            #print(table(names(df.small)==sample.ids))
            
            # Save into list
            .list[[i]] <- df.small
            
        }
        
        df.imputed <- do.call(cbind.data.frame, .list)
        
        # Impute values for population-specific events
        set.seed(seed)
        df.imputed[is.na(df.imputed)] <- runif(n=sum(is.na(df.imputed)), min=0, max=1)
        
        # Remove non-expressed population-specific events
        #. <- apply(df.small, 1, function(x) { sum(is.na(x)) })
        #index <- which(. == 0)
        #df.imputed <- df.imputed[index, ]
        #df.feature <- df.feature[which(df.feature$tran_id %in% row.names(df.imputed)), ]
        
        # Check alignment
        table(names(df.imputed)==names(df))
        
        # Fix alignment
        df.imputed <- df.imputed[, names(df)]
        table(names(df.imputed)==names(df))
        
        # Save as new object
        df <- df.imputed
    
    }
    
    ##############################################
    
    # Reduce dimension
    res.pca <- PCA(as.data.frame(t(df)), scale.unit=TRUE, ncp=20, graph=FALSE)
    
    if(retrieve.non.outliers==TRUE) {
        
        # Retrieve coordinates
        . <- as.data.frame(res.pca$ind$coord)
        
        # Subset PC1
        . <- .[which(.$Dim.1 > pc1.min & .$Dim.1 < pc1.max), ]
        
        # Subset PC2
        . <- .[which(.$Dim.2 > pc2.min & .$Dim.2 < pc2.max), ]
        
        # Save non-outliers
        MarvelObject$PCA$PSI$sample.ids.non.outliers <- row.names(.)
        
    }
    
    # Scatterplot
        # Definition
        data <- as.data.frame(res.pca$ind$coord)
        x <- data[,1]
        y <- data[,2]
        z <- df.pheno$cell.type.label
        maintitle <- paste(nrow(df), " splicing events", sep="")
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
                axis.text=element_text(size=10, colour="black"),
                legend.title=element_text(size=8),
                legend.text=element_text(size=8)
                )
    
    # Elbow plot
        # Definition
        if(length(features) < 20) {
            
            data.2 <- as.data.frame(get_eigenvalue(res.pca)[c(1:length(features)),])
            
        } else {
            
            data.2 <- as.data.frame(get_eigenvalue(res.pca)[c(1:20),])
            
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
                axis.text.x=element_text(size=10, colour="black"),
                axis.text.y=element_text(size=10, colour="black"),
                legend.title=element_text(size=8),
                legend.text=element_text(size=8)
                )
                          
    # Save to new slot
    MarvelObject$PCA$PSI$Results <- res.pca
    MarvelObject$PCA$PSI$Plot <- plot
    MarvelObject$PCA$PSI$Plot.Elbow <- plot.2
    
    return(MarvelObject)
        
}
