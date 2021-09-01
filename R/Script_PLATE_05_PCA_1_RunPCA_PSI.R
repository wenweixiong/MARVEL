#' @title Principle Component Analysis for Splicing Data
#'
#' @description
#' \code{RunPCA} performs principle component analysis on splicing data.
#'
#' @details
#' This function performs principle component analysis on splicing data and visualise cells on a reducted dimension space, i.e. 2D.
#'
#' @param MarvelObject S3 object generated from \code{CreateMarvelObject} or \code{ComputePSI} function.
#' @param cell.type.columns Character string. To indicate which columns in the \code{MarvelObject$SplicePheno} slot to refer to when filtering samples (cells) for analysis.
#' @param cell.type.variables List of character string. To indicate which specific variables of the corresponding columns to keep the samples (cells). This should be same length as the \code{cell.type.columns} argument.
#' @param cell.type.labels Character string. To indicate the cell group labels. Should be same length as the number of items in \code{cell.type.columns} and \code{cell.type.variables} lists.
#' @param n.cells Numeric value. The minimum no. of cells expressing the splicing event to be included for analysis.
#' @param features Character string. Vector of \code{tran_id} for analysis. Should match \code{tran_id} column of \code{MarvelObject$ValidatedSpliceFeature}.
#' @param point.size Numeric value. Size of data points on reduced dimension space.
#' @param event.type Character string. Indicate which splicing event type to include for analysis. Can take value \code{"SE"}, \code{"MXE"}, \code{"RI"}, \code{"A5SS"}, or \code{"A3SS"} which represents skipped-exon (SE), mutually-exclusive exons (MXE), retained-intron (RI), alternative 5' splice site (A5SS), and alternative 3' splice site (A3SS), respectively.
#' @param seed Numeric value. Ensures imputed values for NA PSIs are reproducible.
#' @param retrieve.non.outliers Logical. If set to \code{TRUE}, this function will retrieve \code{sample.id} of non-outliers based on the intial PCA. Define the non-outliers based on the initial PCA coordinates. Use in conjunction with arguments \code{pc1.min}, \code{pc1.max}, \code{pc2.min}, and \code{pc2.max}.
#' @param pc1.min Numeric value. When \code{retrieve.non.outliers} set to \code{TRUE}. To indicate the PC1 value above which to retrieve the sample IDs. Use in conjunction with \code{pc1.max}, \code{pc2.min}, and \code{pc2.max}.
#' @param pc1.max Numeric value. When \code{retrieve.non.outliers} set to \code{TRUE}. To indicate the PC1 value below which to retrieve the sample IDs. Use in conjunction with \code{pc1.min}, \code{pc2.min}, and \code{pc2.max}.
#' @param pc2.min Numeric value. When \code{retrieve.non.outliers} set to \code{TRUE}. To indicate the PC2 value above which to retrieve the sample IDs. Use in conjunction with \code{pc1.max}, \code{pc1.max}, and \code{pc2.max}.
#' @param pc2.max Numeric value. When \code{retrieve.non.outliers} set to \code{TRUE}. To indicate the PC2 value below which to retrieve the sample IDs. Use in conjunction with \code{pc1.min}, \code{pc2.min}, and \code{pc2.min}.
#' @param remove.outliers Logical. If set to \code{TRUE}, re-run PCA by only including non-outliers. Use after running the function with \code{retrieve.non.outliers} set to \code{TRUE}.

#' @export
#'
#' @return An object of class S3 containing all the original slots as inputted by the user in addition to one new slot named \code{MarvelObject$PCA$PSI}. Contains 2D scatterplot in \code{MarvelObject$PCA$PSI$Plot} and the corresponding x- and y-coordinates for each sample in \code{MarvelObject$PCA$PSI$Results}. Additional slow named \code{MarvelObject$PCA$PSI$sample.ids.non.outliers} containing non-outlier sample ID will be returned when \code{retrieve.non.outliers} set to \code{TRUE}.
#'
#' @author Sean Wen <sean.wenwx@gmail.com>
#'
#' @importFrom plyr join
#' @import stats
#' @import methods
#' @import FactoMineR
#' @import factoextra
#' @import ggplot2
#'
#' @examples
#' # Load input
#' marvel <- readRDS(system.file("extdata/Data", "MarvelObject.rds", package="MARVEL"))
#'
#' # Define events to reduce dimension
#' . <- do.call(rbind.data.frame, marvel$SpliceFeatureValidated)
#' tran_ids <- .$tran_id
#'
#' # Run example
#' marvel <- RunPCA.PSI(MarvelObject=marvel,
#'                  cell.type.columns=list(c("cell.type"), c("cell.type")),
#'                  cell.type.variables=list(list("iPSC"), list("Endoderm")),
#'                  cell.type.labels=c("iPSC", "Endoderm"),
#'                  n.cells=2,
#'                  features=tran_ids,
#'                  point.size=0.5,
#'                  event.type=c("SE", "MXE", "RI", "A5SS", "A3SS"),
#'                  seed=1
#'                  )
#'
#' # Check output
#' marvel$PCA$PSI$Results
#' marvel$PCA$PSI$Plot
#' marvel$PCA$PSI$Plot.Elbow


RunPCA.PSI <- function(MarvelObject, cell.type.columns, cell.type.variables, cell.type.labels, n.cells, features,
                       point.size=0.5, event.type, seed,
                       retrieve.non.outliers=FALSE, pc1.min=NULL, pc1.max=NULL, pc2.min=NULL, pc2.max=NULL,
                       remove.outliers=FALSE
                       ) {

    # Define arguments
    df <- do.call(rbind.data.frame, MarvelObject$PSI)
    df.pheno <- MarvelObject$SplicePheno
    df.feature <- do.call(rbind.data.frame, MarvelObject$SpliceFeatureValidated)
    cell.type.columns <- cell.type.columns
    cell.type.variables <- cell.type.variables
    cell.type.labels <- cell.type.labels
    n.cells <- n.cells
    features <- features
    point.size <- point.size
    seed <- seed
    event.type <- event.type
    retrieve.non.outliers <- retrieve.non.outliers
    pc1.min <- pc1.min
    pc1.max <- pc1.max
    pc2.min <- pc2.min
    pc2.max <- pc2.max
    remove.outliers <- remove.outliers
    
    # Example arguments
    #df <- do.call(rbind.data.frame, marvel$PSI)
    #df.pheno <- marvel$SplicePheno
    #df.feature <- do.call(rbind.data.frame, marvel$SpliceFeatureValidated)
    #cell.type.columns=list(c("cell.type"), c("cell.type"))
    #cell.type.variables=list(list("iPSC"), list("Endoderm"))
    #cell.type.labels=c("iPSC", "Unknown")
    #n.cells <- 2
    #features <- tran_ids
    #seed <- 1
    #event.type <- c("SE", "MXE", "RI", "A5SS", "A3SS")
    #point.size <- 1
    #retrieve.non.outliers <- TRUE
    #pc1.min <- -20
    #pc1.max <- 9
    #pc2.min <- -5
    #pc2.max <- 5
    
    # Create row names for matrix
    row.names(df) <- df$tran_id
    df$tran_id <- NULL
    
    # Retrieve sample ids for each group
    .list.2 <- list()
    
    for(j in 1:length(cell.type.columns)) {
        
        cell.type.column <- cell.type.columns[[j]]
        cell.type.variable <- cell.type.variables[[j]]
        
        .list <- list()
        
        for(i in 1:length(cell.type.column)) {
        
            .list[[i]] <- df.pheno[which(df.pheno[[cell.type.column[i]]] %in% cell.type.variable[[i]]), "sample.id"]
        
        }
        
        sample.ids <- Reduce(intersect, .list)
        
        .list.2[[j]] <- data.frame("sample.id"=sample.ids, cell.type.label=cell.type.labels[j], stringsAsFactors=FALSE)
        
    }
    
    md <- do.call(rbind.data.frame, .list.2)
       
    # Set factor levels
    md$cell.type.label <- factor(md$cell.type.label, levels=cell.type.labels)
    
    # Subset relevant cells
    df.pheno <- df.pheno[which(df.pheno$sample.id %in% md$sample.id), ]
    df.pheno <- join(df.pheno, md, by="sample.id", type="left")
    df <- df[, which(names(df) %in% df.pheno$sample.id)]

    # Subset relevant event type
    df.feature <- df.feature[which(df.feature$event_type %in% event.type), ]
    df <- df[df.feature$tran_id, ]
    
    # Subset events with sufficient cells
    . <- apply(df, 1, function(x) {sum(!is.na(x))})
    index.keep <- which(. >= n.cells)
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
    set.seed(seed)
    df[is.na(df)] <- runif(n=sum(is.na(df)), min=0, max=1)

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
        
        # Plot
        plot <- ggplot() +
            geom_point(data, mapping=aes(x=x, y=y, color=z), size=point.size, alpha=0.75) +
            labs(title=maintitle, x=xtitle, y=ytitle, color=legendtitle) +
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
