#' @title Plot Percent Spliced-in (PSI) Values
#'
#' @description
#' \code{PlotValues.PSI} plots percent spliced-in (PSI) values across different groups of cells.
#'
#' @details
#' This function plots percent spliced-in (PSI) across different groups of cells. Violin plot is used for PSI values.
#'
#' @param MarvelObject S3 object generated from \code{CreateMarvelObject} or \code{ComputePSI} function.
#' @param cell.type.columns List. To indicate which columns in the \code{MarvelObject$SplicePheno} slot to refer to when filtering samples (cells) for analysis. Each item in the list defines the columns for each group to plot.
#' @param cell.type.variables List. To indicate which specific variables of the corresponding columns to keep the samples (cells). Each item in the list should be the same length as the corresponding item in the \code{cell.type.columns} list.
#' @param cell.type.labels Character string. To indicate the cell group labels on the x-axis. Should be same length as the number of items in \code{cell.type.columns} and \code{cell.type.variables} lists.
#' @param feature Character string. \code{tran_id} for plotting. Should match \code{tran_id} column of \code{MarvelObject$ValidatedSpliceFeature} slot.
#' @param maintitle Character string. Column to use as plot main title as per \code{MarvelObject$ValidatedSpliceFeature}. Default is \code{"gene_short_name"} column.
#' @param xlabels.size Numeric value. Size of x-axis labels as per \code{ggplot2} function. Default is 8.
#' @param n.cells.jitter.threshold Numeric value. Maximum number of cells for jitter points. Cells are randomly downsampled to show on jitter plot. Useful when there are large number of cells so that individual jitter points do not overcrowd the violin plot.
#' @param n.cells.jitter.seed Numeric value. Cells downsampled are reproducible.
#' @param n.cells Numeric value. The minimum no. of cells expressing the splicing event to be included for analysis. Please refer to \code{AssignModality} function help page for more details.
#' @param sigma.sq Numeric value. The variance threshold below which the included/excluded modality will be defined as primary sub-modality, and above which it will be defined as dispersed sub-modality. Please refer to \code{AssignModality} function help page for more details. Default is 0.001.
#' @param bimodal.adjust Logical. When set to \code{TRUE}, MARVEL will identify false bimodal modalities and reassign them as included/excluded modality. Please refer to \code{AssignModality} function help page for more details.
#' @param seed Numeric value. Ensure the \code{fitdist} function returns the same values for alpha and beta paramters each time this function is executed using the same random number generator. Please refer to \code{AssignModality} function help page for more details.
#' @param modality.column Character string. Can take the value \code{"modality"}, \code{"modality.var"} or \code{"modality.bimodal.adj"}. Please refer to \code{AssignModality} function help page for more details. Default is \code{"modality.bimodal.adj"}.
#' @param scale.y Logical value. If set to \code{FALSE}, the y-axis of will not be scaled to 0 to 1, but will instead follow the range of the data values. Useful when most PSI values are extremely small (< 0.02) or big (> 0.98). Default is \code{TRUE}.
#' @param n.cells.jitter.threshold Numeric value. Maximum number of cells for jitter points. Cells are randomly downsampled to show on jitter plot. Useful when there are large number of cells so that individual jitter points do not overcrowd the violin plot. Specified together with \code{n.cells.jitter.seed}
#' @param n.cells.jitter.seed Numeric value. Cells downsampled are reproducible. Specified together with \code{n.cells.jitter.threshold}.
#'
#' @export
#'
#' @return An object of class S3 containing all the original slots as inputted by the user in addition to one new slot named \code{MarvelObject$adhocPlot$PSI}.
#'
#' @author Sean Wen <sean.wenwx@gmail.com>
#'
#' @importFrom plyr join
#' @import stats
#' @import methods
#' @import ggplot2
#'
#' @examples
#' # Load input
#' marvel <- readRDS(system.file("extdata/Data", "MarvelObject.rds", package="MARVEL"))
#'
#' # Define event to plot
#' tran_id <- marvel$DE$PSI$Table$tran_id[1]
#'
#' # Run example
#' marvel <- PlotValues.PSI(MarvelObject=marvel,
#'                          cell.type.columns=list(c("cell.type"), c("cell.type")),
#'                          cell.type.variables=list(list("iPSC"), list("Endoderm")),
#'                          cell.type.labels=c("iPSC", "Endoderm"),
#'                          feature=tran_id,
#'                          xlabels.size=8,
#'                          n.cells=2,
#'                          bimodal.adjust=TRUE,
#'                          seed=1,
#'                          )
#'
#' # Check output
#' marvel$adhocPlot$PSI

PlotValues.PSI <- function(MarvelObject, cell.type.columns, cell.type.variables, cell.type.labels, feature, maintitle="gene_short_name", xlabels.size=8, n.cells.jitter.threshold=NULL, n.cells.jitter.seed=NULL, n.cells, sigma.sq=0.001, bimodal.adjust, seed, modality.column="modality.bimodal.adj", scale.y=TRUE) {
    
    # Define arguments
    df <- do.call(rbind.data.frame, MarvelObject$PSI)
    df.pheno <- MarvelObject$SplicePheno
    df.feature <- do.call(rbind.data.frame, MarvelObject$SpliceFeatureValidated)
    cell.type.columns <- cell.type.columns
    cell.type.variables <- cell.type.variables
    cell.type.labels <- cell.type.labels
    feature <- feature
    maintitle <- maintitle
    xlabels.size <- xlabels.size
    n.cells <- n.cells
    sigma.sq <- sigma.sq
    bimodal.adjust <- bimodal.adjust
    seed <- seed
    modality.column <- modality.column
    scale.y <- scale.y
        
    # Example arguments
    #df <- do.call(rbind.data.frame, marvel$PSI)
    #df.pheno <- marvel$SplicePheno
    #df.feature <- do.call(rbind.data.frame, marvel$SpliceFeatureValidated)
    #cell.type.columns <- list(c("sample.type", "cell.type"), c("sample.type", "cell.type"))
    #cell.type.variables <- list(list("Single Cell", "iPSC"), list("Single Cell", "Endoderm"))
    #cell.type.labels <- c("Group 1", "Group 2")
    #feature="chr13:43059394:43059714:+@chr13:43062190:43062295"
    #maintitle="gene_short_name"
    #xlabels.size=8.5
    #n.cells.jitter.threshold=50
    #n.cells.jitter.seed=1
    #n.cells=25
    #sigma.sq=0.001
    #bimodal.adjust=TRUE
    #seed=1
    #modality.column="modality.bimodal.adj"
    #scale.y=TRUE
    
    # Subset relevant feature
    df.feature <- df.feature[which(df.feature$tran_id == feature), , drop=FALSE]
    df <- df[which(df$tran_id == feature), ]
    
    # Create ad hoc marvel object
    s3 <- list()
    class(s3) <- "Marvel"
    s3$PSI <- list(df)
    s3$SplicePheno <- df.pheno
    s3$SpliceFeatureValidated <- list(df.feature)

    # Assign modalities
    mod <- NULL
    
    for(i in 1:length(cell.type.columns)) {
                
        . <- AssignModality(MarvelObject=s3,
                       cell.type.columns=cell.type.columns[[i]],
                       cell.type.variables=cell.type.variables[[i]],
                       n.cells=n.cells,
                       sigma.sq=sigma.sq,
                       bimodal.adjust=TRUE,
                       seed=1
                       )
                       
        mod[i] <- .$Modality$Results[, modality.column]
        
    
    }
    
    modality <- data.frame("cell.type.label"=cell.type.labels, "modality"=mod, stringsAsFactors=FALSE)
    
    #######################################################################################
    
    # Create row names for matrix
    row.names(df) <- df$tran_id
    df$tran_id <- NULL
    
    # Subset relevant feature
    df.feature.small <- df.feature[which(df.feature$tran_id ==  feature), , drop=FALSE]
    df.small <- df[feature, , drop=FALSE]
    df.small <- as.data.frame(t(df.small))
    names(df.small) <- "psi"
    df.small$sample.id <- row.names(df.small)
    row.names(df.small) <- NULL
    
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
    
    # Annotate group labels
    df.small <- join(df.small, md, by="sample.id", type="left")
    
    # Remove un-defined samples
    df.small <- df.small[!is.na(df.small$cell.type.label), ]
    
    # Compute statistics
        # n cells per cell type
        n.cells <- tapply(df.small$psi, df.small$cell.type.label, function(x) {sum(!is.na(x))})
        n.cells <- data.frame("cell.type.label"=names(n.cells), "freq"=n.cells, stringsAsFactors=FALSE)
        row.names(n.cells) <- NULL
        
        n.cells$label <- paste(n.cells$cell.type.label, "\n", "(n=", n.cells$freq, ")", "\n", modality$modality, sep="")
        #n.cells$label <- gsub(".", "\n", n.cells$label, fixed=TRUE)
        
        # Average
        ave <- tapply(df.small$psi, df.small$cell.type, function(x) {mean(x, na.rm=TRUE)})
        ave <- data.frame("cell.type.label"=names(ave), "average"=ave, stringsAsFactors=FALSE)
        row.names(ave) <- NULL
        
    #######################################################################################
    
    if(is.null(n.cells.jitter.threshold)) {
    
        # Violin plot
            # Definition
            data <- df.small
            x <- data$cell.type.label
            y <- data$psi
            z <- data$cell.type.label
            maintitle <- df.feature.small[, maintitle]
            ytitle <- "PSI"
            xtitle <- ""
            xlabels <- n.cells$label
        
            # Plot
            if(scale.y==TRUE) {
                
                plot <- ggplot() +
                    geom_violin(data, mapping=aes(x=x, y=y, fill=z), color="gray", scale="width") +
                    geom_jitter(data, mapping=aes(x=x, y=y), position=position_jitter(width=0.1, height=0), size=0.001) +
                    stat_summary(data, mapping=aes(x=x, y=y), geom="point", fun="mean", fill="red", col="black", size=2, shape=23) +
                    scale_x_discrete(labels=xlabels) +
                    scale_y_continuous(breaks=seq(0, 1, by=0.25), limits=c(0, 1)) +
                    labs(title=maintitle, x=xtitle, y=ytitle) +
                    theme(panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(),
                        panel.background = element_blank(),
                        panel.border=element_blank(),
                        plot.title=element_text(hjust = 0.5, size=12),
                        plot.subtitle=element_text(hjust = 0.5, size=12),
                        axis.line.y.left = element_line(color="black"),
                        axis.line.x = element_line(color="black"),
                        axis.title=element_text(size=12),
                        axis.text.x=element_text(size=xlabels.size, colour="black"),
                        axis.text.y=element_text(size=10, colour="black"),
                        legend.position="none"
                        )
            
            } else {
                
                plot <- ggplot() +
                    geom_violin(data, mapping=aes(x=x, y=y, fill=z), color="gray", scale="width") +
                    geom_jitter(data, mapping=aes(x=x, y=y), position=position_jitter(width=0.1, height=0), size=0.001) +
                    stat_summary(data, mapping=aes(x=x, y=y), geom="point", fun="mean", fill="red", col="black", size=2, shape=23) +
                    scale_x_discrete(labels=xlabels) +
                    labs(title=maintitle, x=xtitle, y=ytitle) +
                    theme(panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(),
                        panel.background = element_blank(),
                        panel.border=element_blank(),
                        plot.title=element_text(hjust = 0.5, size=12),
                        plot.subtitle=element_text(hjust = 0.5, size=12),
                        axis.line.y.left = element_line(color="black"),
                        axis.line.x = element_line(color="black"),
                        axis.title=element_text(size=12),
                        axis.text.x=element_text(size=xlabels.size, colour="black"),
                        axis.text.y=element_text(size=10, colour="black"),
                        legend.position="none"
                        )
                
            }
                    
        # Save into new slow
        MarvelObject$adhocPlot$PSI <- plot
    
    } else if(n.cells.jitter.threshold != 0) {
        
        # Add ids to trace entries
        df.small$id <- c(1:nrow(df.small))
        
        # Downsample data points for jitter plot
        df.small$psi.downsampled <- df.small$psi
        set.seed(n.cells.jitter.seed)
        
        cell.type.labels <- levels(df.small$cell.type.label)
        
        .list <- list()
        
        for(i in 1:length(cell.type.labels)) {
            
            . <- df.small[which(df.small$cell.type.label==cell.type.labels[i]), ]
                        
            if(nrow(.) > n.cells.jitter.threshold) {
                
                set.seed(n.cells.jitter.seed)
                
                index <- sample(1:nrow(.), size=n.cells.jitter.threshold)
                
                .$psi.downsampled[-index] <- NA
            
                .list[[i]] <- .
                
            } else {
                
                .list[[i]] <- .
                
            }
            
        }
        
        df.small <- do.call(rbind.data.frame, .list)
        df.small <- df.small[order(df.small$id), , drop=FALSE]
        df.small$id <- NULL
        
        # Violin plot
            # Definition
            data <- df.small
            x <- data$cell.type.label
            y <- data$psi
            z <- data$cell.type.label
            maintitle <- df.feature.small[, maintitle]
            ytitle <- "PSI"
            xtitle <- ""
            xlabels <- n.cells$label
        
            x.jitter <- data$cell.type
            y.jitter <- data$psi.downsampled

            # Plot
            
            if(scale.y==TRUE) {
                
                plot <- ggplot() +
                    geom_violin(data, mapping=aes(x=x, y=y, fill=z), color="gray", scale="width") +
                    geom_jitter(data, mapping=aes(x=x.jitter, y=y.jitter), position=position_jitter(width=0.1, height=0), size=0.001) +
                    stat_summary(data, mapping=aes(x=x, y=y), geom="point", fun="mean", fill="red", col="black", size=2, shape=23) +
                    scale_x_discrete(labels=xlabels) +
                    scale_y_continuous(breaks=seq(0, 1, by=0.25), limits=c(0, 1)) +
                    labs(title=maintitle, x=xtitle, y=ytitle) +
                    theme(panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(),
                        panel.background = element_blank(),
                        panel.border=element_blank(),
                        plot.title=element_text(hjust = 0.5, size=12),
                        plot.subtitle=element_text(hjust = 0.5, size=12),
                        axis.line.y.left = element_line(color="black"),
                        axis.line.x = element_line(color="black"),
                        axis.title=element_text(size=10),
                        axis.text=element_text(size=10),
                        axis.text.x=element_text(size=xlabels.size, colour="black"),
                        axis.text.y=element_text(size=10, colour="black"),
                        legend.position="none",
                        legend.title=element_text(size=10),
                        legend.text=element_text(size=10)
                        )
                    
            } else {
                
                plot <- ggplot() +
                    geom_violin(data, mapping=aes(x=x, y=y, fill=z), color="gray", scale="width") +
                    geom_jitter(data, mapping=aes(x=x.jitter, y=y.jitter), position=position_jitter(width=0.1, height=0), size=0.001) +
                    stat_summary(data, mapping=aes(x=x, y=y), geom="point", fun="mean", fill="red", col="black", size=2, shape=23) +
                    scale_x_discrete(labels=xlabels) +
                    labs(title=maintitle, x=xtitle, y=ytitle) +
                    theme(panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(),
                        panel.background = element_blank(),
                        panel.border=element_blank(),
                        plot.title=element_text(hjust = 0.5, size=12),
                        plot.subtitle=element_text(hjust = 0.5, size=12),
                        axis.line.y.left = element_line(color="black"),
                        axis.line.x = element_line(color="black"),
                        axis.title=element_text(size=12),
                        axis.text.x=element_text(size=xlabels.size, colour="black"),
                        axis.text.y=element_text(size=10, colour="black"),
                        legend.position="none"
                        )
                
                
                
            }
            
        # Save into new slow
        MarvelObject$adhocPlot$PSI <- plot
                
    } else {
        
        # Violin plot
            # Definition
            data <- df.small
            x <- data$cell.type.label
            y <- data$psi
            z <- data$cell.type.label
            maintitle <- df.feature.small[, maintitle]
            ytitle <- "PSI"
            xtitle <- ""
            xlabels <- n.cells$label
            
            # Plot
            if(scale.y==TRUE) {
                
                plot <- ggplot() +
                    geom_violin(data, mapping=aes(x=x, y=y, fill=z), color="gray", scale="width") +
                    #geom_jitter(data, mapping=aes(x=x, y=y), position=position_jitter(width=0.1, height=0), size=0.001) +
                    stat_summary(data, mapping=aes(x=x, y=y), geom="point", fun="mean", fill="red", col="black", size=2, shape=23) +
                    scale_x_discrete(labels=xlabels) +
                    scale_y_continuous(breaks=seq(0, 1, by=0.25), limits=c(0, 1)) +
                    labs(title=maintitle, x=xtitle, y=ytitle) +
                    theme(panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(),
                        panel.background = element_blank(),
                        panel.border=element_blank(),
                        plot.title=element_text(hjust = 0.5, size=12),
                        plot.subtitle=element_text(hjust = 0.5, size=12),
                        axis.line.y.left = element_line(color="black"),
                        axis.line.x = element_line(color="black"),
                        axis.title=element_text(size=12),
                        axis.text.x=element_text(size=xlabels.size, colour="black"),
                        axis.text.y=element_text(size=10, colour="black"),
                        legend.position="none"
                        )
                        
            } else {
                
                
                
            }

        # Save into new slow
        MarvelObject$adhocPlot$PSI <- plot
        
    }
    
    return(MarvelObject)
    
}
