#' @title Plot percent spliced-in (PSI) values
#'
#' @description Violin plot of percent spliced-in (PSI) values across different groups of cells.
#'
#' @param MarvelObject Marvel object. S3 object generated from \code{TransformExpValues} function.
#' @param cell.group.list List of character strings. Each element of the list is a vector of cell IDs corresponding to a cell group. The name of the element will be the cell group label.
#' @param feature Character string. Coordinates of splicing event to plot.
#' @param maintitle Character string. Column to use as plot main title as per \code{MarvelObject$ValidatedSpliceFeature}. Default is \code{"gene_short_name"} column.
#' @param xlabels.size Numeric value. Size of x-axis labels as per \code{ggplot2} function. Default is 8.
#' @param max.cells.jitter Numeric value. Maximum number of cells for jitter points. Cells are randomly downsampled to show on jitter plot. Useful when there are large number of cells so that individual jitter points do not overcrowd the violin plot.
#' @param max.cells.jitter.seed Numeric value. Cells downsampled are reproducible.
#' @param min.cells Numeric value. The minimum no. of cells expressing the splicing event to be included for analysis. Please refer to \code{AssignModality} function help page for more details.
#' @param sigma.sq Numeric value. The variance threshold below which the included/excluded modality will be defined as primary sub-modality, and above which it will be defined as dispersed sub-modality. Please refer to \code{AssignModality} function help page for more details. Default is 0.001.
#' @param bimodal.adjust Logical. When set to \code{TRUE}, MARVEL will identify false bimodal modalities and reassign them as included/excluded modality. Please refer to \code{AssignModality} function help page for more details.
#' @param seed Numeric value. Ensure the \code{fitdist} function returns the same values for alpha and beta paramters each time this function is executed using the same random number generator. Please refer to \code{AssignModality} function help page for more details.
#' @param modality.column Character string. Can take the value \code{"modality"}, \code{"modality.var"} or \code{"modality.bimodal.adj"}. Please refer to \code{AssignModality} function help page for more details. Default is \code{"modality.bimodal.adj"}.
#' @param scale.y.log Logical value. Only applicable when \code{level} set to \code{"splicing"}. If set to \code{TRUE}, the y-axis of will log10-scaled. Useful when most PSI values are extremely small (< 0.02) or big (> 0.98). Default is \code{FALSE}.
#' @param cell.group.colors Character string. Vector of colors for the cell groups specified for PCA analysis using \code{cell.type.columns}, \code{cell.type.variable}, and \code{cell.type.labels}. If not specified, default \code{ggplot2} colors will be used.
#'
#' @return An object of class S3 with new slot \code{MarvelObject$adhocPlot$PSI}.
#'
#' @importFrom plyr join
#' @import stats
#' @import methods
#' @import ggplot2
#' @import scales
#' @importFrom grDevices hcl
#'
#' @export

PlotValues.PSI <- function(MarvelObject, cell.group.list, feature, maintitle="gene_short_name", xlabels.size=8, max.cells.jitter=10000, max.cells.jitter.seed=1, min.cells=25, sigma.sq=0.001, bimodal.adjust=TRUE, seed=1, modality.column="modality.bimodal.adj", scale.y.log=FALSE, cell.group.colors=NULL) {
    
    # Define arguments
    df <- do.call(rbind.data.frame, MarvelObject$PSI)
    df.pheno <- MarvelObject$SplicePheno
    df.feature <- do.call(rbind.data.frame, MarvelObject$SpliceFeatureValidated)
    cell.group.list <- cell.group.list
    feature <- feature
    maintitle <- maintitle
    xlabels.size <- xlabels.size
    min.cells <- min.cells
    sigma.sq <- sigma.sq
    bimodal.adjust <- bimodal.adjust
    seed <- seed
    modality.column <- modality.column
    scale.y.log <- scale.y.log
    cell.group.colors <- cell.group.colors
    
    # Example arguments
    #MarvelObject <- marvel
    #df <- do.call(rbind.data.frame, MarvelObject$PSI)
    #df.pheno <- MarvelObject$SplicePheno
    #df.feature <- do.call(rbind.data.frame, MarvelObject$SpliceFeatureValidated)
    #cell.group.list <- cell.group.list
    #feature=tran_ids[1]
    #maintitle="gene_short_name"
    #xlabels.size=8.5
    #max.cells.jitter=10000
    #max.cells.jitter.seed=1
    #min.cells=0
    #sigma.sq=0.001
    #bimodal.adjust=TRUE
    #seed=1
    #modality.column="modality.bimodal.adj"
    #scale.y.log <- TRUE
    #cell.group.colors <- c("white", "grey")
    
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
    
    for(i in 1:length(cell.group.list)) {
                
        . <- AssignModality(MarvelObject=s3,
                            sample.ids=cell.group.list[[i]],
                            min.cells=min.cells,
                            sigma.sq=sigma.sq,
                            bimodal.adjust=TRUE,
                            seed=seed
                            )
                       
        mod[i] <- .$Modality$Results[, modality.column]
        
    
    }
    
    modality <- data.frame("cell.type.label"=names(cell.group.list), "modality"=mod, stringsAsFactors=FALSE)
    
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
    
    # Annotate group labels
    df.small <- join(df.small, md, by="sample.id", type="left")
    
    # Remove un-defined samples
    df.small <- df.small[!is.na(df.small$cell.type.label), ]
    
    # Remove missing psi (low coverage)
    df.small <- df.small[!is.na(df.small$psi), ]
    
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
        
    ######################## DOWNSAMPLE JITTER DATA POINTS #############################
    
    # Add ids to trace entries
    df.small$id <- c(1:nrow(df.small))
    
    # Downsample data points for jitter plot
    df.small$psi.downsampled <- df.small$psi
    set.seed(max.cells.jitter.seed)
    
    cell.type.labels <- levels(df.small$cell.type.label)
    
    .list <- list()
    
    for(i in 1:length(cell.type.labels)) {
        
        . <- df.small[which(df.small$cell.type.label==names(cell.group.list)[i]), ]
                    
        if(nrow(.) > max.cells.jitter) {
            
            set.seed(max.cells.jitter.seed)
            
            index <- sample(1:nrow(.), size=max.cells.jitter)
            
            .$psi.downsampled[-index] <- NA
        
            .list[[i]] <- .
            
        } else {
            
            .list[[i]] <- .
            
        }
        
    }
    
    df.small <- do.call(rbind.data.frame, .list)
    df.small <- df.small[order(df.small$id), , drop=FALSE]
    df.small$id <- NULL
    
    #################################### PLOT ######################################
    
    # Definition
    data <- df.small
    x <- data$cell.type.label
    y <- data$psi * 100
    z <- data$cell.type.label
    maintitle <- df.feature.small[, maintitle]
    ytitle <- "PSI"
    xtitle <- ""
    xlabels <- n.cells$label

    data.2 <- df.small[which(!is.na(df.small$psi.downsampled)), ]
    x.jitter <- data.2$cell.type
    y.jitter <- data.2$psi.downsampled * 100
    
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
    if(scale.y.log==FALSE) {
        
        plot <- ggplot() +
            geom_violin(data, mapping=aes(x=x, y=y, fill=z), color="gray", scale="width") +
            geom_jitter(data.2, mapping=aes(x=x.jitter, y=y.jitter), position=position_jitter(width=0.1, height=0), size=0.001) +
            stat_summary(data, mapping=aes(x=x, y=y), geom="point", fun="mean", fill="red", col="black", size=2, shape=23) +
            scale_fill_manual(values=cols) +
            scale_x_discrete(labels=xlabels) +
            scale_y_continuous(breaks=seq(0, 100, by=25), limits=c(0, 100)) +
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
            
    } else if(scale.y.log==TRUE){
        
        # Check skew direction
        y.mean <- mean(y, na.rm=TRUE)
        
        if(y.mean < 50) {
            
            y <- (y + 1)
            y[y>100] <- 100
            
            y.jitter <- (y.jitter + 1)
            y.jitter[y.jitter>100] <- 100
                            
            ytitle <- "PSI (log10 Scale)"
            
            plot <- ggplot() +
                geom_violin(data, mapping=aes(x=x, y=y, fill=z), color="gray", scale="width") +
                geom_jitter(data.2, mapping=aes(x=x.jitter, y=y.jitter), position=position_jitter(width=0.1, height=0), size=0.001) +
                stat_summary(data, mapping=aes(x=x, y=y), geom="point", fun="mean", fill="red", col="black", size=2, shape=23) +
                scale_fill_manual(values=cols) +
                scale_x_discrete(labels=xlabels) +
                scale_y_log10(breaks=c(1, 2.5, 5, 10, 100), labels=c(0, 2.5, 5, 10, 100), limits=c(1, 100)) +
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
                            
        } else if(y.mean >= 0.5){
            
            y <- 100 - y
            y <- (y + 1)
            y[y>100] <- 100
            
            y.jitter <- 100 - y.jitter
            y.jitter <- (y.jitter + 1)
            y.jitter[y.jitter>100] <- 100
                            
            ytitle <- "PSI (log10 Scale)"
            
            reverselog_trans <- function(base = exp(1)) {
              trans <- function(x) -log(x, base)
              inv <- function(x) base^(-x)
              trans_new(paste0("reverselog-", format(base)), trans, inv,
                        log_breaks(base = base),
                        domain = c(1e-100, Inf))
            }
            
            plot <- ggplot() +
                geom_violin(data, mapping=aes(x=x, y=y, fill=z), color="gray", scale="width") +
                geom_jitter(data.2, mapping=aes(x=x.jitter, y=y.jitter), position=position_jitter(width=0.1, height=0), size=0.001) +
                stat_summary(data, mapping=aes(x=x, y=y), geom="point", fun="mean", fill="red", col="black", size=2, shape=23) +
                scale_fill_manual(values=cols) +
                scale_x_discrete(labels=xlabels) +
                scale_y_continuous(trans = reverselog_trans(10), breaks=c(1, 2.5, 5, 10, 100), labels=c(100, 97.5, 95, 90, 0), limits=c(100,1)) +
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
        
    }

    # Save into new slot
    MarvelObject$adhocPlot$PSI <- plot
            
    # Return final object
    return(MarvelObject)
    
}
