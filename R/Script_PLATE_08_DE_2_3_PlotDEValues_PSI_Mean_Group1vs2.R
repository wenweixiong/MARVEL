#' @title Plot differential splicing analysis results based on mean PSI difference
#'
#' @description Scatterplot of differential splicing analysis results based on mean PSI difference between 2 groups of cells. x-axis represents the mean PSI values of cell group 1. y-axis represents the mean PSI values of cell group 2.
#'
#' @param MarvelObject Marvel object. S3 object generated from \code{CompareValues} function.
#' @param method Character string. The statistical method used for differential splicing analysis.
#' @param pval Numeric value. Adjusted p-value below which the splcing event are considered as statistically significant and will consequently be color-annotated on the plot.
#' @param delta Numeric value. The positive (and negative) value specified above (and below) which the splicing events are considered to be statistically significant and will consequently be color-annotated on the plot.
#' @param point.size Numeric value. The point size for the data points. Default value is \code{1}.
#' @param anno Logical value. If set to \code{TRUE}, the specific gene names will be annotated on the plot. Speficified together with \code{anno.tran_id}.
#' @param anno.tran_id Vector of character strings. When \code{anno} set to \code{TRUE}, the coordinates of the splicing events to be annotated on the plot.
#' @param label.size Numeric value. Only applicable if \code{anno} set to TRUE. Size of the gene name labels.
#' @param xlabel.size Numeric value. Font size of the xtick labels. Default is \code{8}.
#' @param point.alpha Numeric value. Transpancy of data points. Default is \code{1}.
#' @param event.types Vector of character string(s). The specific splicing event to plot. May take any one or more of the following values \code{"SE"}, \code{"MXE"}, \code{"RI"}, \code{"A5SS"}, \code{"A3SS"}, \code{"AFE"}, and  \code{"ALE"}.
#' @param event.types.colors Vector of character string(s). Customise colors as per splicing event type specified in \code{event.types} option. Should be of same length as \code{event.types} option.
#'
#' @return An object of class S3 containing with new slot \code{MarvelObject$DE$PSI$Plot[["method"]]}.
#'
#' @importFrom plyr join
#' @import ggplot2
#' @import scales
#'
#' @export
#'
#' @examples
#' marvel.demo <- readRDS(system.file("extdata/data", "marvel.demo.rds", package="MARVEL"))
#'
#' marvel.demo <- PlotDEValues.PSI.Mean.g2vsg1(MarvelObject=marvel.demo,
#'                                      method="ad",
#'                                      pval=0.10,
#'                                      delta=5
#'                                      )
#'
#' # Check output
#' marvel.demo$DE$PSI$Plot
#' marvel.demo$DE$PSI$Summary

PlotDEValues.PSI.Mean.g2vsg1 <- function(MarvelObject, method, pval, delta=5, point.size=1, xlabel.size=8, anno=FALSE, anno.tran_id=NULL, label.size=2.5, point.alpha=1.0, event.types=c("SE", "MXE", "RI", "A5SS", "A3SS", "AFE", "ALE"), event.types.colors=NULL) {

    # Define arguments
    MarvelObject <- MarvelObject
    method <- method
    df <- MarvelObject$DE$PSI$Table[[method]]
    pval <- pval
    delta <- delta
    anno <- anno
    anno.tran_id <- anno.tran_id
    label.size <- label.size
    point.size <- point.size
    label.size <- label.size
    xlabel.size <- xlabel.size
    event.types <- event.types
    event.types.colors <- event.types.colors
    
    # Example arguments
    #MarvelObject <- marvel
    #method <- c("ad", "dts")
    #df <- MarvelObject$DE$PSI$Table[[method]]
    #pval <- c(0.10, 0.10)
    #delta <- 10
    #anno <- FALSE
    #anno.tran_id <- c("chr6:90561622:90561667:-@chr6:90560214|90560234:90560076", "chrX:119629318:119629508:-@chrX:119625379|119625396:119625335")
    #label.size <- 2.5
    #point.size <- 1
    #point.alpha <- 1
    #event.types <- c("SE", "MXE", "RI", "A5SS", "A3SS", "AFE", "ALE")
    #event.types.colors <- event.types.colors
    
    ##############################################
    
    # Tabulate sig events
    .list <- list()
    
    for(i in 1:length(method)) {
    
        # Subset relevent splicing DE results
        de.psi <- MarvelObject$DE$PSI$Table[[method[i]]]
        
        # Subset sig events
        index <- which(abs(de.psi$mean.diff) > delta & de.psi$p.val.adj < pval[i] & de.psi$outlier==FALSE)
        de.psi <- de.psi[index, ]
        
        # Subset gene metadata
        cols <- c("tran_id", "event_type", "gene_id", "gene_short_name", "mean.g1", "mean.g2")
        de.psi <- de.psi[, cols]
        
        # Save into list
        .list[[i]] <- de.psi
        
    }
    
    df <- do.call(rbind.data.frame, .list)
    df <- unique(df)
    
    df$sig <- ifelse(df$mean.g2 > df$mean.g1, "up", "down")
    
    # Append non-sig events
    de.psi.2 <- MarvelObject$DE$PSI$Table[[method[1]]]
    de.psi.2 <- de.psi.2[-which(de.psi.2$tran_id %in% df$tran_id), ]
    cols <- c("tran_id", "event_type", "gene_id", "gene_short_name", "mean.g1", "mean.g2")
    de.psi.2 <- de.psi.2[, cols]
    de.psi.2$sig <- "n.s."
    
    df <- rbind.data.frame(df, de.psi.2)
    
    # Subset relevant event types
    df <- df[which(df$event_type %in% event.types), ]
    
    # Set factor levels
    df$sig <- factor(df$sig, levels=c("up", "down", "n.s."))
    df <- df[order(df$sig, decreasing=TRUE), ]
    
    # Indicate color scheme
    sig.up <- which(df$sig=="up")
    sig.down <- which(df$sig=="down")

    if(length(sig.up) != 0 & length(sig.down) != 0) {

        col.breaks <- c("red", "blue", "gray")
        
    } else if(length(sig.up) != 0 & length(sig.down) == 0) {

        col.breaks <- c("red", "gray")
        
    } else if(length(sig.up) == 0 & length(sig.down) != 0) {

        col.breaks <- c("blue", "gray")

    } else if(length(sig.up) == 0 & length(sig.down) == 0) {

        col.breaks <- "gray"
        
    }
    
    # Create labels
    if(anno==TRUE) {
        
        if(is.null(event.types.colors[1])) {
            
            # Create labels
            df$label <- NA
            df$gene_short_name.event_type <- paste(df$gene_short_name, " (", df$event_type, ")", sep="")
                        
            index <- which(df$tran_id %in% anno.tran_id)
            df$label[index] <- df$gene_short_name.event_type[index]
            
            # Definition
            data <- df
            x <- data$mean.g1
            y <- data$mean.g2
            z <- data$sig
            label <- data$label
            maintitle <- ""
            xtitle <- "Cell group 1 (mean PSI)"
            ytitle <- "Cell group 2 (mean PSI)"
            legendtitle <- ""
            
            xmin <- 0 ; xmax <- 100 ; xinterval <- 25
            ymin <- 0 ; ymax <- 100 ; yinterval <- 25
            
            # Plot
            plot <- ggplot() +
                geom_point(data, mapping=aes(x=x, y=y, color=z), size=point.size, alpha=point.alpha) +
                geom_abline(intercept=c(delta*-1, delta), size=0.25, linetype="dashed", color="black") +
                ggrepel::geom_text_repel(data, mapping=aes(x=x, y=y, label=label), max.overlaps = Inf, box.padding = 0.5, size=label.size, max.time = 1, max.iter = 1e5, segment.alpha=0.5, segment.size=0.1, min.segment.length = 0) +
                scale_color_manual(values=col.breaks) +
                #ggrepel::geom_text_repel(data, mapping=aes(x=x, y=y, label=label), max.overlaps = Inf, box.padding = 1.0, size=2, max.time = 1, max.iter = 1e5, segment.alpha=0.5, segment.size=0.1, min.segment.length = 0) +
                scale_x_continuous(breaks=seq(xmin, xmax, by=xinterval), limits=c(xmin, xmax)) +
                scale_y_continuous(breaks=seq(ymin, ymax, by=yinterval), limits=c(ymin, ymax)) +
                labs(title=maintitle, x=xtitle, y=ytitle, color=legendtitle) +
                theme(panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.background = element_blank(),
                      panel.border=element_blank(),
                      plot.title=element_text(size=12, hjust=0.5),
                      axis.line=element_line(colour = "black"),
                      axis.title=element_text(size=10),
                      axis.text=element_text(size=10, colour="black"),
                      legend.position="none",
                      legend.title=element_text(size=8),
                      legend.text=element_text(size=8),
                      legend.key = element_blank()
                      )
                      
        } else {
            
            # Define color scheme
                # Create reference table for color scheme
                event.type.color.df <- data.frame("event.types"=c(event.types, "n.s."),
                                                  "event.types.colors"=c(event.types.colors, "gray"),
                                                  stringsAsFactors=FALSE
                                                  )
                
                # Censor n.s. events
                df$event_type[which(df$sig=="n.s.")] <- "n.s."
                
                # Subset events found
                event.type.color.df <- event.type.color.df[which(event.type.color.df$event.types %in% unique(df$event_type)), ]
                
                # Define final colors
                df$event_type <- factor(df$event_type, levels=event.type.color.df$event.types)
                col.breaks <- event.type.color.df$event.types.colors
                
                # Reorder for asthetic purpose
                df <- df[order(df$event_type, decreasing=TRUE), ]
                
            # Create labels
            df$label <- NA
            df$gene_short_name.event_type <- paste(df$gene_short_name, " (", df$event_type, ")", sep="")
                        
            index <- which(df$tran_id %in% anno.tran_id)
            df$label[index] <- df$gene_short_name.event_type[index]
            
            # Definition
            data <- df
            x <- data$mean.g1
            y <- data$mean.g2
            z <- data$event_type
            label <- data$label
            maintitle <- ""
            xtitle <- "Cell group 1 (mean PSI)"
            ytitle <- "Cell group 2 (mean PSI)"
            legendtitle <- ""
            
            xmin <- 0 ; xmax <- 100 ; xinterval <- 25
            ymin <- 0 ; ymax <- 100 ; yinterval <- 25
            
            # Plot
            plot <- ggplot() +
                geom_point(data, mapping=aes(x=x, y=y, color=z), size=point.size, alpha=point.alpha) +
                geom_abline(intercept=c(delta*-1, delta), size=0.25, linetype="dashed", color="black") +
                ggrepel::geom_text_repel(data, mapping=aes(x=x, y=y, label=label), max.overlaps = Inf, box.padding = 0.5, size=label.size, max.time = 1, max.iter = 1e5, segment.alpha=0.5, segment.size=0.1, min.segment.length = 0) +
                scale_color_manual(values=col.breaks) +
                #ggrepel::geom_text_repel(data, mapping=aes(x=x, y=y, label=label), max.overlaps = Inf, box.padding = 1.0, size=2, max.time = 1, max.iter = 1e5, segment.alpha=0.5, segment.size=0.1, min.segment.length = 0) +
                scale_x_continuous(breaks=seq(xmin, xmax, by=xinterval), limits=c(xmin, xmax)) +
                scale_y_continuous(breaks=seq(ymin, ymax, by=yinterval), limits=c(ymin, ymax)) +
                labs(title=maintitle, x=xtitle, y=ytitle, color=legendtitle) +
                theme(panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.background = element_blank(),
                      panel.border=element_blank(),
                      plot.title=element_text(size=12, hjust=0.5),
                      axis.line=element_line(colour = "black"),
                      axis.title=element_text(size=10),
                      axis.text=element_text(size=10, colour="black"),
                      #legend.position="none",
                      legend.title=element_text(size=8),
                      legend.text=element_text(size=8),
                      legend.key = element_blank()
                      ) +
                guides(color = guide_legend(override.aes=list(size=2, ncol=1)))
                        
        }
                
    ########################################################################
    
    } else {
        
        if(is.null(event.types.colors[1])) {

            # Definition
            data <- df
            x <- data$mean.g1
            y <- data$mean.g2
            z <- data$sig
            maintitle <- ""
            xtitle <- "Cell group 1 (mean PSI)"
            ytitle <- "Cell group 2 (mean PSI)"
            legendtitle <- ""
            
            xmin <- 0 ; xmax <- 100 ; xinterval <- 25
            ymin <- 0 ; ymax <- 100 ; yinterval <- 25
            
            # Plot
            plot <- ggplot() +
                geom_point(data, mapping=aes(x=x, y=y, color=z), size=point.size, alpha=point.alpha) +
                geom_abline(intercept=c(delta*-1, delta), size=0.25, linetype="dashed", color="black") +
                scale_color_manual(values=col.breaks) +
                #ggrepel::geom_text_repel(data, mapping=aes(x=x, y=y, label=label), max.overlaps = Inf, box.padding = 1.0, size=2, max.time = 1, max.iter = 1e5, segment.alpha=0.5, segment.size=0.1, min.segment.length = 0) +
                scale_x_continuous(breaks=seq(xmin, xmax, by=xinterval), limits=c(xmin, xmax)) +
                scale_y_continuous(breaks=seq(ymin, ymax, by=yinterval), limits=c(ymin, ymax)) +
                labs(title=maintitle, x=xtitle, y=ytitle, color=legendtitle) +
                theme(panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.background = element_blank(),
                      panel.border=element_blank(),
                      plot.title=element_text(size=12, hjust=0.5),
                      axis.line=element_line(colour = "black"),
                      axis.title=element_text(size=10),
                      axis.text=element_text(size=10, colour="black"),
                      legend.position="none",
                      legend.title=element_text(size=8),
                      legend.text=element_text(size=8),
                      legend.key = element_blank()
                      )
                  
        } else {
            
            # Define color scheme
                # Create reference table for color scheme
                event.type.color.df <- data.frame("event.types"=c(event.types, "n.s."),
                                                  "event.types.colors"=c(event.types.colors, "gray"),
                                                  stringsAsFactors=FALSE
                                                  )
                
                # Censor n.s. events
                df$event_type[which(df$sig=="n.s.")] <- "n.s."
                
                # Subset events found
                event.type.color.df <- event.type.color.df[which(event.type.color.df$event.types %in% unique(df$event_type)), ]
                
                # Define final colors
                df$event_type <- factor(df$event_type, levels=event.type.color.df$event.types)
                col.breaks <- event.type.color.df$event.types.colors
                
                # Reorder for asthetic purpose
                df <- df[order(df$event_type, decreasing=TRUE), ]
                
            # Definition
            data <- df
            x <- data$mean.g1
            y <- data$mean.g2
            z <- data$event_type
            maintitle <- ""
            xtitle <- "Cell group 1 (mean PSI)"
            ytitle <- "Cell group 2 (mean PSI)"
            legendtitle <- ""
            
            xmin <- 0 ; xmax <- 100 ; xinterval <- 25
            ymin <- 0 ; ymax <- 100 ; yinterval <- 25
            
            # Plot
            plot <- ggplot() +
                geom_point(data, mapping=aes(x=x, y=y, color=z), size=point.size, alpha=point.alpha) +
                geom_abline(intercept=c(delta*-1, delta), size=0.25, linetype="dashed", color="black") +
                scale_color_manual(values=col.breaks) +
                #ggrepel::geom_text_repel(data, mapping=aes(x=x, y=y, label=label), max.overlaps = Inf, box.padding = 1.0, size=2, max.time = 1, max.iter = 1e5, segment.alpha=0.5, segment.size=0.1, min.segment.length = 0) +
                scale_x_continuous(breaks=seq(xmin, xmax, by=xinterval), limits=c(xmin, xmax)) +
                scale_y_continuous(breaks=seq(ymin, ymax, by=yinterval), limits=c(ymin, ymax)) +
                labs(title=maintitle, x=xtitle, y=ytitle, color=legendtitle) +
                theme(panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.background = element_blank(),
                      panel.border=element_blank(),
                      plot.title=element_text(size=12, hjust=0.5),
                      axis.line=element_line(colour = "black"),
                      axis.title=element_text(size=10),
                      axis.text=element_text(size=10, colour="black"),
                      #legend.position="none",
                      legend.title=element_text(size=8),
                      legend.text=element_text(size=8),
                      legend.key = element_blank()
                      )  +
                guides(color = guide_legend(override.aes=list(size=2, ncol=1)))
            
        }
    
    }
    
    ########################################################################

    # Summary
    tbl <- as.data.frame(table(df$sig))
    names(tbl) <- c("sig", "freq")
    
    ##############################################
    
    # Save to new slot
    MarvelObject$DE$PSI$Plot <- plot
    MarvelObject$DE$PSI$Summary <- tbl
    
    # Return final object
    return(MarvelObject)
        
}
