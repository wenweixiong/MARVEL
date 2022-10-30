#' @title Tabulate differentially spliced splicing event
#'
#' @description Tabulates the percentage or absoluate number of significant splicing events for each splicing type.
#'
#' @param MarvelObject Marvel object. S3 object generated from \code{CompareValues} function.
#' @param method Character string. The statistical method used for differential splicing analysis.
#' @param psi.pval Numeric value. Adjusted p-value below which the splicing event is considered differentially spliced and included for tabulation.
#' @param psi.mean.diff Numeric value. The minimum absolute differences in PSI values between the two cell groups above which the splicing event is considered differentially spliced and included for tabulation.
#' @param ylabels.size Numeric value. Size of the xtick labels. Default is \code{8}.
#' @param barlabels.size Numeric value. Size of the labels above each bar. Default is \code{3}
#' @param x.offset Numeric value. The values on the x-axis to offset by. Useful when right margin overshadow the numbers above the bars. Default value is \code{0}.
#' @param direction.color Character strings. Vector of length 2 to specify the colors for significanly down- and up-regulated splicing events. Default is \code{NULL}, which corresponds to default \code{ggplot2} color scheme.
#' @param mode Character strings. When set to \code{"percentage"} (default), percentage of significant splicing events over total splicing events detected will be tabulate. When set to \code{absolute}, the number of significant splicing events will be tabulated.
#'
#' @return An object of class S3 with new slots \code{MarvelObject$DE$AbsASE$Table} and \code{MarvelObject$DE$AbsASE$Plot}.
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
#' marvel.demo <- PctASE(MarvelObject=marvel.demo,
#'                       method="ad",
#'                       psi.pval=0.1,
#'                       psi.mean.diff=0
#'                       )

PctASE <- function(MarvelObject, method, psi.pval, psi.mean.diff, ylabels.size=8, barlabels.size=3, x.offset=0, direction.color=NULL, mode="percentage") {
    
    # Define arguments
    method <- method
    psi.pval <- psi.pval
    psi.mean.diff <- psi.mean.diff
    ylabels.size <- ylabels.size
    barlabels.size <- barlabels.size
    x.offset <- x.offset
    direction.color <- direction.color
    mode <- mode
    
    # Example arguments
    #MarvelObject <- marvel
    #method <- c("ad", "dts")
    #psi.pval <- c(0.10, 0.10)
    #psi.mean.diff <- 0
    #ylabels.size <- 8
    #barlabels.size <- 3
    #mode <- "percentage"
    #x.offset <- 1
    #direction.color <- c("blue", "red3")
    
    # Subset sig events
    results.list <- list()
    
    for(i in 1:length(method)) {
        
        # Retrieve DE result table
        results <- MarvelObject$DE$PSI$Table[[method[i]]]
        
        # Subset sig events
        index.up <- which(results$p.val.adj < psi.pval[i] & results$mean.diff > psi.mean.diff & results$outliers==FALSE)
        index.down <- which(results$p.val.adj < psi.pval[i] & results$mean.diff < (psi.mean.diff * -1) & results$outliers==FALSE)
        results.small <- results[c(index.up, index.down), ]
        
        # Subset relevant columns
        cols <- c("tran_id", "event_type", "mean.diff")
        results.small <- results.small[, cols]
        
        # Save into list
        results.list[[i]] <- results.small
        
    }
    
    results <- unique(do.call(rbind.data.frame, results.list))
    results <- unique(results)
    results$direction <- ifelse(results$mean.diff > psi.mean.diff, "up", "down")
    results.sig <- results
    
    # Retrieve all expressed events
    results <- MarvelObject$DE$PSI$Table[[method[1]]]
    
    # Annotate sig events
    results <- join(results, results.sig[,c("tran_id", "direction")], by="tran_id", type="left")
    
    ################################################################################
    
    # Define events included for DE analysis
    event_types <- intersect(c("SE", "MXE", "RI", "A5SS", "A3SS", "AFE", "ALE"), unique(results$event_type))
    
    # Tabulate n sig, total events
    n.total <- NULL
    n.sig.up <- NULL
    n.sig.down <- NULL

    for(i in 1:length(event_types)) {

        # Subset event type
        results.small <- results[which(results$event_type==event_types[i]), ]
        
        # Tabulate n
        n.total[i] <- nrow(results.small)
        n.sig.up[i] <- sum(results.small$direction=="up", na.rm=TRUE)
        n.sig.down[i] <- sum(results.small$direction=="down", na.rm=TRUE)

    }

    results <- data.frame("event_type"=event_types,
                          "n.total"=n.total,
                          "n.sig.up"=n.sig.up,
                          "pct.sig.up"= round((n.sig.up/n.total) * 100, digits=1),
                          "n.sig.down"=n.sig.down,
                          "pct.sig.down"= round((n.sig.down/n.total) * 100, digits=1),
                          stringsAsFactors=FALSE
                          )

    # Set factor levels
    results$event_type <- factor(results$event_type, levels=rev(event_types))
    results <- results[order(results$event_type), ]
    
    # Save into new slots
    MarvelObject$DE$AbsASE$Table <- results
      
    # Indicate n total events in x-labels
    xlabels <- paste(results$event_type, " (n=", results$n.total, ")", sep="")

    # Format data frame for barchat
        # %
        results.small <- results[,c("event_type", "pct.sig.up", "pct.sig.down")]
        results.small <- reshape2::melt(data=results.small, id.vars="event_type", measure.vars=c("pct.sig.up", "pct.sig.down"))
        names(results.small) <- c("event_type", "direction", "pct")
        results.small$direction <- gsub("pct.", "", results.small$direction, fixed=TRUE)
        results.small.1 <- results.small
        
        # n
        results.small <- results[,c("event_type", "n.sig.up", "n.sig.down")]
        results.small <- reshape2::melt(data=results.small, id.vars="event_type", measure.vars=c("n.sig.up", "n.sig.down"))
        names(results.small) <- c("event_type", "direction", "n")
        results.small$direction <- gsub("n.", "", results.small$direction, fixed=TRUE)
        results.small.2 <- results.small
        
        # Merge
        results.small.1$id <- paste(results.small.1$event_type, "_", results.small.1$direction, sep="")
        results.small.2$id <- paste(results.small.2$event_type, "_", results.small.2$direction, sep="")
        results <- join(results.small.1, results.small.2[,c("id", "n")], by="id", type="left")
        results$id <- NULL
        
    # Set factor levels
    results$direction <- factor(results$direction, levels=c("sig.down", "sig.up"), labels=c("down", "up"))

    if(sum(results$pct)==0) {
        
        message(paste("No splicing events statistcally significant at psi.pval and/or psi.mean.diff threshold specified"))
        
        return(MarvelObject)
        
    }
    
    ###########################################################################
        
    if(mode=="percentage") {
        
        # Barplot
            # Definition
            data <- results
            x <- data$event_type
            y <- data$pct
            y2 <- data$n
            z <- data$direction
            maintitle <- ""
            xtitle <- ""
            ytitle <- "% Sig ASE"
            legendtitle <- "Direction"

            ymin <- 0 ; ymax <- max(y) + x.offset
            
            # Color scheme
            if(is.null(direction.color[1])) {
            
                gg_color_hue <- function(n) {
                  hues = seq(15, 375, length = n + 1)
                  hcl(h = hues, l = 65, c = 100)[1:n]
                }
                
                n = length(levels(z))
                cols = gg_color_hue(n)
            
            } else {
                
                cols <- direction.color
                
            }
            
            # Plot
            plot <- ggplot() +
                geom_bar(data=data, aes(x=x, y=y, fill=z), stat="identity", color="black", position=position_dodge(), width=0.8) +
                geom_text(data=data, mapping=aes(x=x, y=y, fill=z, label=y2), position=position_dodge(width=1.05), vjust=0.25, hjust=-0.25, size=barlabels.size) +
                scale_y_continuous(limits=c(ymin, ymax)) +
                scale_fill_manual(values=cols) +
                scale_x_discrete(labels=xlabels) +
                labs(title=maintitle, x=xtitle, y=ytitle, fill=legendtitle) +
                theme(panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.background = element_blank(),
                    panel.border=element_blank(),
                    plot.title=element_text(hjust = 0.5, size=12),
                    plot.subtitle=element_text(hjust = 0.5, size=12),
                    axis.line.y.left = element_line(color="black"),
                    axis.line.x = element_line(color="black"),
                    axis.title=element_text(size=12),
                    axis.text=element_text(size=12),
                    axis.text.x=element_text(size=10, colour="black"),
                    axis.text.y=element_text(size=ylabels.size, colour="black"),
                    legend.title=element_text(size=8),
                    legend.text=element_text(size=8)
                    ) +
                coord_flip()
         
        # Retrieve pct
        results.small <- results[,c("event_type", "direction", "pct")]
        results.small <- results.small[order(results.small$event_type, decreasing=TRUE, results.small$direction), ]
        
        # Save into new slots
        MarvelObject$DE$PctASE$Table <- results.small
        MarvelObject$DE$PctASE$Plot <- plot
        
    } else if(mode=="absolute") {
        
        # Barplot
            # Definition
            data <- results
            x <- data$event_type
            y <- data$n
            y2 <- data$n
            z <- data$direction
            maintitle <- ""
            xtitle <- ""
            ytitle <- "ASE (n)"
            legendtitle <- "Direction"

            ymin <- 0 ; ymax <- max(y) + x.offset
       
           # Color scheme
           if(is.null(direction.color[1])) {
           
               gg_color_hue <- function(n) {
                 hues = seq(15, 375, length = n + 1)
                 hcl(h = hues, l = 65, c = 100)[1:n]
               }
               
               n = length(levels(z))
               cols = gg_color_hue(n)
           
           } else {
               
               cols <- direction.color
               
           }
            # Plot
            plot <- ggplot() +
                geom_bar(data=data, aes(x=x, y=y, fill=z), stat="identity", color="black", position=position_dodge(), width=0.8) +
                geom_text(data=data, mapping=aes(x=x, y=y, fill=z, label=y2), position=position_dodge(width=1.05), vjust=0.25, hjust=-0.25, size=barlabels.size) +
                scale_fill_manual(values=cols) +
                scale_y_continuous(limits=c(ymin, ymax)) +
                scale_x_discrete(labels=xlabels) +
                labs(title=maintitle, x=xtitle, y=ytitle, fill=legendtitle) +
                theme(panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.background = element_blank(),
                    panel.border=element_blank(),
                    plot.title=element_text(hjust = 0.5, size=12),
                    plot.subtitle=element_text(hjust = 0.5, size=12),
                    axis.line.y.left = element_line(color="black"),
                    axis.line.x = element_line(color="black"),
                    axis.title=element_text(size=12),
                    axis.text=element_text(size=12),
                    axis.text.x=element_text(size=10, colour="black"),
                    axis.text.y=element_text(size=ylabels.size, colour="black"),
                    legend.title=element_text(size=8),
                    legend.text=element_text(size=8)
                    ) +
                coord_flip()
         
        # Retrieve n
        results.small <- results[,c("event_type", "direction", "n")]
        results.small <- results.small[order(results.small$event_type, decreasing=TRUE, results.small$direction), ]
         
        # Save into new slots
        MarvelObject$DE$AbsASE$Plot <- plot
        
    }
    
    return(MarvelObject)
        
}
