#' @title Tabulate the number of expressed splicing events
#'
#' @description Tabulates and plots the number of expressed splicing events for each splicing event category for a specified cell group.
#'
#' @param MarvelObject Marvel object. S3 object generated from \code{TransformExpValues} function.
#' @param sample.ids Vector of character strings. Sample IDs that constitute the cell group.
#' @param min.cells Numeric value. Minimum number of cells expressing the splicing event for the event to be included for tabulation. A splicing event is defined as expressed when it has a non-missing PSI value.
#' @param event.group.colors Vector of character strings. Colors for the event groups. If not specified, default \code{ggplot2} colors will be used.
#'
#' @return An object of class S3 with new slots \code{MarvelObject$N.Events$Table} and \code{MarvelObject$N.Events$Plot}.
#'
#' @importFrom plyr join
#' @import methods
#' @import ggplot2
#' @import scales
#' @importFrom grDevices hcl
#'
#' @export
#'
#' @examples
#' marvel.demo <- readRDS(system.file("extdata/data", "marvel.demo.rds", package="MARVEL"))
#'
#' # Define cell group for analysis
#' df.pheno <- marvel.demo$SplicePheno
#' sample.ids <- df.pheno[which(df.pheno$cell.type=="iPSC"), "sample.id"]
#'
#' # Tabulate expressed events
#' marvel.demo <- CountEvents(MarvelObject=marvel.demo,
#'                            sample.ids=sample.ids,
#'                            min.cells=5,
#'                            event.group.colors=NULL
#'                            )
#'
#' # Check outputs
#' marvel.demo$N.Events$Table
#' marvel.demo$N.Events$Plot

CountEvents <- function(MarvelObject, sample.ids, min.cells, event.group.colors=NULL) {

    # Define arguments
    psi <- do.call(rbind.data.frame, MarvelObject$PSI)
    psi.feature <- do.call(rbind.data.frame, MarvelObject$SpliceFeatureValidated)
    psi.pheno <- MarvelObject$SplicePheno
    sample.ids <- sample.ids
    min.cells <- min.cells
    event.group.colors <- event.group.colors
    
    # Example arguments
    #MarvelObject <- marvel
    #psi <- do.call(rbind.data.frame, MarvelObject$PSI)
    #psi.feature <- do.call(rbind.data.frame, MarvelObject$SpliceFeatureValidated)
    #psi.pheno <- MarvelObject$SplicePheno
    #sample.ids <- sample.id
    #min.cells <- 1
    #event.group.colors <- NULL
    
    ####################################################
    
    # Generate row names
    row.names(psi) <- psi$tran_id
    psi$tran_id <- NULL
    
    # Subset relevant cells
    psi.pheno <- psi.pheno[psi.pheno$sample.id %in% sample.ids, , drop=FALSE]
    psi <- psi[, which(names(psi) %in% psi.pheno$sample.id), drop=FALSE]
    
    if(nrow(psi)==0) {
        
        warning("No cells identified")
        return(MarvelObject)
        
    }
    
    # Compute number of expressed events
    event_types <- unique(psi.feature$event_type)
    
    n.events <- NULL
    
    for(i in 1:length(event_types)) {
        
        # Subset relevant event
        psi.feature.small <- psi.feature[which(psi.feature$event_type==event_types[i]), , drop=FALSE]
        psi.small <- psi[psi.feature.small$tran_id, , drop=FALSE]
        
        # Subset expressed events
        . <- apply(psi.small, 1, function(x) {sum(!is.na(x))})
        . <- sum(. >= min.cells)
        n.events[i] <- .
        
    }
    
    results <- data.frame("event_type"=event_types,
                          "freq"=n.events,
                          stringsAsFactors=FALSE
                          )
    
    # Set factor levels
    levels <- c("SE", "MXE", "RI", "A5SS", "A3SS", "AFE", "ALE")
    results$event_type <- factor(results$event_type, levels=levels)

    # Compute proportions
    results$pct <- results$freq / sum(results$freq) * 100
    results <- results[order(results$pct, decreasing=TRUE), ]
    
    # Define color scheme
        # Retrieve default colors
        gg_color_hue <- function(n) {
            hues = seq(15, 375, length = n + 1)
            hcl(h = hues, l = 65, c = 100)[1:n]
        }
            
        n = length(levels)
        cols = gg_color_hue(n)
        
        color.df <- data.frame("event_type"=levels, "color"=cols, stringsAsFactors=FALSE)
        
        # Subset expressed events
        color.df <- color.df[which(color.df$event_type %in% unique(results$event_type)), ]
        
        # Final color scheme
        if(is.null(event.group.colors[1])) {
        
            cols <- color.df$color
        
        } else {
            
            cols <- event.group.colors
            
        }

    # Plot
        # Compute statistics for plot
        results$fraction <- results$freq / sum(results$freq)
        results$ymax <- cumsum(results$fraction)
        results$ymin = c(0, results$ymax[-length(results$ymax)])
        
        # Definitions
        data <- results
        xmax <- nrow(data) + 1
        xmin <- nrow(data)
        ymax <- data$ymax
        ymin <- data$ymin
        z <- data$event_type
        maintitle <- ""
        xtitle <- ""
        ytitle <- ""
        legendtitle <- "Event Type"
        
        # Plot
        plot <- ggplot() +
            geom_rect(data=data, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=z), color="black") +
            coord_polar(theta="y") +
            xlim(c(2, 4)) +
            scale_fill_manual(values=cols) +
            labs(title=maintitle, x=xtitle, y=ytitle, fill=legendtitle) +
            theme(panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.background = element_blank(),
                panel.border=element_blank(),
                plot.title=element_text(hjust = 0.5, size=15),
                plot.subtitle=element_text(hjust = 0.5, size=15),
                axis.line = element_blank(),
                axis.ticks=element_blank(),
                axis.text=element_blank(),
                legend.title=element_text(size=6),
                legend.text=element_text(size=6)
                )
  
    ####################################################
    
    # Update slots
    MarvelObject$N.Events$Table <- results[, c("event_type", "freq", "pct")]
    MarvelObject$N.Events$Plot <- plot
    
    # Return final object
    return(MarvelObject)
            
}
