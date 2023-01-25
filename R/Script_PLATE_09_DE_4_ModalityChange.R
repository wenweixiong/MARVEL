#' @title Classify modality changes
#'
#' @description Classifies the type of modality change for each splicing event that has taken place between 2 groups of cells.
#'
#' @param MarvelObject Marvel object. S3 object generated from \code{CompareValues} function.
#' @param method Character string. The statistical method used for differential splicing analysis.
#' @param psi.pval Numeric value. Adjusted p-value below which the splicing event is considered differentially spliced and included for modality analysis.
#' @param psi.delta Numeric value. The absolute difference between the means PSI values of cell group 1 and 2, above which, the splicing event is considered differentially spliced and included for modality analysis.
#'
#' @return An object of class S3 with new slots \code{MarvelObject$DE$Modality$Table}, \code{MarvelObject$DE$Modality$Plot}, and \code{MarvelObject$DE$Modality$Plot.Stats}.
#'
#' @importFrom plyr join
#' @importFrom stats var
#' @import methods
#' @import ggplot2
#'
#' @export
#'
#' @examples
#' marvel.demo <- readRDS(system.file("extdata/data", "marvel.demo.rds", package="MARVEL"))
#'
#' marvel.demo <- ModalityChange(MarvelObject=marvel.demo,
#'                               method="ad",
#'                               psi.pval=0.1,
#'                               psi.delta=0
#'                               )
#'
#' # Check outputs
#' head(marvel.demo$DE$Modality$Table)
#' marvel.demo$DE$Modality$Plot
#' marvel.demo$DE$Modality$Plot.Stats

ModalityChange <- function(MarvelObject, method, psi.pval, psi.delta=0) {
    
    # Define arguments
    method <- method
    psi.pval <- psi.pval
    psi.delta <- psi.delta
    
    # Example arguments
    #MarvelObject <- marvel
    #method <- c("ad", "dts")
    #psi.pval <- c(0.10, 0.10)
    #psi.delta <- 5
    
    # Annotate modality change
    for(i in 1:length(method)) {
        
        # Retrieve DE result table
        results <- MarvelObject$DE$PSI$Table[[method[i]]]
        
        # Indicate modality change
            # Restricted
            results$modality.change <- NA
            index.na <- which(is.na(results$modality.change))
            index.new <- which(results$modality.bimodal.adj.g1==results$modality.bimodal.adj.g2)
            index <- intersect(index.na, index.new)
            results$modality.change[index] <- "Restricted"
            
            # Implicit
                # Included
                index.na <- which(is.na(results$modality.change))
                index.new <- intersect(grep("^Included", results$modality.bimodal.adj.g1),
                                       grep("^Included", results$modality.bimodal.adj.g2)
                                       )
                index <- intersect(index.na, index.new)
                results$modality.change[index] <- "Implicit"
                
                # Included
                index.na <- which(is.na(results$modality.change))
                index.new <- intersect(grep("^Excluded", results$modality.bimodal.adj.g1),
                                       grep("^Excluded", results$modality.bimodal.adj.g2)
                                       )
                index <- intersect(index.na, index.new)
                results$modality.change[index] <- "Implicit"
        
            # Explicit
            index.na <- which(is.na(results$modality.change))
            index.new <- intersect(which(!is.na(results$modality.bimodal.adj.g1)), which(!is.na(results$modality.bimodal.adj.g2)))
            index <- intersect(index.na, index.new)
            results$modality.change[index] <- "Explicit"
            
            # Save into original slot
            MarvelObject$DE$PSI$Table[[method[i]]] <- results
            
    }
    
    # Subset sig events
    results.list <- list()
    
    for(i in 1:length(method)) {
        
        # Retrieve DE result table
        results <- MarvelObject$DE$PSI$Table[[method[i]]]
        
        # Subset sig events
        index <- which(results$p.val.adj < psi.pval[i] & abs(results$mean.diff) > psi.delta & results$outlier==FALSE)
        results.small <- results[index, ]
        
        # Subset relevant columns
        cols <- c("tran_id", "modality.change")
        results.small <- results.small[, cols]
        
        # Save into list
        results.list[[i]] <- results.small
        
    }
    
    results <- unique(do.call(rbind.data.frame, results.list))

    # Check for missing modality change assignment
    modality.missing <- sum(is.na(results$modality.change))
    
    if(modality.missing==0) {
        
        message("All modality change assignment SUCCESSFUL")
        
    } else {
        
        message("Missing modality assignment, please contact author")
        
    }
    
    # Check of conflicting modality change
    tbl <- as.data.frame(table(results$tran_id))
    freq.max <- max(tbl$Freq)
    
    if(freq.max==1) {
        
        message("Modality change category CONSISTENT across all statistical test")
        
    } else if(freq.max >= 2){
        
        message("Modality change category NOT CONSISTENT across all statistical test")
        
    }
    
    # Doughnut plot
        # Tabulate freq
        . <- as.data.frame(table(results$modality.change), stringsAsFactors=FALSE)
        names(.) <- c("modality.change", "freq")
        .$pct <- .$freq / sum(.$freq) * 100
        
        # Set factor levels
        .$modality.change <- factor(.$modality.change, levels=c("Explicit", "Implicit", "Restricted"))
        . <- .[order(.$modality.change), ]
        
        # Compute statistics for plot
        .$fraction <- .$freq / sum(.$freq)
        .$ymax <- cumsum(.$fraction)
        .$ymin = c(0, .$ymax[-length(.$ymax)])
        
        # Definitions
        data <- .
        xmax <- nrow(data) + 1
        xmin <- nrow(data)
        ymax <- data$ymax
        ymin <- data$ymin
        z <- data$modality.change
        maintitle <- ""
        xtitle <- ""
        ytitle <- ""
        legendtitle <- "Modality Change"
        
        # Plot
        plot <- ggplot() +
            geom_rect(data=data, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=z), color="black") +
            coord_polar(theta="y") +
            xlim(c(2, 4)) +
            #scale_fill_manual(values=colors) +
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
                legend.title=element_text(size=9),
                legend.text=element_text(size=9)
                )
    
    # Annotate event metadata
    df.feature <- MarvelObject$DE$PSI$Table[[method[1]]]
    cols <- c("tran_id", "event_type", "gene_id", "gene_short_name", "gene_type", "modality.bimodal.adj.g1", "modality.bimodal.adj.g2")
    df.feature <- df.feature[,cols]
    
    results <- join(results, df.feature, by="tran_id", type="left")
    col.1 <- "modality.change"
    col.2 <- setdiff(names(results), col.1)
    results <- results[,c(col.2, col.1)]
    
    # Save into new slow
    MarvelObject$DE$Modality$Table <- results
    MarvelObject$DE$Modality$Plot <- plot
    MarvelObject$DE$Modality$Plot.Stats <- .[,c("modality.change", "freq", "pct")]
    
    return(MarvelObject)
        
}
