#' @title Modality Proportion for a Given Event Type(s).
#'
#' @description
#' \code{PropModality.Doughnut} tabulates and plots the proportion of each modality for a specified splicing event type(s).
#'
#' @details
#' This function tabulates and plots the proportion of each modality for a specified splicing event type(s).
#'
#' @param MarvelObject S3 object generated from \code{AssignModality} function.
#' @param modality.column Character string. Can take the value \code{"modality"}, \code{"modality.var"} or \code{"modality.bimodal.adj"}. Please refer to \code{AssignModality} function help page for more details.
#' @param modality.type Character string. \code{basic} indicates that only the main modalities (included, excluded, bimodal, middle, multimodal) are analysed. Sub-modalities (primary and dispersed) will be merged. \code{complete} indicates that both main and sub-modalities are analysed. Sub-modalities will not be merged.
#' @param event.type Character string. To indicate which event type to analyse. Can take the value \code{"SE"}, \code{"MXE"}, \code{"RI"}, \code{"A5SS"} or \code{"A3SS"}. Specify \code{"all"} to include all event types.
#'
#' @export
#'
#' @return An object of class S3 containing all the original slots as inputted by the user in addition to one new slot name \code{$Modality$Prop$DoughnutChart}.
#'
#' @author Sean Wen <sean.wenwx@gmail.com>
#'
#' @importFrom fitdistrplus fitdist
#' @import methods
#' @import ggplot2
#' @importFrom plyr join
#'
#' @examples
#' # Load input
#' marvel <- readRDS(system.file("extdata/Data", "MarvelObject.rds", package="MARVEL"))
#'
#' # Run example
#' marvel <- PropModality.Doughnut(MarvelObject=marvel,
#'                                 modality.column="modality.bimodal.adj",
#'                                 modality.type="basic",
#'                                 event.type=c("SE", "MXE", "RI", "A5SS", "A3SS")
#'                                 )
#'
#' # Check output
#' marvel$Modality$Prop$DoughnutChart$Table
#' marvel$Modality$Prop$DoughnutChart$Plot

PropModality.Doughnut <- function(MarvelObject, modality.column, modality.type, event.type) {

    # Define arguments
    df.feature <- MarvelObject$Modality$Results
    modality.column <- modality.column
    modality.type <- modality.type
    event.type <- event.type
        
    # Subset relevant modality column
    df.feature <- df.feature[,c("event_type", modality.column)]
    names(df.feature)[which(names(df.feature)==modality.column)] <- "modality"
    
    # Subset relevant event type
    df.feature <- df.feature[which(df.feature$event_type %in% event.type), ]
    
    # Merge sub-modalities
    if(modality.type=="basic") {

        df.feature$modality[grep("^Included", df.feature$modality)] <- "Included"
        df.feature$modality[grep("^Excluded", df.feature$modality)] <- "Excluded"
    
    }
    
    # Compute proportions
    . <- as.data.frame(table(df.feature$modality), stringsAsFactors=FALSE)
    names(.) <- c("modality", "freq")
    .$pct <- .$freq / sum(.$freq) * 100
    
    # Set factor levels
    if(modality.type=="basic") {
        
        # Set factor levels
        .$modality <- factor(.$modality, levels=c("Included", "Excluded", "Bimodal", "Middle", "Multimodal"))
        . <- .[order(.$modality), ]
        
        
    } else {
        
        # Set factor levels
        .$modality <- factor(.$modality, levels=c("Included.Primary", "Included.Dispersed", "Excluded.Primary", "Excluded.Dispersed", "Bimodal", "Middle", "Multimodal"))
        . <- .[order(.$modality), ]
        
    }
    
    # Save to new slot
    MarvelObject$Modality$Prop$DoughnutChart$Table <- .

    # Plot
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
        z <- data$modality
        maintitle <- ""
        xtitle <- ""
        ytitle <- ""
        legendtitle <- "Modality"
        
        if(modality.type=="basic") {
            
            colors <- c("red", "blue", "green", "yellow", "purple")
            legend.text.size <- 9
            
        } else {
            
            colors <- c("red", "indianred1", "blue", "deepskyblue1", "green", "yellow", "purple")
            legend.text.size <- 6
            
        }
        
        # Plot
        plot <- ggplot() +
            geom_rect(data=data, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=z), color="black") +
            coord_polar(theta="y") +
            xlim(c(2, 4)) +
            scale_fill_manual(values=colors) +
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
                legend.title=element_text(size=legend.text.size),
                legend.text=element_text(size=legend.text.size)
                )
                
        # Save to new slot
        MarvelObject$Modality$Prop$DoughnutChart$Plot <- plot
    
    return(MarvelObject)
            
}
