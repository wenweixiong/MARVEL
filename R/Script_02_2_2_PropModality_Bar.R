#' @title Modality Proportion across Event Types.
#'
#' @description
#' \code{PropModality.Bar} compares proportion of each modality across specified splicing event types.
#'
#' @details
#' This function compares proportion of each modality across specified splicing event types.
#'
#' @param MarvelObject S3 object generated from \code{AssignModality} function.
#' @param modality.column Character string. Can take the value \code{"modality"}, \code{"modality.var"} or \code{"modality.bimodal.adj"}. Please refer to \code{AssignModality} function help page for more details.
#' @param modality.type Character string. \code{basic} indicates that only the main modalities (included, excluded, bimodal, middle, multimodal) are analysed. Sub-modalities (primary and dispersed) will be merged. \code{extended} indicates that both main and sub-modalities are analysed. Sub-modalities will not be merged.
#' @param event.type Character string. To indicate which event type to analyse. Can take the value \code{"SE"}, \code{"MXE"}, \code{"RI"}, \code{"A5SS"} or \code{"A3SS"}. Specify \code{"all"} to include all event types.
#' @param prop.test Character string. Only applicable when \code{across.event.type} set to \code{TRUE}. \code{chisq} Chi-squared test used to compare the proportion of modalities across the different event splicing type. \code{fisher} Fisher test used to compare the proportion of modalities across the different splicing event type.
#' @param prop.adj Character string. Only applicable when \code{across.event.type} set to \code{TRUE}. Adjust p-values generated from \code{prop.test} for multiple testing. Options available as per \code{p.adjust} function.
#' @export
#' @return An object of class S3 containing all the original slots as inputted by the user in addition to one new slot named \code{$Modality$Prop$BarChart} slot.
#' @author Sean Wen <sean.wenwx@gmail.com>
#' @importFrom fitdistrplus fitdist
#' @import methods
#' @import ggplot2
#' @importFrom plyr join
#' @examples
#' marvel <- readRDS(system.file("extdata/Data", "MarvelObject.rds", package="MARVEL"))
#'
#' marvel <- PropModality.Bar(MarvelObject=marvel,
#'                            modality.column="modality.bimodal.adj",
#'                            modality.type="basic",
#'                            event.type="all",
#'                            prop.test="fisher",
#'                            prop.adj="fdr"
#'                           )
#'
#' marvel$Modality$Prop$BarChart$Plot

PropModality.Bar <- function(MarvelObject, modality.column, modality.type, event.type, prop.test, prop.adj) {
    
    # Define arguments
    df.feature <- MarvelObject$Modality$Results
    modality.column <- modality.column
    modality.type <- modality.type
    event.type <- event.type
    prop.test <- prop.test
    prop.adj <- prop.adj

    if(event.type[1]=="all") {
        
        event.type <- c("SE", "MXE", "RI", "A3SS", "A5SS")
        
    }
    
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
    
    # Tabulate modality by event type
    event_types <- event.type

    .list <- list()

    for(i in 1:length(event_types)) {

        # Subset relevant event type
        . <- df.feature[which(df.feature$event_type==event_types[i]), ]

        # Compute %
        . <- as.data.frame(table(.$modality), stringsAsFactors=FALSE)
        names(.) <- c("modality", "freq")
        .$pct <- .$freq / sum(.$freq) * 100

        # Set factor levels
        if(modality.type=="basic") {

            .$modality <- factor(.$modality, levels=c("Included", "Excluded", "Bimodal", "Middle", "Multimodal"))
            . <- .[order(.$modality), ]
            
        } else {
            
            .$modality <- factor(.$modality, levels=c("Included.Primary", "Included.Dispersed", "Excluded.Primary", "Excluded.Dispersed", "Bimodal", "Middle", "Multimodal"))
            . <- .[order(.$modality), ]
            
        }
        
        # Indicate event type
        .$event_type <- event_types[i]
        
        # Save into list
        .list[[i]] <- .
        
    }
    
    . <- do.call(rbind.data.frame, .list)
    
    # Set factor levels
    .$event_type <- factor(.$event_type, levels=event.type)

    # Save to new slot
    MarvelObject$Modality$Prop$BarChart$Table <- .
    
    # Plot
        # Definitions
        data <- .
        x <- data$modality
        y <- data$pct
        z <- data$event_type
        maintitle <- ""
        ytitle <- "%"
        xtitle <- ""
        #fivenum(y) ; ymin <- 0 ; ymax <- 76 ; yinterval <- 10
        legendtitle <- "Event Type"
        
        if(modality.type=="basic") {
        
            axis.text.x.size <- 10
            
        } else {
            
            axis.text.x.size <- 7
        }
        
        # Plot
        plot <- ggplot() +
            geom_bar(data=data, aes(x=x, y=y, fill=z), stat="identity", color="black", position=position_dodge(), width=0.8) +
            #scale_y_continuous(breaks=seq(ymin, ymax, by=yinterval), limits=c(ymin, ymax)) +
            labs(title=maintitle, x=xtitle, y=ytitle, fill=legendtitle) +
            theme(panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.background = element_blank(),
                panel.border=element_blank(),
                plot.title=element_text(hjust = 0.5, size=15),
                plot.subtitle=element_text(hjust = 0.5, size=15),
                axis.line.y.left = element_line(color="black"),
                axis.line.x = element_line(color="black"),
                axis.title=element_text(size=15),
                axis.text=element_text(size=10),
                axis.text.x=element_text(size=axis.text.x.size, colour="black"),
                axis.text.y=element_text(size=10, colour="black")
                )

        # Save to new slot
        MarvelObject$Modality$Prop$BarChart$Plot <- plot
        
    # x^2 test
        # Compute the sum by event type
        total <- tapply(.$freq, .$event_type, sum)
        total <- data.frame(total)
        total$event_type <- row.names(total)
        
        # Tabulate freq (others)
        mod <- levels(.$modality)
        
        p.val <- NULL
        
        for(i in 1:length(mod)) {
        
            .small <- .[which(.$modality==mod[i]),]
            
            if(nrow(.small) >= 2) {
                
                
                .small <- join(.small, total, by="event_type", type="left")
                .small$freq.others <- .small$total - .small$freq
                
                if(prop.test=="chisq") {
                    
                    p.val[i] <- chisq.test(.small[,c("freq", "freq.others")])$p.value
                    
                } else {
                    
                    p.val[i] <- fisher.test(.small[,c("freq", "freq.others")])$p.value
                    
                }
            
            
            } else {
                
                p.val[i] <- NA
                
            }
                        
        }

        results <- data.frame("modality"=mod, "p.val"=p.val, stringsAsFactors=FALSE)
        
        # Adjuste for multiple testing
        results$p.val.adj <- p.adjust(results$p.val, method=prop.adj, n=length(results$p.val))

        # Save to new slot
        MarvelObject$Modality$Prop$BarChart$Test <- results
        
        
    return(MarvelObject)
            
}
