#' @title Modality proportion broken down by event type
#'
#' @description Tabulates and plots the proportion of each modality broken down by splicing event type.
#'
#' @param MarvelObject Marvel object. S3 object generated from \code{AssignModality} function.
#' @param modality.column Character string. Can take the value \code{"modality"}, \code{"modality.var"} or \code{"modality.bimodal.adj"}. Please refer to \code{AssignModality} function help page for more details.
#' @param modality.type Character string. \code{basic} indicates that only the main modalities (included, excluded, bimodal, middle, multimodal) are analysed. Sub-modalities (primary and dispersed) will be merged. \code{extended} indicates that both main and sub-modalities are analysed. Sub-modalities will not be merged.
#' @param event.type Character string. To indicate which event type to analyse. Can take the value \code{"SE"}, \code{"MXE"}, \code{"RI"}, \code{"A5SS"} or \code{"A3SS"}. Specify \code{"all"} to include all event types.
#' @param xlabels.size Numeric value. Size of x-axis labels as per \code{ggplot2} function. Default is 8.
#' @param zoom Logical value. If set to \code{TRUE}, users can specify the range of the y-axis using \code{yinterval} argument. Useful when scrutinasing low-frequency event types, e.g. middle and multimodal.
#' @param yinterval Logical value. Only applicable when \code{zoom} is set to \code{TRUE}.
#' @param prop.test Character string. Only applicable when \code{across.event.type} set to \code{TRUE}. \code{chisq} Chi-squared test used to compare the proportion of modalities across the different event splicing type. \code{fisher} Fisher test used to compare the proportion of modalities across the different splicing event type.
#' @param prop.adj Character string. Only applicable when \code{across.event.type} set to \code{TRUE}. Adjust p-values generated from \code{prop.test} for multiple testing. Options available as per \code{p.adjust} function.
#'
#' @return An object of class S3 containing new slots \code{MarvelObject$Modality$Prop$BarChart$Table} and \code{MarvelObject$Modality$Prop$BarChart$Stats}.
#'
#' @importFrom plyr join
#' @importFrom stats chisq.test fisher.test p.adjust p.adjust.methods
#' @import methods
#' @import ggplot2
#'
#' @export
#'
#' @examples
#' marvel.demo <- readRDS(system.file("extdata/data", "marvel.demo.rds", package="MARVEL"))
#'
#' marvel.demo <- PropModality.Bar(MarvelObject=marvel.demo,
#'                                 modality.column="modality.bimodal.adj",
#'                                 modality.type="extended",
#'                                 event.type=c("SE", "MXE", "RI", "A5SS", "A3SS", "AFE", "ALE"),
#'                                 prop.test="fisher",
#'                                 prop.adj="fdr"
#'                                 )
#'
#' # Check outputs
#' head(marvel.demo$Modality$Prop$BarChart$Table)
#' marvel.demo$Modality$Prop$BarChart$Plot
#' marvel.demo$Modality$Prop$BarChart$Stats

PropModality.Bar <- function(MarvelObject, modality.column, modality.type, event.type, xlabels.size=8, zoom=FALSE, yinterval=NULL, prop.test, prop.adj) {
    
    # Define arguments
    df.feature <- MarvelObject$Modality$Results
    modality.column <- modality.column
    modality.type <- modality.type
    event.type <- event.type
    xlabels.size <- xlabels.size
    prop.test <- prop.test
    prop.adj <- prop.adj
    zoom <- zoom
    yinterval <- yinterval
    
    # Example arguments
    #df.feature <- marvel$Modality$Results
    #modality.column <- "modality.bimodal.adj"
    #modality.type <- "extended"
    #event.type <- c("SE", "MXE", "RI", "A5SS", "A3SS")
    #xlabels.size <- 8
    #prop.test <- "chisq"
    #prop.adj <- "fdr"
    #zoom <- TRUE
    #yinterval <- c(0, 2.5)
        
    # Subset relevant modality column
    df.feature <- df.feature[,c("event_type", modality.column)]
    names(df.feature)[which(names(df.feature)==modality.column)] <- "modality"
    
    # Subset relevant event type
    event.type <- intersect(event.type, unique(df.feature$event_type))
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

            levels <- intersect(c("Included", "Excluded", "Bimodal", "Middle", "Multimodal"), unique(.$modality))
            .$modality <- factor(.$modality, levels=levels)
            . <- .[order(.$modality), ]
            
        } else if(modality.type=="extended"){
            
            levels <- intersect(c("Included.Primary", "Included.Dispersed", "Excluded.Primary", "Excluded.Dispersed", "Bimodal", "Middle", "Multimodal"), unique(.$modality))
            labels <- gsub(".", "\n", levels, fixed=TRUE)
            labels <- gsub("Primary", "(Primary)", labels, fixed=TRUE)
            labels <- gsub("Dispersed", "(Dispersed)", labels, fixed=TRUE)
            .$modality <- factor(.$modality, levels=levels, labels=labels)
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
    
    # Create dummy entry for missing modalities
    event_types <- levels(.$event_type)
    modalities <- levels(.$modality)
    
    for(i in 1:length(event_types)) {
        
        .small <- .[which(.$event_type==event_types[i]), ]
        
        if(length(modalities)==nrow(.small)) {
            
            .list[[i]] <- .small
            
        } else {
            
            modality.missing <- setdiff(modalities, .small$modality)
            modality.missing.df <- data.frame("modality"=modality.missing,
                                              "freq"=0,
                                              "pct"=0,
                                              "event_type"=event_types[i],
                                              stringsAsFactors=FALSE
                                              )
            .list[[i]] <- rbind.data.frame(.small, modality.missing.df)
            
        }
        
    }
    
    . <- do.call(rbind.data.frame, .list)
    . <- .[order(.$event_type, .$modality),]
    
    # Save to new slot
    temp <- .
    temp$modality <- gsub("\n", " ", temp$modality)
    
    MarvelObject$Modality$Prop$BarChart$Table <- temp
    
    # Plot
    if(zoom==FALSE) {
        
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
                
        # Plot
        plot <- ggplot() +
            geom_bar(data=data, aes(x=x, y=y, fill=z), stat="identity", color="black", position=position_dodge(), width=0.8) +
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
                axis.text.x=element_text(size=xlabels.size, colour="black"),
                axis.text.y=element_text(size=10, colour="black"),
                legend.title=element_text(size=8),
                legend.text=element_text(size=8)
                )

        # Save to new slot
        MarvelObject$Modality$Prop$BarChart$Plot <- plot
        
    } else {
        
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
                
        # Plot
        plot <- ggplot() +
            geom_bar(data=data, aes(x=x, y=y, fill=z), stat="identity", color="black", position=position_dodge(), width=0.8) +
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
                axis.text.x=element_text(size=xlabels.size, colour="black"),
                axis.text.y=element_text(size=10, colour="black"),
                legend.title=element_text(size=8),
                legend.text=element_text(size=8)
                ) +
                coord_cartesian(ylim=yinterval)
                
        # Save to new slot
        MarvelObject$Modality$Prop$BarChart$Plot <- plot
        
    }
        
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

    # Remove line breaks 
    results$modality <- gsub("\n", " ", results$modality, fixed=TRUE)
    
    # Save to new slot
    MarvelObject$Modality$Prop$BarChart$Stats <- results
        
    return(MarvelObject)
            
}
