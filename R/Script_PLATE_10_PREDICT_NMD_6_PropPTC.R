#' @title Tabulate proportion of transcripts with PTC
#'
#' @description Tabulates and plots the proportion of transcripts with PTC for each splicing event type.
#'
#' @param MarvelObject Marvel object. S3 object generated from \code{FindPTC} function.
#' @param xlabels.size Numeric value. Size of the x-axis tick labels. Default is 8.
#' @param show.NovelSJ.NoCDS Logical value. If set to \code{TRUE} transcripts not analysed for premature terminal codon (PTC), e.g. non-protein-coding transcripts are tabulated and plotted.
#' @param prop.test Character string. \code{chisq} Chi-squared test used to compare the proportion of transcripts with PTC across the different event splicing type. \code{fisher} Fisher test used to compare the proportion of transcripts with PTC across the different splicing event type.
#'
#' @return An object of class S3 with new slots \code{MarvelObject$NMD$PTC.Prop$Table}, \code{MarvelObject$NMD$PTC.Prop$Plot}, and \code{MarvelObject$NMD$PTC.Prop$Plot.Stats}.
#'
#' @importFrom plyr join
#' @importFrom stats aggregate chisq.test fisher.test
#' @import ggplot2
#' @import scales
#'
#' @export
#'
#' @examples
#' marvel.demo <- readRDS(system.file("extdata/data", "marvel.demo.rds", package="MARVEL"))
#'
#' marvel.demo <- PropPTC(MarvelObject=marvel.demo,
#'                        xlabels.size=8,
#'                        show.NovelSJ.NoCDS=TRUE,
#'                        prop.test="fisher"
#'                        )
#'
#' # Check outputs
#' head(marvel.demo$NMD$PTC.Prop$Table)
#' marvel.demo$NMD$PTC.Prop$Plot
#' marvel.demo$NMD$PTC.Prop$Plot.Stats

PropPTC <- function(MarvelObject, xlabels.size=8, show.NovelSJ.NoCDS=TRUE, prop.test) {

    # Define arguments
    df <- MarvelObject$NMD$Prediction
    xlabels.size <- xlabels.size
    show.NovelSJ.NoCDS <- show.NovelSJ.NoCDS
    prop.test <- prop.test
    
    # Example arguments
    #MarvelObject <- marvel
    #df <- marvel$NMD$Prediction
    #xlabels.size <- 8
    #show.NovelSJ.NoCDS <- TRUE
    #prop.test <- "fisher"
    
    # Set factor levels
    levels <- intersect(c("SE", "RI", "A5SS", "A3SS"), unique(df$event_type))
    df$event_type <- factor(df$event_type, levels=levels)
    
    # Recode PTC
        # PTC
        . <- "yes_distance of PTC to final SJ > 50bp"
        df$NMD[which(df$NMD %in% "yes_distance of PTC to final SJ > 50bp")] <- ">=1 transcript with PTC"
        
        # No PTC
        . <- c("no_distance of PTC to final SJ <= 50bp",
               "no_PTC not found",
               "splicing event located outside of ORF"
               )
        df$NMD[which(df$NMD %in% .)] <- "No transcripts with PTC"
        
        # Non-protein coding
        . <- c("non-protein-coding transcript",
                "transcript with no START and/or STOP codon"
                )
        df$NMD[which(df$NMD %in% .)] <- "Non-coding transcripts only"
        
        # Novel SJ
        . <- c("No transcripts with matching SJ")
        df$NMD[which(df$NMD %in% .)] <- "Novel SJ"
        
        # Set factor levels
        if(show.NovelSJ.NoCDS==TRUE) {
            
            df$NMD <- factor(df$NMD, levels=c("Novel SJ", "Non-coding transcripts only", "No transcripts with PTC", ">=1 transcript with PTC"), labels=c("Novel SJ", "No CDS", "No PTC", "PTC"))
            
        } else if(show.NovelSJ.NoCDS==FALSE) {
            
            df <- df[which(df$NMD != "Non-coding transcripts only"), ]
            df$NMD <- factor(df$NMD, levels=c("No transcripts with PTC", ">=1 transcript with PTC"), labels=c("No PTC", "PTC"))
            
        }

    # X-tab
    xtab <- as.data.frame(table(df$event_type, df$NMD))
    names(xtab) <- c("event_type", "NMD", "Freq")
    
    # Create x-labels
    n <- aggregate(Freq ~ event_type, data=xtab, sum)
    xlabels <- paste(n$event_type, "\n(n=", n$Freq, ")", sep="")
    
    # Barplot
        # Definition
        #data <- df[which(df$NMD != "Non-CD"), ]
        data <- xtab
        x <- data$event_type
        y <- data$Freq
        z <- data$NMD
        maintitle <- ""
        xtitle <- ""
        ytitle <- "% Isoform"
        xlabels <- xlabels
        legendtitle <- "Isoform\nType"
        
        # Plot
        plot <- ggplot() +
            geom_bar(data, mapping=aes(x=x, y=y, fill=z), position="fill", stat="identity", color="black") +
            scale_x_discrete(labels=xlabels) +
            scale_y_continuous(labels=percent_format()) +
            #scale_fill_manual(values=c("lightblue", "red", "green"), name=legend.title) +
            labs(title=NULL, x=NULL, y=ytitle, fill=legendtitle) +
            theme(panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.background = element_blank(),
                panel.border=element_blank(),
                plot.title=element_text(hjust = 0.5, size=15),
                plot.subtitle=element_text(hjust = 0.5, size=15),
                axis.line.y.left = element_line(color="black"),
                axis.line.x = element_line(color="black"),
                axis.title=element_text(size=12),
                axis.text.x=element_text(size=xlabels.size, colour="black"),
                axis.text.y=element_text(size=10, colour="black"),
                legend.title=element_text(size=8),
                legend.text=element_text(size=8)
                )
                
    # Categorical test
        # X-tabulate
        xtab <- table(df$event_type, df$NMD)
        
        # X-square
        if(prop.test=="chisq") {
            
            pval <- chisq.test(xtab)$p.value
            
        } else if(prop.test=="fisher"){
            
            pval <- fisher.test(xtab)$p.value
            
        }
        
        # Convert n to %
        xtab <- as.data.frame.matrix(xtab)
        xtab <- apply(xtab, 1, function(x){
            
            paste(x, " (", round(x/sum(x)*100, digits=1), "%)", sep="")
            
        })
        xtab <- as.data.frame(xtab)
        xtab <- as.data.frame(t(xtab))
        names(xtab) <- levels(df$NMD)
        . <- data.frame("event_type"=row.names(xtab), stringsAsFactors=FALSE)
        xtab <- cbind.data.frame(., xtab)
        xtab$pval <- c(pval, rep("", times=nrow(xtab)-1))
        
    # Save to new slot
    MarvelObject$NMD$PTC.Prop$Table <- df
    MarvelObject$NMD$PTC.Prop$Plot <- plot
    MarvelObject$NMD$PTC.Prop$Plot.Stats <- xtab
 
    # Return MARVEL object
    return(MarvelObject)
        
}
