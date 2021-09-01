#' @title Plot Differential Splicing Analysis Results
#'
#' @description
#' \code{PlotDEValues.PSI.Distance} plots the differential splicing analysis results.
#'
#' @details
#' This function plot the differential splicing analysis results and can annotate the user-defined significant data points and gene names. Typically used when the distribution of percent spliced-in (PSI) values of 2 groups of cells were previously tested using Kolmogorov-Smirnov, Kuiper or Anderson-Darling test.
#'
#' @param MarvelObject S3 object generated from \code{CompareValues} function.
#' @param de.p.val.adj Numeric value. Adjusted p-value below which the splcing events are considered as statistically significant and will consequently be color-annotated on the plot.
#' @param anno Logical value. If set to \code{TRUE}, the specific gene names will be annotated on the plot. Speficified together with \code{n.top}.
#' @param n.top Numeric value. Only applicable when \code{anno} set to \code{TRUE}. Top n differentially expressed splicing events are annotated on the plot.
#' @param label.size Numeric value. Only applicable if \code{anno} set to TRUE. Size of the gene name labels.
#'
#' @export
#'
#' @return An object of class S3 containing all the original slots as inputted by the user in addition to one new slot named \code{MarvelObject$DE$PSI}.
#'
#' @author Sean Wen <sean.wenwx@gmail.com>
#'
#' @importFrom plyr join
#' @import stats
#' @import ggplot2
#' @import ggrepel
#' @import scales
#'
#' @examples
#' # Load input
#' marvel <- readRDS(system.file("extdata/Data", "MarvelObject.rds", package="MARVEL"))
#'
#' # Run example
#' marvel <- PlotDEValues.PSI.Distance(MarvelObject=marvel,
#'                                     de.p.val.adj=0.10,
#'                                     anno=FALSE
#'                                     )
#'
#' # Check output
#' marvel$DE$PSI$Plot


PlotDEValues.PSI.Distance <- function(MarvelObject, de.p.val.adj, anno=FALSE, n.top=NULL, label.size=2.5) {

    # Define arguments
    df <- MarvelObject$DE$PSI$Table
    de.p.val.adj <- de.p.val.adj
    anno <- anno
    n.top <- n.top
    label.size <- label.size
    
    # Example arguments
    #df <- marvel$DE$PSI
    #de.p.val.adj <- 0.10
    #anno <- TRUE
    #n.top <- 5
    #label.size <- 2.5
    
    # Indicate sig events and direction
    df$sig <- ifelse(df$p.val.adj < de.p.val.adj, "sig", "n.s.")
    df$sig <- factor(df$sig, levels=c("sig", "n.s."))
    
    # Indicate direction on D statistic
   # df$statistic <- ifelse(df$mean.diff < 0, df$statistic * -1, df$statistic)
    
    # Create labels
    if(anno==TRUE) {
        
        #df$label <- ifelse((df$statistic > anno.x.above &
                           # -log10(df$p.val.adj) > anno.y.above),
                           # paste(df$gene_short_name, ".", df$event_type, sep=""),
                           # NA
                           # )
                           
        df$label <- paste(df$gene_short_name, ".", df$event_type, sep="")
        df$label[-c(1:n.top)] <- NA
                           
       # Definition
       data <- df
       x <- data$statistic + 1
       y <- -log10(data$p.val.adj)
       z <- data$sig
       label <- data$label
       maintitle <- ""
       xtitle <- "D Statistic + 1"
       ytitle <- "-log10(p-value)"

       ymin <- 0 ; ymax <- max(y) + 5 ; yinterval <- 5
           
       sig <- sum(data$sig=="sig")
       n.s. <- sum(data$sig=="n.s.")

       if(sig != 0 & n.s. != 0) {

           col.breaks <- c("red", "gray")
           
       } else if(sig != 0 & n.s. == 0) {

           col.breaks <- "red"
           
       } else if(sig == 0 & n.s. != 0) {

           col.breaks <- "gray"

       }

       # Plot
       plot <- ggplot() +
                  geom_point(data, mapping=aes(x=x, y=y, color=z), shape=20, alpha = 0.75) +
                  geom_text_repel(data, mapping=aes(x=x, y=y, label=label), max.overlaps = Inf, box.padding = 1.0, size=label.size, max.time = 1, max.iter = 1e5, segment.alpha=0.5, segment.size=0.1, min.segment.length = 0) +
                  scale_colour_manual(values=col.breaks) +
                  scale_x_log10(breaks=c(1, 10, 25, 50, 100), limits=c(1, 200)) +
                  scale_y_continuous(breaks=seq(ymin, ymax, by=yinterval), limits=c(ymin, ymax)) +
                  labs(title=maintitle, x=xtitle, y=ytitle) +
                   theme(panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(),
                       panel.background = element_blank(),
                       panel.border=element_blank(),
                       plot.title=element_text(hjust = 0.5, size=15),
                       plot.subtitle=element_text(hjust = 0.5, size=15),
                       axis.line.y.left = element_line(color="black"),
                       axis.line.x = element_line(color="black"),
                       axis.title=element_text(size=12),
                       axis.text=element_text(size=12),
                       axis.text.x=element_text(size=10, colour="black"),
                       axis.text.y=element_text(size=10, colour="black"),
                       legend.position="none",
                       legend.title=element_text(size=8),
                       legend.text=element_text(size=8)
                       )
    
    } else {
        
        # Definition
        data <- df
        x <- data$statistic + 1
        y <- -log10(data$p.val.adj)
        z <- data$sig
        maintitle <- ""
        xtitle <- "D Statistic + 1"
        ytitle <- "-log10(p-value)"

        ymin <- 0 ; ymax <- max(y) ; yinterval <- 5
            
        sig <- sum(data$sig=="sig")
        n.s. <- sum(data$sig=="n.s.")

        if(sig != 0 & n.s. != 0) {

            col.breaks <- c("red", "gray")
            
        } else if(sig != 0 & n.s. == 0) {

            col.breaks <- "red"
            
        } else if(sig == 0 & n.s. != 0) {

            col.breaks <- "gray"

        }

        # Plot
        plot <- ggplot() +
                   geom_point(data, mapping=aes(x=x, y=y, color=z), shape=20, alpha = 0.75) +
                   scale_colour_manual(values=col.breaks) +
                   scale_x_log10(breaks=c(1, 10, 25, 50, 100), limits=c(1, 200)) +
                   scale_y_continuous(breaks=seq(ymin, ymax, by=yinterval), limits=c(ymin, ymax)) +
                   labs(title=maintitle, x=xtitle, y=ytitle) +
                    theme(panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(),
                        panel.background = element_blank(),
                        panel.border=element_blank(),
                        plot.title=element_text(hjust = 0.5, size=15),
                        plot.subtitle=element_text(hjust = 0.5, size=15),
                        axis.line.y.left = element_line(color="black"),
                        axis.line.x = element_line(color="black"),
                        axis.title=element_text(size=12),
                        axis.text=element_text(size=12),
                        axis.text.x=element_text(size=10, colour="black"),
                        axis.text.y=element_text(size=10, colour="black"),
                        legend.position="none",
                        legend.title=element_text(size=8),
                        legend.text=element_text(size=8)
                        )
    
    }
                          
    # Save to new slot
    MarvelObject$DE$PSI$Plot <- plot
  
    return(MarvelObject)
        
}
