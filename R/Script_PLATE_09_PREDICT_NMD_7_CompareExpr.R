#' @title Compares gene expression changes based on nonsense-mediated decay (NMD) status
#'
#' @description Compares gene expression changes based on NMD status for each splicing event type.
#'
#' @param MarvelObject Marvel object. S3 object generated from \code{FindPTC} function.
#' @param xlabels.size Numeric value. Size of the x-axis tick labels. Default is 8.
#'
#' @return An object of class S3 new slots \code{MarvelObject$NMD$NMD.Expr$Table}, \code{MarvelObject$NMD$NMD.Expr$Plot}, and \code{MarvelObject$NMD$NMD.Expr$Plot.Stats}.
#'
#' @importFrom plyr join
#' @import ggplot2
#'
#' @export
#'
#' @examples
#' marvel.demo <- readRDS(system.file("extdata/data", "marvel.demo.rds", package="MARVEL"))
#'
#' marvel.demo <- CompareExpr(MarvelObject=marvel.demo)
#'
#' # Check outputs
#' head(marvel.demo$NMD$NMD.Expr$Table)
#' marvel.demo$NMD$NMD.Expr$Plot
#' marvel.demo$NMD$NMD.Expr$Plot.Stats

CompareExpr <- function(MarvelObject, xlabels.size=8) {

    # Define arguments
    df <- MarvelObject$NMD$Prediction
    de.gene <- MarvelObject$DE$Exp.Spliced$Table
    xlabels.size <- xlabels.size
    
    # Example arguments
    #MarvelObject <- marvel
    #df <- MarvelObject$NMD$Prediction
    #de.gene <- MarvelObject$DE$Exp.Spliced$Table
    #xlabels.size <- 8

    # Set factor levels
    levels <- intersect(c("SE", "RI", "A5SS", "A3SS"), unique(df$event_type))
    df$event_type <- factor(df$event_type, levels=levels)
        
    # Collapse by gene for each event type
    .list <- list()
    
    for(j in 1:length(levels)) {
        
        # Subset relevant event type
        df.small <- df[which(df$event_type==levels[j]), ]
        
        # Collapse by gene
        gene_ids <- unique(df.small$gene_id)
        
        nmd <- NULL
        
        for(i in 1:length(gene_ids)) {
        
            # Subset relevant gene
            df.small. <- df.small[which(df.small$gene_id==gene_ids[i]), ]
            
            # Check for NMD
            index.nmd.yes <- grep("^yes_", df.small.$NMD)

            index.nmd.no <- unique(c(
                                    grep("no_distance of PTC to final SJ <= 50bp", df.small.$NMD, fixed=TRUE),
                                    grep("no_PTC not found", df.small.$NMD, fixed=TRUE),
                                    grep("splicing event located outside of ORF", df.small.$NMD, fixed=TRUE),
                                    grep("No transcripts with matching SJ", df.small.$NMD, fixed=TRUE)
                                    ))
                                    
            index.npc <- unique(c(
                             grep("non-protein-coding transcript", df.small.$NMD, fixed=TRUE),
                             grep("transcript with no START and/or STOP codon", df.small.$NMD, fixed=TRUE)
                         ))
                            
            
            if(length(index.nmd.yes) >= 1) {
                
                nmd[i] <- ">=1 transcript with PTC"
                
            } else if(length(index.nmd.yes)==0 & length(index.nmd.no) >= 1) {
                
                
                nmd[i] <- "No transcripts with PTC"
                
            } else if(length(index.npc) >= 1) {
                
                nmd[i] <- "Non-coding transcripts only"
            
            }
        
        }
    
        # Tabulate results
        .list[[j]] <- data.frame("gene_id"=gene_ids, "event_type"=levels[j], "NMD"=nmd, stringsAsFactors=FALSE)
    
    }
    
    # Tabulate results
    results <- do.call(rbind.data.frame, .list)
    
    # Set factor levels
    results$NMD <- factor(results$NMD, levels=c("Non-coding transcripts only", "No transcripts with PTC", ">=1 transcript with PTC"), labels=c("FALSE", "FALSE", "TRUE"))
    results$event_type <- factor(results$event_type, levels=levels)
    
    # Annotate gene log2fc
    results <- join(results, de.gene[,c("gene_id", "gene_short_name", "log2fc")], by="gene_id", type="left")
    cols <- c("gene_id", "gene_short_name", "log2fc", "event_type", "NMD")
    results <- results[,cols]
    
    # Add sample size to xlabels
        # non-NMD
        . <- results[which(results$NMD=="FALSE"), ]
        . <- as.data.frame(table(.$event_type))
        names(.) <- c("event_type", "n.nmd.false")
        nmd.false <- .
        
        # NMD
        . <- results[which(results$NMD=="TRUE"), ]
        . <- as.data.frame(table(.$event_type))
        names(.) <- c("event_type", "n.nmd.true")
        nmd.true <- .
        
        # Merge
        n <- join(nmd.false, nmd.true, by="event_type", type="full")
        n[is.na(n)] <- 0
        xlabels <- paste(n$event_type, "\n(n=", n$n.nmd.false, ",", n$n.nmd.true, ")", sep="")
    
    # Boxplot
        # Definition
        data <- results
        x <- data$event_type
        y <- data$log2fc
        z <- data$NMD
        maintitle <- ""
        ytitle <- "log2FC (Gene Expression)"
        xtitle <- ""
        legendtitle <- "NMD"
        #xlabels <- n.cells$label
        #fivenum(y) ; ymin <- -0.5 ; ymax <- 1.0 ; yinterval <- 0.25

        # Plot
        plot <- ggplot() +
            geom_boxplot(data, mapping=aes(x=x, y=y, fill=z), outlier.size=0.1) +
            #geom_jitter(data, mapping=aes(x=x, y=y), position=position_jitter(width=0.1, height=0), size=0.001) +
            #stat_summary(data, mapping=aes(x=x, y=y), geom="point", fun="mean", fill="red", col="black", size=2, shape=23) +
            scale_x_discrete(labels=xlabels) +
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
                axis.title=element_text(size=12),
                axis.text.x=element_text(size=xlabels.size, colour="black"),
                axis.text.y=element_text(size=10, colour="black"),
                legend.title=element_text(size=8),
                legend.text=element_text(size=8)
                )
        
    # Pair-wise comparison
    nmd.null.cells <- NULL
    nmd.cells <- NULL
    nmd.null.log2fc.mean <- NULL
    nmd.log2fc.mean <- NULL
    p.val <- NULL
    
    for(i in 1:length(levels)) {
        
        # Subset relevant event type
        results.small <- results[which(results$event_type==levels[i]), ]
        
        # Retrieve values
        nmd.null <- results.small[which(results.small$NMD=="FALSE"), "log2fc"]
        nmd <- results.small[which(results.small$NMD=="TRUE"), "log2fc"]
        
        # Tabulate sample size
        nmd.null.cells[i] <- length(nmd.null)
        nmd.cells[i] <- length(nmd)
                    
        # Tabulate mean
        nmd.null.log2fc.mean[i] <- mean(nmd.null)
        nmd.log2fc.mean[i] <- mean(nmd)
        
        if(length(nmd.null) != 0 & length(nmd) != 0) {
            
            # Wilcox
            p.val[i] <- wilcox.test(nmd.null, nmd)$p.value
        
        } else {
            
            p.val[i] <- NA
            
        }
        
    }
        
    stats <- data.frame(event_type=levels,
                        "nmd.null.cells"=nmd.null.cells,
                        "nmd.cells"=nmd.cells,
                        "nmd.null.log2fc.mean"=nmd.null.log2fc.mean,
                        "nmd.log2fc.mean"=nmd.log2fc.mean,
                        "p.val"=p.val,
                        stringsAsFactors=FALSE
                        )
                        
    # Save to new slot
    MarvelObject$NMD$NMD.Expr$Table <- results
    MarvelObject$NMD$NMD.Expr$Plot <- plot
    MarvelObject$NMD$NMD.Expr$Plot.Stats <- stats
 
    # Return MARVEL object
    return(MarvelObject)
        
}
