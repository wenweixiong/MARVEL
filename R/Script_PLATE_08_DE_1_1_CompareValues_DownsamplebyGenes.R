#' @title Downsample cells based on number of expressed genes
#'
#' @description Match the cells from both cell groups based on the number of expressed genes.
#'
#' @param MarvelObject Marvel object. S3 object generated from \code{TransformExpValues} function.
#' @param cell.group.g1 Vector of character strings. Cell IDs corresponding to Group 1 (reference group).
#' @param cell.group.g2 Vector of character strings. Cell IDs corresponding to Group 2.
#' @param sliding.window Numeric values. The matching cell to look for with the number of genes detected within plus minus of this specified value.
#'
#' @return An object of class S3 new slots \code{MarvelObject$DE$Exp$Downsampled.Data$Table}, \code{MarvelObject$DE$Exp$Downsampled.Data$Plot.Before.Downsample}, and \code{MarvelObject$DE$Exp$Downsampled.Data$Plot.After.Downsample}.
#'
#' @importFrom plyr join
#' @import stats
#' @import methods
#' @import utils
#' @export

DownsampleByGenes <- function(MarvelObject, cell.group.g1, cell.group.g2, sliding.window=50) {

    # Define arguments
    df <- MarvelObject$Exp
    df.pheno <- MarvelObject$SplicePheno
    df.feature <- MarvelObject$GeneFeature
    cell.group.g1 <- cell.group.g1
    cell.group.g2 <- cell.group.g2
    sliding.window <- sliding.window
    
    # Define arguments
    #df <- marvel$Exp
    #df.pheno <- marvel$SplicePheno
    #df.feature <- marvel$GeneFeature
    #cell.group.g1 <- cell.group.g1
    #cell.group.g2 <- cell.group.g2
    #sliding.window <- 50
    
    ######################################################
    
    # Create row names for matrix
    row.names(df) <- df$gene_id
    df$gene_id <- NULL
    
    # Retrieve sample ids
        # Group 1
        sample.ids.1 <- cell.group.g1
        
        # Group 2
        sample.ids.2 <- cell.group.g2
    
    # Compute n genes detected
        # Group 1
        df.small <- df[, sample.ids.1]
        . <- apply(df.small, 2, function(x) {sum(x >=1)})
        results.1 <- data.frame("cell.group"="cell.group.g1",
                                "sample.id"=sample.ids.1,
                                "n.genes.detected"=as.numeric(.),
                                stringsAsFactors=FALSE
                                )
                                
        # Group 1
        df.small <- df[, sample.ids.2]
        . <- apply(df.small, 2, function(x) {sum(x >=1)})
        results.2 <- data.frame("cell.group"="cell.group.g2",
                                "sample.id"=sample.ids.2,
                                "n.genes.detected"=as.numeric(.),
                                stringsAsFactors=FALSE
                                )

        # Merge
        results <- rbind.data.frame(results.1, results.2)

    # Boxplot
        # Definition
        data <- results
        x <- data$cell.group
        y <- data$n.genes.detected
        z <- data$cell.group
        maintitle <- ""
        ytitle <- "Genes detected (n)"
        xtitle <- ""
        
        # Create x-axis labels
        . <- as.data.frame(table(x))
        xlabels.1 <- paste(.[,1], "\n(n=", .[,2], ")", sep="")
        
        # Plot
        plot <- ggplot() +
            geom_boxplot(data, mapping=aes(x=x, y=y, fill=z), outlier.color="black") +
            #geom_jitter(data, mapping=aes(x=x, y=y), position=position_jitter(width=0.1, height=0), size=0.1, alpha=0.2) +
            #stat_summary(data, mapping=aes(x=x, y=y), geom="point", fun="mean", fill="red", col="black", size=2, shape=23) +
            scale_x_discrete(labels=xlabels.1) +
            labs(title=maintitle, x=xtitle, y=ytitle) +
            theme(panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.background = element_blank(),
                panel.border=element_blank(),
                plot.title=element_text(hjust = 0.5, size=12),
                plot.subtitle=element_text(hjust = 0.5, size=12),
                axis.line.y.left = element_line(color="black"),
                axis.line.x = element_line(color="black"),
                axis.title=element_text(size=12),
                axis.text.x=element_text(size=6, colour="black"),
                axis.text.y=element_text(size=10, colour="black"),
                legend.position="none"
                )
                
    ################################################################
    ########################## DOWNSAMPLE ##########################
    ################################################################
    
    # Find match based on smaller cell group
    if(nrow(results.1) < nrow(results.2)){
        
        sample.ids.2 <- list()
        
        match.ratio <- round(nrow(results.2)/nrow(results.1), digits=0)

        for(i in 1:nrow(results.1)) {
            
            query <- results.1[i, "n.genes.detected"]
            query.lower <- query - sliding.window
            query.upper <- query + sliding.window
            
            index <- which(results.2$n.genes.detected >= query.lower & results.2$n.genes.detected <= query.upper)
            sample.ids.match <- results.2$sample.id[index]
            
            if(i == 1) {
                
                sample.ids.2[[1]] <- sample.ids.match[c(1:match.ratio)]
                
            } else {
                
                sample.ids.match <- setdiff(sample.ids.match, unlist(sample.ids.2))
                sample.ids.2[[i]] <- sample.ids.match[c(1:match.ratio)]
                
            }
            
        }
        
        index <- sapply(sample.ids.2, function(x) {sum(is.na(x))==match.ratio})
        results.1.small <- results.1[!index, ]
        
        index <- which(results.2$sample.id %in% unlist(sample.ids.2))
        results.2.small <- results.2[index, ]
        
    } else if(nrow(results.2) < nrow(results.1)){
        
        sample.ids.1 <- list()

        match.ratio <- round(nrow(results.1)/nrow(results.2), digits=0)

        for(i in 1:nrow(results.2)) {
            
            query <- results.2[i, "n.genes.detected"]
            query.lower <- query - sliding.window
            query.upper <- query + sliding.window
            
            index <- which(results.1$n.genes.detected >= query.lower & results.1$n.genes.detected <= query.upper)
            sample.ids.match <- results.1$sample.id[index]
            
            if(i == 1) {
                
                sample.ids.1[[1]] <- sample.ids.match[c(1:match.ratio)]
                
            } else {
                
                sample.ids.match <- setdiff(sample.ids.match, unlist(sample.ids.1))
                sample.ids.1[[i]] <- sample.ids.match[c(1:match.ratio)]
                
            }
            
        }

        index <- sapply(sample.ids.1, function(x) {sum(is.na(x))==match.ratio})
        results.2.small <- results.2[!index, ]

        index <- which(results.1$sample.id %in% unlist(sample.ids.1))
        results.1.small <- results.1[index, ]
        
    }
    
    results.small <- rbind.data.frame(results.1.small, results.2.small)
    
    # Track progress
    print(paste(nrow(results.1.small), " cells downsampled from ", nrow(results.1), " cells for Group1 1", sep=""))
    print(paste(nrow(results.2.small), " cells downsampled from ", nrow(results.2), " cells for Group1 1", sep=""))
    
    # Boxplot
        # Definition
        data.2 <- results.small
        x.2 <- data.2$cell.group
        y.2 <- data.2$n.genes.detected
        z.2 <- data.2$cell.group
        maintitle <- ""
        ytitle <- "Genes detected (n)"
        xtitle <- ""
        
        # Create x-axis labels
        . <- as.data.frame(table(x.2))
        xlabels.2 <- paste(.[,1], "\n(n=", .[,2], ")", sep="")
        
        # Plot
        plot.2 <- ggplot() +
            geom_boxplot(data.2, mapping=aes(x=x.2, y=y.2, fill=z.2), outlier.color="black") +
            #geom_jitter(data, mapping=aes(x=x, y=y), position=position_jitter(width=0.1, height=0), size=0.1, alpha=0.2) +
            #stat_summary(data, mapping=aes(x=x, y=y), geom="point", fun="mean", fill="red", col="black", size=2, shape=23) +
            scale_x_discrete(labels=xlabels.2) +
            labs(title=maintitle, x=xtitle, y=ytitle) +
            theme(panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.background = element_blank(),
                panel.border=element_blank(),
                plot.title=element_text(hjust = 0.5, size=12),
                plot.subtitle=element_text(hjust = 0.5, size=12),
                axis.line.y.left = element_line(color="black"),
                axis.line.x = element_line(color="black"),
                axis.title=element_text(size=12),
                axis.text.x=element_text(size=6, colour="black"),
                axis.text.y=element_text(size=10, colour="black"),
                legend.position="none"
                )
                
    ######################################################
  
    # Save into new slot
    MarvelObject$DE$Exp$Downsampled.Data$Table <- results.small
    MarvelObject$DE$Exp$Downsampled.Data$Plot.Before.Downsample <- plot
    MarvelObject$DE$Exp$Downsampled.Data$Plot.After.Downsample <- plot.2
    
    # Return final object
    return(MarvelObject)
        
}
