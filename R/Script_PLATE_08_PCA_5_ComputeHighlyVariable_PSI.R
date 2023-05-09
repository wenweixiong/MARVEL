#' @title Identify highly variable splicing events
#'
#' @description Identify highly variable splicing events by fitting spline model. The vector splicing events returned may be used for downstream dimension reduction analysis with the \code{RunPCA} function.
#'
#' @param MarvelObject Marvel object. S3 object generated from \code{ComputePSI} function.
#' @param sample.ids Character strings. Specific cells to plot.
#' @param cell.group.column Character string. The name of the sample metadata column in which the variables will be used included for analysis
#' @param cell.group.order Character string. The the variables under the sample metadata column specified in \code{cell.group.column} to be included for analysis
#' @param min.cells Numeric value. The minimum no. of cells expressing the splicing event, above which, the splicing event will be retained for analysis.
#'
#' @return An object of class S3 containing with new slots \code{ MarvelObject$VariableSplicing$tran_ids} and  \code{MarvelObject$VariableSplicing$plot}
#'
#' @importFrom plyr join
#' @importFrom stats rnorm runif sd
#' @import methods
#' @import ggplot2
#' @importFrom grDevices hcl
#'
#' @export
#'

IdentifyVariableEvents <- function(MarvelObject, sample.ids=NULL,
                                   cell.group.column, cell.group.order,
                                   min.cells=25
                                   ) {

    # Define arguments
    df <- do.call(rbind.data.frame, MarvelObject$PSI)
    df.pheno <- MarvelObject$SplicePheno
    df.feature <- do.call(rbind.data.frame, MarvelObject$SpliceFeatureValidated)
    sample.ids <- sample.ids
    cell.group.column <- cell.group.column
    cell.group.order <- cell.group.order
    min.cells <- min.cells
    
    # Example arguments
    #MarvelObject <- marvel
    #df <- do.call(rbind.data.frame, MarvelObject$PSI)
    #df.pheno <- MarvelObject$SplicePheno
    #df.feature <- do.call(rbind.data.frame, MarvelObject$SpliceFeatureValidated)
    #sample.ids <- NULL
    #cell.group.column <- "cell.type"
    #cell.group.order <- c("0-hrs", "24-hrs", "48-hrs", "72-hrs")
    #min.cells <- 25
    
    ######################################################################
        
    # Create row names for matrix
    row.names(df) <- df$tran_id
    df$tran_id <- NULL
    
    # Rename cell group label
    df.pheno$pca.cell.group.label <- df.pheno[[cell.group.column]]
    
    # Subset relevant cells: overall
    if(!is.null(sample.ids[1])) {
        
        df.pheno <- df.pheno[which(df.pheno$sample.id %in% sample.ids), ]
        
    }
    
    # Subset relevant cells
        # Check if cell group order is defined
        if(is.null(cell.group.order[1])) {
            
            cell.group.order <- unique(df.pheno$pca.cell.group.label)
            
        }

        # Cell group
        index <- which(df.pheno$pca.cell.group.label %in% cell.group.order)
        df.pheno <- df.pheno[index, ]
        
        # Report progress
        print(paste(length(index), " of ", ncol(df), " cells specified were found and retained", sep=""))
        
        # Subset matrix
        df <- df[, df.pheno$sample.id]

    # Set factor levels
    levels <- intersect(cell.group.order, unique(df.pheno$pca.cell.group.label))
    df.pheno$pca.cell.group.label <- factor(df.pheno$pca.cell.group.label, levels=levels)

    # Subset events with sufficient cells
    . <- apply(df, 1, function(x) {sum(!is.na(x))})
    index.keep <- which(. >= min.cells)
    df <- df[index.keep, ]
    
    df.feature <- df.feature[which(df.feature$tran_id %in% row.names(df)), ]
    
    print(paste(length(index.keep), " of ", length(.), " splicing events expressed in at least ", min.cells, " cells and are retained", sep=""))
    
    # Compute features
        # Mean
        ave <- apply(df, 1, function(x) {mean(x, na.rm=TRUE)})
    
        # Standard deviation
        std <- apply(df, 1, function(x) {sd(x, na.rm=TRUE)})
        
        # Save into table
        results <- data.frame("tran_id"=row.names(df),
                              "mean"=ave,
                              "sd"=std,
                              stringsAsFactors=FALSE
                              )
        
        # Smoothed SD
        model <- mgcv::gam(sd ~ s(mean, bs="ps", sp=0.6), data=results)
        pred <- predict(model, results, type="link", se.fit=TRUE)
        
        results$sd_pred <- pred$fit
        results$sd_pred_ci_lower <- pred$fit - (2 * pred$se.fit)
        results$sd_pred_ci_upper <- pred$fit + (2 * pred$se.fit)
    
    # Indicate highly varible events
    results$variable <- ifelse(results$sd > results$sd_pred, TRUE, FALSE)
    
    print(paste(sum(results$variable==TRUE), " of ", nrow(results), " splicing events identified as highly variable", sep=""))
    
    # Scatter + lineplot
        # Definition
        data <- results
        x <- data$mean * 100
        y <- data$sd
        y2 <- data$sd_pred
        z <- data$variable
        maintitle <- paste(sum(results$variable==TRUE), " variable events", sep="")
        xtitle <- paste("Mean PSI")
        ytitle <- paste("SD PSI")
        legendtitle <- "Variable splicing events"
                
        # Plot
        plot <- ggplot() +
            geom_point(data, mapping=aes(x=x, y=y, color=z), size=0.1, alpha=0.5) +
            geom_line(data, mapping=aes(x=x, y=y2), color="blue") +
            scale_color_manual(values=c("black", "red")) +
            labs(title=maintitle, x=xtitle, y=ytitle, color=legendtitle) +
            theme(panel.grid.major=element_blank(),
                panel.grid.minor=element_blank(),
                panel.background=element_blank(),
                plot.title = element_text(size=12, hjust=0.5),
                axis.line=element_line(colour = "black"),
                axis.title=element_text(size=12),
                axis.text.x=element_text(size=10, colour="black"),
                axis.text.y=element_text(size=10, colour="black"),
                legend.title=element_text(size=8),
                legend.text=element_text(size=8)
                )  +
        guides(color = guide_legend(override.aes=list(size=2, alpha=1, stroke=0.1), ncol=1))
    
    # Retrieve variable splicing ids
    index <- which(results$variable==TRUE)
    tran_ids <- results[index, "tran_id"]
    
    ##############################################
     
    # Save to new slot
    MarvelObject$VariableSplicing$tran_ids <- tran_ids
    MarvelObject$VariableSplicing$plot <- plot
    
    return(MarvelObject)
        
}
