#' @title Plot Percent Spliced-in (PSI) Values
#'
#' @description
#' \code{PlotValues.PSI} plots percent spliced-in (PSI) values across different groups of cells.
#'
#' @details
#' This function plots percent spliced-in (PSI) across different groups of cells. Violin plot is used for PSI values.
#'
#' @param MarvelObject S3 object generated from \code{CreateMarvelObject} or \code{ComputePSI} function.
#' @param cell.types Character string. To indicate which groups of cells that will be used for plotting. Group names should match those in \code{cell.type} column of \code{$SplicePheno} slot.
#' @param feature Character string. \code{tran_id} for plotting. Should match \code{tran_id} column of \code{$ValidatedSpliceFeature} slot.
#' @param maintitle Character string. Column to use as plot main title as per \code{ValidatedSpliceFeature}.
#' @param xlabels.size Numeric value. Size of x-axis labels as per \code{ggplot2} function.
#' @param n.cells Numeric value. The minimum no. of cells expressing the splicing event to be included for analysis.
#' @param sigma.sq Numeric value. The variance threshold below which the included/excluded modality will be defined as primary sub-modality, and above which it will be defined as dispersed sub-modality. Please refer to \code{AssignModality} function help page for more details.
#' @param bimodal.adjust Logical. When set to \code{TRUE}, MARVEL will identify false bimodal modalities and reassign them as included/excluded modality. Please refer to \code{AssignModality} function help page for more details.
#' @param seed Numeric value. Ensure the \code{fitdist} function returns the same values for alpha and beta paramters each time this function is executed using the same random number generator. Please refer to \code{AssignModality} function help page for more details.
#' @param modality.column Character string. Can take the value \code{"modality"}, \code{"modality.var"} or \code{"modality.bimodal.adj"}. Please refer to \code{AssignModality} function help page for more details.
#' @export
#' @return An object of class S3 containing all the original slots as inputted by the user in addition to one new slot named \code{$adhocPlot$PSI}.
#' @author Sean Wen <sean.wenwx@gmail.com>
#' @importFrom plyr join
#' @import stats
#' @import methods
#' @import ggplot2
#' @examples
#' marvel <- readRDS(system.file("extdata/Data", "MarvelObject.rds", package="MARVEL"))
#'
#' feature <- marvel$SpliceFeature$SE$tran_id[1]
#'
#' marvel <- PlotValues.PSI(MarvelObject=marvel,
#'                          cell.types=c("iPSC", "Endoderm"),
#'                          feature=feature,
#'                          maintitle="gene_short_name",
#'                          xlabels.size=12,
#'                          n.cells=25,
#'                          sigma.sq=0.001,
#'                          bimodal.adjust=TRUE,
#'                          seed=1,
#'                          modality.column="modality.bimodal.adj"
#'                          )
#'
#' marvel$adhocPlot$PSI

PlotValues.PSI <- function(MarvelObject, cell.types, feature, maintitle, xlabels.size, n.cells, sigma.sq, bimodal.adjust, seed, modality.column) {
    
    # Define arguments
    df <- do.call(rbind.data.frame, MarvelObject$PSI)
    df.pheno <- MarvelObject$SplicePheno
    df.feature <- do.call(rbind.data.frame, MarvelObject$SpliceFeatureValidated)
    cell.types <- cell.types
    feature <- feature
    maintitle <- maintitle
    xlabels.size <- xlabels.size
    n.cells <- n.cells
    sigma.sq <- sigma.sq
    bimodal.adjust <- bimodal.adjust
    seed <- seed
    modality.column <- modality.column
    
    #MarvelObject <- marvel
    #df <- do.call(rbind.data.frame, marvel$PSI)
    #df.pheno <- marvel$SplicePheno
    #df.feature <- do.call(rbind.data.frame, marvel$SpliceFeatureValidated)
    #cell.types <- c("iPSC", "Endoderm")
    #feature <- "chr3:16264769:16264969:-@chr3:16264155:16264229:-@chr3:16257061:16260829"
    #maintitle <- "gene_short_name"
    #n.cells <- 25
    #sigma.sq <- 0.001
    #bimodal.adjust <- TRUE
    #seed <- 1
    #modality.column <- "modality.bimodal.adj"
    
    # Subset relevant feature
    df.feature <- df.feature[which(df.feature$tran_id == feature), ]
    df <- df[which(df$tran_id == feature), ]
    
    # Create ad hoc marvel object
    s3 <- list()
    class(s3) <- "Marvel"
    s3$PSI <- list(df)
    s3$SplicePheno <- df.pheno
    s3$SpliceFeatureValidated <- list(df.feature)

    # Assign modalities
    mod <- NULL
    
    for(i in 1:length(cell.types)) {
                
        . <- AssignModality(MarvelObject=s3,
                       cell.type=cell.types[i],
                       n.cells=n.cells,
                       sigma.sq=sigma.sq,
                       bimodal.adjust=TRUE,
                       seed=1
                       )
                       
        mod[i] <- .$Modality$Results[,modality.column]
        
    
    }
    
    modality <- data.frame("cell.type"=cell.types, "modality"=mod, stringsAsFactors=FALSE)
    
    #######################################################################################
    
    # Create row names for matrix
    row.names(df) <- df$tran_id
    df$tran_id <- NULL
    
    # Subset overlapping samples in matrix and pheno file
    df <- df[, which(names(df) %in% df.pheno$sample.id)]

    # Subset relevant feature
    df.feature.small <- df.feature[which(df.feature$tran_id ==  feature), ]
    df.small <- df[feature, ]
    df.small <- as.data.frame(t(df.small))
    names(df.small) <- "psi"
    df.small$sample.id <- row.names(df.small)
    row.names(df.small) <- NULL
    
    # Subset relevant cell types
    df.small <- join(df.small, df.pheno[,c("sample.id", "cell.type")], by="sample.id", type="left")
    df.small <- df.small[which(df.small$cell.type %in% cell.types), ]
    df.small$cell.type <- factor(df.small$cell.type, levels=cell.types)
    
    # Compute statistics
        # n cells per cell type
        n.cells <- tapply(df.small$psi, df.small$cell.type, function(x) {sum(!is.na(x))})
        n.cells <- data.frame("cell.type"=names(n.cells), "freq"=n.cells, stringsAsFactors=FALSE)
        row.names(n.cells) <- NULL
        
        n.cells$label <- paste(n.cells$cell.type, "\n", "(n=", n.cells$freq, ")", "\n", modality$modality, sep="")
        #n.cells$label <- gsub(".", "\n", n.cells$label, fixed=TRUE)
        
        # Average
        ave <- tapply(df.small$psi, df.small$cell.type, function(x) {mean(x, na.rm=TRUE)})
        ave <- data.frame("cell.type"=names(ave), "average"=ave, stringsAsFactors=FALSE)
        row.names(ave) <- NULL
        
    # Keep complete entries
    df.small <- df.small[which(!is.na(df.small$cell.type)), ]
    df.small <- df.small[which(!is.na(df.small$psi)), ]
    
    # Violin plot
        # Definition
        data <- df.small
        x <- data$cell.type
        y <- data$psi
        z <- data$cell.type
        maintitle <- df.feature.small[, maintitle]
        ytitle <- "PSI"
        xtitle <- ""
        xlabels <- n.cells$label
        fivenum(y) ; ymin <- 0 ; ymax <- 1 ; yinterval <- 0.25

        # Plot
        plot <- ggplot() +
            geom_violin(data, mapping=aes(x=x, y=y, fill=z), color="gray", scale="width") +
            geom_jitter(data, mapping=aes(x=x, y=y), position=position_jitter(width=0.1, height=0), size=0.001) +
            stat_summary(data, mapping=aes(x=x, y=y), geom="point", fun="mean", fill="red", col="black", size=2, shape=23) +
            scale_x_discrete(labels=xlabels) +
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
                axis.title=element_text(size=15),
                axis.text=element_text(size=15),
                axis.text.x=element_text(size=xlabels.size, colour="black"),
                axis.text.y=element_text(size=15, colour="black"),
                legend.position="none",
                legend.title=element_text(size=12),
                legend.text=element_text(size=12)
                )

    # Save into new slow
    MarvelObject$adhocPlot$PSI <- plot
    
    return(MarvelObject)
    
}
