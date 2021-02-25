#' @title Classify Modality Changes
#'
#' @description
#' \code{ModalityChange} Classifies the type of modality change for each splicing event that has taken place between 2 groups of cells.
#'
#' @details
#' This function classifies the type of modality change for each splicing event that has taken place between 2 groups of cells. Explicit: When modality changes between one of the five main modalities, e.g. included to multimodal. Implicit: When modality changes between primary and dispersed sub-modalities, e.g. included-primary to included-dispersed. Restricted: No modality change, e.g. included to included.
#'
#' @param MarvelObject S3 object generated from \code{CompareValues} function.
#' @param psi.de.sig Numeric value. Adjusted p-value below which the splicing event is considered differentially spliced and included for isoform switching analysis.
#' @param cell.types Character string. To indicate which 2 groups of cells that will be used for differential splicing analysis. Group names should match those in \code{cell.type} column of \code{$SplicePheno} slot.
#' @param n.cells Numeric value. The minimum no. of cells expressing the splicing event for the event to be included for differential splicing analysis. Please refer to \code{AssignModality} function help page for more details.
#' @param sigma.sq Numeric value. The variance threshold below which the included/excluded modality will be defined as primary sub-modality, and above which it will be defined as dispersed sub-modality. Please refer to \code{AssignModality} function help page for more details.
#' @param bimodal.adjust Logical. When set to \code{TRUE}, MARVEL will identify false bimodal modalities and reassign them as included/excluded modality. Please refer to \code{AssignModality} function help page for more details.
#' @param seed Numeric value. Ensure the \code{fitdist} function returns the same values for alpha and beta paramters each time this function is executed using the same random number generator. Please refer to \code{AssignModality} function help page for more details.
#' @param modality.column Character string. Can take the value \code{"modality"}, \code{"modality.var"} or \code{"modality.bimodal.adj"}. Please refer to \code{AssignModality} function help page for more details.
#' @export
#' @return An object of class S3 containing all the original slots as inputted by the user in addition to three new slots named \code{$DE$Modality*}. \code{$DE$Modality} Original data frame generated from \code{CompareValues} function with an additional columns to indicate the type of modality changes that have taken place between the 2 groups of cells. \code{$DE$ModalityProp} Tabulated proportion for each type of modality change. \code{$DE$ModalityPlot} Doughnut plot representing the values in \code{$DE$ModalityProp}.
#' @author Sean Wen <sean.wenwx@gmail.com>
#' @importFrom plyr join
#' @import stats
#' @import methods
#' @import ggplot2
#' @examples
#' marvel <- readRDS(system.file("extdata/Data", "MarvelObject.rds", package="MARVEL"))
#'
#' marvel <- ModalityChange(MarvelObject=marvel,
#'                          psi.de.sig=0.05,
#'                          cell.types=c("iPSC", "Endoderm"),
#'                          n.cells=25,
#'                          sigma.sq=0.001,
#'                          bimodal.adjust=TRUE,
#'                          seed=1,
#'                          modality.column="modality.bimodal.adj"
#'                          )
#'
#' marvel$DE$ModalityProp
#' marvel$DE$ModalityPlot

ModalityChange <- function(MarvelObject, psi.de.sig, cell.types, n.cells, sigma.sq, bimodal.adjust, seed, modality.column) {
    
    # Define arguments
    de <- MarvelObject$DE$PSI
    df <- do.call(rbind.data.frame, MarvelObject$PSI)
    df.pheno <- MarvelObject$SplicePheno
    df.feature <- do.call(rbind.data.frame, MarvelObject$SpliceFeatureValidated)
    psi.de.sig <- psi.de.sig
    cell.types <- cell.types
    n.cells <- n.cells
    sigma.sq <- sigma.sq
    bimodal.adjust <- bimodal.adjust
    seed <- seed
    modality.column <- modality.column
    
    #MarvelObject <- marvel
    #de <- marvel$DE$PSI
    #df <- do.call(rbind.data.frame, marvel$PSI)
    #df.pheno <- marvel$SplicePheno
    #df.feature <- do.call(rbind.data.frame, marvel$SpliceFeatureValidated)
    #psi.de.sig <- 0.05
    #cell.types <- c("iPSC", "Endoderm")
    #n.cells <- 25
    #sigma.sq <- 0.001
    #bimodal.adjust <- TRUE
    #seed <- 1
    #modality.column <- "modality.bimodal.adj"

    # Subset overlapping samples in matrix and pheno file
    df <- df[, c(1, which(names(df) %in% df.pheno$sample.id))]
    
    # Subset relevant feature
    tran_ids <- de[which(de$p.val.adj < psi.de.sig), "tran_id"]
    df.feature <- df.feature[which(df.feature$tran_id %in% tran_ids), ]
    df <- df[which(df$tran_id %in% df.feature$tran_id), ]
    
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
                       
        mod[[i]] <- .$Modality$Results
        
    
    }
    
    results <- data.frame("tran_id"=mod[[1]]$tran_id,
                      "modality.g1"=mod[[1]][,modality.column],
                      "modality.g2"=mod[[2]][,modality.column],
                      stringsAsFactors=FALSE
                      )

    # Indicate modality change
        # Restricted
        results$modality.change <- NA
        index.na <- which(is.na(results$modality.change))
        index.new <- which(results$modality.g1==results$modality.g2)
        index <- intersect(index.na, index.new)
        results$modality.change[index] <- "Restricted"
        
        # Implicit
            # Included
            index.na <- which(is.na(results$modality.change))
            index.new <- intersect(grep("^Included", results$modality.g1),
                                   grep("^Included", results$modality.g2)
                                   )
            index <- intersect(index.na, index.new)
            results$modality.change[index] <- "Implicit"
            
            # Included
            index.na <- which(is.na(results$modality.change))
            index.new <- intersect(grep("^Excluded", results$modality.g1),
                                   grep("^Excluded", results$modality.g2)
                                   )
            index <- intersect(index.na, index.new)
            results$modality.change[index] <- "Implicit"
    
        # Explicit
        index.na <- which(is.na(results$modality.change))
        index <- index.na
        results$modality.change[index] <- "Explicit"

    # Annotate DE data frame
    de <- join(de, results, by="tran_id", type="left")
    de <- de[which(!is.na(de$modality.change)), ]
    
    # Doughnut plot
        # Tabulate freq
        . <- as.data.frame(table(de$modality.change), stringsAsFactors=FALSE)
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
    
    # Save into new slow
    MarvelObject$DE$Modality <- de
    MarvelObject$DE$ModalityPlot <- plot
    MarvelObject$DE$ModalityProp <- .[,c("modality.change", "freq", "pct")]
    
    return(MarvelObject)
    
}
