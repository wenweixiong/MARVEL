#' @title Annotate reduced dimension space with PSI values
#'
#' @description Annotates reduced dimension space, e.g., UMAP and tSNE, with PSI values.
#'
#' @param MarvelObject Marvel object. S3 object generated from \code{CheckAlignment.10x} function.
#' @param cell.ids Vector of character strings. Specific set of cells to plot.
#' @param coord.intron Character string. Coordinates of splice junction whose expression will be plotted.
#' @param min.gene.count Numeric value. Minimum raw gene count, above which, the PSI value will be calculate for the cell. Default is \code{3}.
#' @param point.size Numeric value. Size of data points. Default is \code{1}.
#' @param log2.transform Logical value. If set to \code{TRUE}, PSI values will be log2-transformed. Useful for highlighting small changes in PSI values between cell groups. Default is \code{FALSE}.
#' @param color.gradient Vector of character strings. Colors to indicate low, moderate, and high expression. Default is \code{c("grey90","blue","red")}.
#' @param type Character string. Type of reduced dimension space. Options are \code{"umap"} and \code{"tsne"}.
#'
#' @return An object of class S3 with new slot \code{MarvelObject$adhocPlot$PCA$PSI}.
#'
#' @importFrom plyr join
#' @import ggplot2
#' @import Matrix
#'
#' @export
#'
#' @examples
#'
#' marvel.demo.10x <- readRDS(system.file("extdata/data",
#'                                "marvel.demo.10x.rds",
#'                                package="MARVEL")
#'                                )
#'
#' # Define cell groups
#'     # Retrieve sample metadata
#'     sample.metadata <- marvel.demo.10x$sample.metadata
#'
#'     # iPSC
#'     index <- which(sample.metadata$cell.type=="iPSC")
#'     cell.ids.1 <- sample.metadata[index, "cell.id"]
#'     length(cell.ids.1)
#'
#'     # Cardio day 10
#'     index <- which(sample.metadata$cell.type=="Cardio day 10")
#'     cell.ids.2 <- sample.metadata[index, "cell.id"]
#'     length(cell.ids.2)
#'
#'     # Save into list
#'     cell.group.list <- list("iPSC"=cell.ids.1,
#'                             "Cardio d10"=cell.ids.2
#'                             )
#'
#' # Plot expression
#' marvel.demo.10x <- PlotValues.PCA.PSI.10x(
#'                         MarvelObject=marvel.demo.10x,
#'                         coord.intron="chr1:23693914:23694659",
#'                         min.gene.count=3,
#'                         log2.transform=FALSE,
#'                         color.gradient=c("grey","cyan","green","yellow","red"),
#'                         type="tsne"
#'                         )
#'
#' # Check output
#' marvel.demo.10x$adhocPlot$PCA$PSI

PlotValues.PCA.PSI.10x <- function(MarvelObject, cell.ids=NULL, coord.intron, min.gene.count=3, point.size=0.1, log2.transform=FALSE, color.gradient=c("grey90","blue","red"), type) {

    # Example aruguments
    MarvelObject <- MarvelObject
    df.coord <- MarvelObject$pca
    sj.metadata <- MarvelObject$sj.metadata
    df.gene.count <- MarvelObject$gene.count.matrix
    df.sj.count <- MarvelObject$sj.count.matrix
    cell.ids <- cell.ids
    coord.intron <- coord.intron
    min.gene.count <- min.gene.count
    point.size <- point.size
    log2.transform <- log2.transform
    color.gradient <- color.gradient
    type <- type
    
    # Example aruguments
    #MarvelObject <- marvel
    #df.coord <- MarvelObject$pca
    #sj.metadata <- MarvelObject$sj.metadata
    #df.gene.count <- MarvelObject$gene.count.matrix
    #df.sj.count <- MarvelObject$sj.count.matrix
    #cell.ids <- cell.ids
    #coord.intron <- "chr12:78865110:78977798"
    #min.gene.count <- 3
    #point.size <- 0.1
    #log2.transform <- TRUE
    #color.gradient <- c("grey","cyan","green","yellow","red")
    #type <- "tsne"
    
    ##########################################################################
    ####################### ANNOTATE SJ COUNTS ###############################
    ##########################################################################
    
    # Subset relevant SJ
    df.sj.count <- df.sj.count[coord.intron, ]
    
    # Tabulate counts
    df.sj.count  <- data.frame("cell.id"=names(df.sj.count),
                               "sj.count"=as.numeric(df.sj.count),
                               stringsAsFactors=FALSE
                               )
    
    # Annotate coordinate table
    df.coord <- join(df.coord, df.sj.count, by="cell.id", type="left")
    
    ##########################################################################
    ####################### ANNOTATE SJ COUNTS ###############################
    ##########################################################################
    
    # Retrieve gene
    gene_short_name <- sj.metadata[which(sj.metadata$coord.intron==coord.intron), "gene_short_name.start"]
    
    # Subset relevant gene
    df.gene.count <- df.gene.count[gene_short_name, ]
    
    # Tabulate counts
    df.gene.count  <- data.frame("cell.id"=names(df.gene.count),
                                 "gene.count"=as.numeric(df.gene.count),
                                  stringsAsFactors=FALSE
                                  )
    
    # Annotate coordinate table
    df.coord <- join(df.coord, df.gene.count, by="cell.id", type="left")
    
    ##########################################################################
    ############################# COMPUTE PSI ################################
    ##########################################################################
    
    df.coord$psi <- df.coord$sj.count / df.coord$gene.count * 100
    
    ##########################################################################
    ####################### SET THRESHOLD FOR COUNTS #########################
    ##########################################################################
    
    # Define margins before removing cells
    xmin <- min(df.coord$x); xmax <- max(df.coord$x)
    ymin <- min(df.coord$y); ymax <- max(df.coord$y)
    
    # Remove unexpressed genes
    #coord.xy.small <- coord.xy.small[which(coord.xy.small$count.exp >= 1), ]

    #if(min.gene.count >= 2) {
    
        # Remove lowly-expressed genes
        #coord.xy.small <- coord.xy.small[which(coord.xy.small$count.exp >= min.gene.count), ]
        
    #}
    
    # Censor lowly expressed gene
    df.coord$psi[which(df.coord$gene.count < min.gene.count)] <- NA
    
    # Remove non-expressing/lowly expressed gene
    df.coord <- df.coord[!is.na(df.coord$psi), ]
    
    ##########################################################################
    ############################## SUBSET CELLS ##############################
    ##########################################################################
    
    if(!is.null(cell.ids[1])) {
        
        df.coord <- df.coord[which(df.coord$cell.id %in% cell.ids), ]
        
        
    }
    
    ##########################################################################
    ################################ PLOT ####################################
    ##########################################################################
    
    # Cap max PSI
    df.coord$psi[which(df.coord$psi > 100)] <- 100
    
    # Reorder by expression
    df.coord <- df.coord[order(df.coord$psi), ]
        
    # Definitions
    data <- df.coord
    x <- data$x
    y <- data$y
    
    if(log2.transform==TRUE) {
        
        z <- log2(data$psi + 1)
        
    } else {
    
        z <- data$psi
    
    }
    
    
    if(type=="umap"){
        
        xtitle <- "UMAP-1"
        ytitle <- "UMAP-2"
        
    } else if(type=="tsne"){
    
        xtitle <- "tSNE-1"
        ytitle <- "tSNE-2"
    
    }
    
    if(log2.transform==TRUE) {
        
        legendtitle <- "log2(PSI)"
        
    } else {
    
        legendtitle <- "PSI"
    
    }
            
    # Plot
    plot <- ggplot() +
        geom_point(data, mapping=aes(x=x, y=y, color=z), size=point.size) +
        scale_color_gradientn(colors=color.gradient) +
        scale_x_continuous(limits=c(xmin, xmax)) +
        scale_y_continuous(limits=c(ymin, ymax)) +
        labs(title=NULL, x=xtitle, y=ytitle, color=legendtitle) +
        theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            panel.border=element_blank(),
            plot.title=element_text(hjust = 0.5, size=15),
            plot.subtitle=element_text(hjust = 0.5, size=15),
            axis.line.y.left = element_line(color="black"),
            axis.line.x = element_line(color="black"),
            axis.title=element_text(size=12),
            axis.ticks.x=element_blank(),
            axis.ticks.y=element_blank(),
            axis.text.x=element_blank(),
            axis.text.y=element_blank(),
            #legend.position="none",
            legend.title=element_text(size=8),
            legend.text=element_text(size=8)
            )

    ##########################################################################

    # Save into new slot
    MarvelObject$adhocPlot$PCA$PSI <- plot

    # Return final object
    return(MarvelObject)

}
