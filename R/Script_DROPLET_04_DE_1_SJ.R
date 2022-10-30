#' @title Differential splice junction analysis
#'
#' @description Performs differential splice junction analysis between two groups of cells.
#'
#' @param MarvelObject Marvel object. S3 object generated from \code{CheckAlignment.10x} function.
#' @param coord.introns Character strings. Specific splice junctions to be included for analysis. Default is \code{NULL}.
#' @param cell.group.g1 Vector of Character strings. Cell IDs corresponding to Group 1 (reference group).
#' @param cell.group.g2 Vector of Character strings. Cell IDs corresponding to Group 2.
#' @param min.pct.cells.genes Numeric value. Minimum percentage of cells in which the gene is expressed for that gene to be included for splice junction expression distribution analysis. Expressed genes defined as genes with non-zero normalised UMI counts. This threshold may be determined from \code{PlotPctExprCells.SJ.10x} function. Default is \code{10}.
#' @param min.pct.cells.sj Numeric value. Minimum percentage of cells in which the splice junction is expressed for that splice junction to be included for splice junction expression distribution analysis. Expressed splice junctions defined as splice junctions with raw UMI counts >= 1. This threshold may be determined from \code{PlotPctExprCells.SJ.10x} function. Default is \code{10}.
#' @param min.gene.norm Numeric value. The average normalised gene expression across the two cell groups above which the splice junction will be included for analysis. Default is \code{1.0}.
#' @param seed Numeric value. Random number generator to be fixed for permutations test and down-sampling.
#' @param n.iterations Numeric value. Number of times to shuffle the cell group labels when building the null distribution. Default is \code{100}.
#' @param downsample Logical value. If set to \code{TRUE}, both cell groups will be down-sampled so that both cell groups will have the same number of cells. The number of cells to downsample will be based on the smallest cell group. Default is \code{FALSE}.
#' @param show.progress Logical value. If set to \code{TRUE} (default), the progress bar will appear.
#'
#' @return An object of class S3 with a new slots \code{MarvelObject$DE$SJ$Table}, \code{MarvelObject$DE$SJ$cell.group.g1}, and \code{MarvelObject$DE$SJ$cell.group.g2}.
#'
#' @importFrom plyr join
#' @import Matrix
#' @importFrom utils txtProgressBar setTxtProgressBar
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
#'     # Group 1 (reference)
#'     index <- which(sample.metadata$cell.type=="iPSC")
#'     cell.ids.1 <- sample.metadata[index, "cell.id"]
#'     length(cell.ids.1)
#'
#'     # Group 2
#'     index <- which(sample.metadata$cell.type=="Cardio day 10")
#'     cell.ids.2 <- sample.metadata[index, "cell.id"]
#'     length(cell.ids.2)
#'
#' # DE
#' marvel.demo.10x <- CompareValues.SJ.10x(
#'                         MarvelObject=marvel.demo.10x,
#'                         cell.group.g1=cell.ids.1,
#'                         cell.group.g2=cell.ids.2,
#'                         min.pct.cells.genes=10,
#'                         min.pct.cells.sj=10,
#'                         min.gene.norm=1.0,
#'                         seed=1,
#'                         n.iterations=100,
#'                         downsample=TRUE,
#'                         show.progress=FALSE
#'                         )
#'
#' # Check output
#' head(marvel.demo.10x$DE$SJ$Table)

CompareValues.SJ.10x <- function(MarvelObject, coord.introns=NULL, cell.group.g1, cell.group.g2, min.pct.cells.genes=10, min.pct.cells.sj=10, min.gene.norm=1.0, seed=1, n.iterations=100, downsample=FALSE, show.progress=TRUE) {
        
    # Define arguments
    MarvelObject <- MarvelObject
    sample.metadata <- MarvelObject$sample.metadata
    sj.metadata <- MarvelObject$sj.metadata
    df.gene.norm <- MarvelObject$gene.norm.matrix
    df.gene.count <- MarvelObject$gene.count.matrix
    df.sj.count <- MarvelObject$sj.count.matrix
    coord.introns.custom <- coord.introns
    cell.group.g1 <- cell.group.g1
    cell.group.g2 <- cell.group.g2
    min.pct.cells.genes <- min.pct.cells.genes
    min.pct.cells.sj <- min.pct.cells.sj
    seed <- seed
    n.iterations <- n.iterations
    downsample <- downsample
    show.progress <- show.progress
    min.gene.norm <- min.gene.norm
    
    # Example arguments
    #MarvelObject <- marvel
    #sample.metadata <- MarvelObject$sample.metadata
    #sj.metadata <- MarvelObject$sj.metadata
    #df.gene.norm <- MarvelObject$gene.norm.matrix
    #df.gene.count <- MarvelObject$gene.count.matrix
    #df.sj.count <- MarvelObject$sj.count.matrix
    #coord.introns.custom <- df$coord.intron
    #cell.group.g1 <- cell.group.g1
    #cell.group.g2 <- cell.group.g2
    #min.pct.cells.genes <- 10
    #min.pct.cells.sj <- 0.1
    #seed <- 1
    #n.iterations <- 100
    #downsample <- FALSE
    #min.gene.norm <- 0.1
    #show.progress <- TRUE
    
    #################################################################
    ################### SUBSET SPECIFC SJs ##########################
    #################################################################
    
    if(!is.null(coord.introns.custom[1])) {
        
        # Subset SJ
            # Find overlap
            overlap <- intersect(coord.introns.custom, row.names(df.sj.count))
            
            # Subset
            df.sj.count <- df.sj.count[overlap,]
            sj.metadata <- sj.metadata[which(sj.metadata$coord.intron %in% overlap), ]
            
        # Subset genes
            # Find overlap
            gene_short_names <- unique(sj.metadata$gene_short_name.start)
            
            # Subset
            df.gene.norm <- df.gene.norm[gene_short_names, ]
            df.gene.count <- df.gene.count[gene_short_names, ]
            
        # Report progress
        message(paste(length(coord.introns.custom), " SJs specified by user", sep=""))
        message(paste(length(overlap), " overlapping SJs found and subset-ed", sep=""))
        message(paste(length(gene_short_names), " corresponding genes found and subset-ed", sep=""))

            
    }
    
    #################################################################
    ######################## DOWNSAMPLE #############################
    #################################################################
    
    if(downsample==TRUE) {
        
        # Set random number generator
        set.seed(seed)
        
        # Find lowest common denominator
        n.cells.downsample <- min(length(cell.group.g1), length(cell.group.g2))
        
        # Downsample
        cell.group.g1 <- sample(cell.group.g1, size=n.cells.downsample, replace=FALSE)
        cell.group.g2 <- sample(cell.group.g2, size=n.cells.downsample, replace=FALSE)
        
    }
    
    # Report progress
    message(paste(length(cell.group.g1), " cells from Group 1 and ", length(cell.group.g2), " cells from Group 2 included", sep=""))
    
    # Subset matrices
    df.gene.norm <- df.gene.norm[, c(cell.group.g1, cell.group.g2)]
    df.gene.count <- df.gene.count[, c(cell.group.g1, cell.group.g2)]
    df.sj.count <- df.sj.count[, c(cell.group.g1, cell.group.g2)]
    
    # Check alignment
    table(colnames(df.gene.norm)==colnames(df.gene.count))
    table(colnames(df.gene.count)==colnames(df.sj.count))

    #################################################################
    ################## SUBSET EXPRESSED GENES (1) ###################
    #################################################################
    
    # Compute num. of cells in which gene is expressed: Group 1
        # Subset cells
        df.gene.norm.small <- df.gene.norm[, cell.group.g1]
        
        # Compute n cells express
        . <- apply(df.gene.norm.small, 1, function(x) { sum(x != 0)})
        . <- data.frame("cell.group"="cell.group.g1",
                        "gene_short_name"=names(.),
                        "n.cells.total"=length(cell.group.g1),
                        "n.cells.expr"=as.numeric(.),
                        "pct.cells.expr"=round(as.numeric(.)/length(cell.group.g1) * 100, digits=2),
                        stringsAsFactors=FALSE
                        )
        
        # Save as new object
        results.g1 <- .
    
    # Compute num. of cells in which gene is expressed: Group 2
        # Subset cells
        df.gene.norm.small <- df.gene.norm[, cell.group.g2]
        
        # Compute n cells express
        . <- apply(df.gene.norm.small, 1, function(x) { sum(x != 0)})
        . <- data.frame("cell.group"="cell.group.g2",
                        "gene_short_name"=names(.),
                        "n.cells.total"=length(cell.group.g2),
                        "n.cells.expr"=as.numeric(.),
                        "pct.cells.expr"=round(as.numeric(.)/length(cell.group.g2) * 100, digits=2),
                        stringsAsFactors=FALSE
                        )
        
        # Save as new object
        results.g2 <- .
        
    # Merge
    results <- rbind.data.frame(results.g1, results.g2)
    results$cell.group <- factor(results$cell.group, levels=c("cell.group.g1", "cell.group.g2"))
    
    # Subset expressed genes
    results <- results[which(results$pct.cells.expr > min.pct.cells.genes), ]
    
    # Subset genes expressing in both cell groups
    gene_short.names.1 <- results[which(results$cell.group=="cell.group.g1"), "gene_short_name"]
    gene_short.names.2 <- results[which(results$cell.group=="cell.group.g2"), "gene_short_name"]
    gene_short_names <- intersect(gene_short.names.1, gene_short.names.2)
    
    # Report progress
    message(paste(length(gene_short.names.1), " genes expressed in cell group 1", sep=""))
    message(paste(length(gene_short.names.2), " genes expressed in cell group 2", sep=""))
    message(paste(length(gene_short_names), " genes expressed in BOTH cell group and retained", sep=""))
    
    #################################################################
    ################## SUBSET EXPRESSED GENES (2) ###################
    #################################################################
    
    # Subset expressed genes from part 1
    df.gene.norm.small <- df.gene.norm[gene_short_names, ]
    
    # Compute combined average
    . <- apply(df.gene.norm.small, 1, function(x) {mean(log2(x + 1))})
    mean.combined.df <- data.frame("gene_short_name"=names(.),
                                   "mean.expr.gene.norm.g1.g2"=as.numeric(.),
                                   stringsAsFactors=FALSE
                                   )
    
    # Subset expressed genes
    index <- which(mean.combined.df$mean.expr.gene.norm.g1.g2 > min.gene.norm)
    mean.combined.df <- mean.combined.df[index, ]
    gene_short_names <- mean.combined.df$gene_short_name
    
    # Report progress
    message(paste(length(gene_short_names), " genes with mean log2(expression + 1) > ", min.gene.norm, " retained", sep=""))
    
    #################################################################
    ###################### SUBSET EXPRESSED SJ ######################
    #################################################################
    
    # Compute num. of cells in which SJ is expressed: Group 1
        # Subset cells
        df.sj.count.small <- df.sj.count[, cell.group.g1]
        
        # Subset expressed genes
        coord.introns <- sj.metadata[which(sj.metadata$gene_short_name.start %in% gene_short_names), "coord.intron"]
        length(coord.introns)
        df.sj.count.small <- df.sj.count.small[coord.introns, ]
        
        # Compute n cells express
        . <- apply(df.sj.count.small, 1, function(x) { sum(x != 0)})
        . <- data.frame("cell.group"="cell.group.g1",
                        "coord.intron"=names(.),
                        "n.cells.total"=length(cell.group.g1),
                        "n.cells.expr"=as.numeric(.),
                        "pct.cells.expr"=round(as.numeric(.)/length(cell.group.g1) * 100, digits=2),
                        stringsAsFactors=FALSE
                        )
        
        # Save as new object
        results.g1 <- .

    # Compute num. of cells in which SJ is expressed: Group 2
        # Subset cells
        df.sj.count.small <- df.sj.count[, cell.group.g2]
        
        # Subset expressed genes
        coord.introns <- sj.metadata[which(sj.metadata$gene_short_name.start %in% gene_short_names), "coord.intron"]
        length(coord.introns)
        df.sj.count.small <- df.sj.count.small[coord.introns, ]
        
        # Compute n cells express
        . <- apply(df.sj.count.small, 1, function(x) { sum(x != 0)})
        . <- data.frame("cell.group"="cell.group.g2",
                        "coord.intron"=names(.),
                        "n.cells.total"=length(cell.group.g2),
                        "n.cells.expr"=as.numeric(.),
                        "pct.cells.expr"=round(as.numeric(.)/length(cell.group.g2) * 100, digits=2),
                        stringsAsFactors=FALSE
                        )
        
        # Save as new object
        results.g2 <- .
            
    # Merge
    results <- rbind.data.frame(results.g1, results.g2)
    results$cell.group <- factor(results$cell.group, levels=c("cell.group.g1", "cell.group.g2"))
    
    # Subset expressed SJ
    results <- results[which(results$pct.cells.expr > min.pct.cells.sj), ]
    
    # Subset SJ expressing in either cell groups
    coord.introns.1 <- results[which(results$cell.group=="cell.group.g1"), "coord.intron"]
    coord.introns.2 <- results[which(results$cell.group=="cell.group.g2"), "coord.intron"]
    coord.introns <- unique(c(coord.introns.1, coord.introns.2))
    
    # Report progress
    message(paste(length(coord.introns.1), " SJ expressed in cell group 1", sep=""))
    message(paste(length(coord.introns.2), " SJ expressed in cell group 2", sep=""))
    message(paste(length(coord.introns), " SJ expressed in EITHER cell groups and retained", sep=""))
    
    # Report final numbers
    n.sj <- length(coord.introns)
    n.genes <- length(unique(sj.metadata[which(sj.metadata$coord.intron %in% coord.introns), "gene_short_name.start"]))
    message(paste("Total of ", n.sj, " SJ from ", n.genes, " genes included for DE analysis", sep=""))
        
    #################################################################
    ############## SUBSET EXPRESSED GENES, SJ, CELLS ################
    #################################################################
    
    # Subset cells, expressed genes, SJs
        # Metadata
        sj.metadata <- sj.metadata[which(sj.metadata$coord.intron %in% coord.introns), ]
    
        # Group 1
        df.sj.count.g1 <- df.sj.count[sj.metadata$coord.intron, cell.group.g1]
        df.gene.count.g1 <- df.gene.count[unique(sj.metadata$gene_short_name.start), cell.group.g1]
        
        # Group 2
        df.sj.count.g2 <- df.sj.count[sj.metadata$coord.intron, cell.group.g2]
        df.gene.count.g2 <- df.gene.count[unique(sj.metadata$gene_short_name.start), cell.group.g2]
        
        # Check alignment
        table(row.names(df.sj.count.g1)==row.names(df.sj.count.g2))
        table(row.names(df.gene.count.g1)==row.names(df.gene.count.g2))

    #################################################################
    ######################### COMPUTE PSI ###########################
    #################################################################
    
    # Compute PSI: Group 1
        # Report progress
        message("Computing PSI for cell group 1...")
        
        # Compute cell group size
        n.cells.total <- ncol(df.sj.count.g1)
        
        # Tabulate SJ metrices
            # Compute % expressed SJ
            n.cells.expr.sj <- apply(df.sj.count.g1, 1, function(x) {sum(x != 0)})
            pct.cells.expr.sj <- round(n.cells.expr.sj/n.cells.total * 100, digits=2)
            
            # Compute total SJ counts
            sj.count.total <- apply(df.sj.count.g1, 1, function(x) {sum(x)})
            
            # Save into data frame
            results.sj <- data.frame("coord.intron"=names(sj.count.total),
                                     "n.cells.total"=n.cells.total,
                                     "n.cells.expr.sj"=n.cells.expr.sj,
                                     "pct.cells.expr.sj"=pct.cells.expr.sj,
                                     "sj.count.total"=sj.count.total,
                                     stringsAsFactors=FALSE
                                     )
            row.names(results.sj) <- NULL
            
            # Annotate gene
            results.sj <- join(results.sj, sj.metadata[, c("coord.intron", "gene_short_name.start")], by="coord.intron", type="left")
            
        # Tabulate gene metrices
            # Compute % expressed genes
            n.cells.expr.gene <- apply(df.gene.count.g1, 1, function(x) {sum(x != 0)})
            pct.cells.expr.gene <- round(n.cells.expr.gene/n.cells.total * 100, digits=2)
            
            # Compute total gene counts
            gene.count.total <- apply(df.gene.count.g1, 1, function(x) {sum(x)})
        
            # Tabulate results
            results.gene <- data.frame("gene_short_name.start"=names(gene.count.total),
                                       "n.cells.expr.gene"=n.cells.expr.gene,
                                       "pct.cells.expr.gene"=pct.cells.expr.gene,
                                       "gene.count.total"=gene.count.total,
                                       stringsAsFactors=FALSE
                                       )
                                       
        # Merge
        results <- join(results.sj, results.gene, by="gene_short_name.start", type="left")
        
        # Compute PSI
        results$psi <- round(results$sj.count.total / results$gene.count.total * 100, digits=2)
        
        # Reorder columns
        cols <- c("coord.intron", "gene_short_name.start", "n.cells.total",
                  "n.cells.expr.sj", "pct.cells.expr.sj",
                  "n.cells.expr.gene", "pct.cells.expr.gene",
                  "sj.count.total", "gene.count.total", "psi"
                  )
        
        results <- results[, cols]
        
        names(results)[which(names(results)=="gene_short_name.start")] <- "gene_short_name"
        
        # Indicate cell group
        names(results)[-which(names(results) %in% c("coord.intron", "gene_short_name"))] <- paste(names(results)[-which(names(results) %in% c("coord.intron", "gene_short_name"))], ".g1", sep="")
     
         # Save as new object
         results.g1 <- results
     
    # Compute PSI: Group 2
        # Report progress
        message("Computing PSI for cell group 2...")
        
        # Compute cell group size
        n.cells.total <- ncol(df.sj.count.g2)
        
        # Tabulate SJ metrices
            # Compute % expressed SJ
            n.cells.expr.sj <- apply(df.sj.count.g2, 1, function(x) {sum(x != 0)})
            pct.cells.expr.sj <- round(n.cells.expr.sj/n.cells.total * 100, digits=2)
            
            # Compute total SJ counts
            sj.count.total <- apply(df.sj.count.g2, 1, function(x) {sum(x)})
            
            # Save into data frame
            results.sj <- data.frame("coord.intron"=names(sj.count.total),
                                     "n.cells.total"=n.cells.total,
                                     "n.cells.expr.sj"=n.cells.expr.sj,
                                     "pct.cells.expr.sj"=pct.cells.expr.sj,
                                     "sj.count.total"=sj.count.total,
                                     stringsAsFactors=FALSE
                                     )
            row.names(results.sj) <- NULL
            
            # Annotate gene
            results.sj <- join(results.sj, sj.metadata[, c("coord.intron", "gene_short_name.start")], by="coord.intron", type="left")
            
        # Tabulate gene metrices
            # Compute % expressed genes
            n.cells.expr.gene <- apply(df.gene.count.g2, 1, function(x) {sum(x != 0)})
            pct.cells.expr.gene <- round(n.cells.expr.gene/n.cells.total * 100, digits=2)
            
            # Compute total gene counts
            gene.count.total <- apply(df.gene.count.g2, 1, function(x) {sum(x)})
        
            # Tabulate results
            results.gene <- data.frame("gene_short_name.start"=names(gene.count.total),
                                       "n.cells.expr.gene"=n.cells.expr.gene,
                                       "pct.cells.expr.gene"=pct.cells.expr.gene,
                                       "gene.count.total"=gene.count.total,
                                       stringsAsFactors=FALSE
                                       )
                                       
        # Merge
        results <- join(results.sj, results.gene, by="gene_short_name.start", type="left")
        
        # Compute PSI
        results$psi <- round(results$sj.count.total / results$gene.count.total * 100, digits=2)
        
        # Reorder columns
        cols <- c("coord.intron", "gene_short_name.start", "n.cells.total",
                  "n.cells.expr.sj", "pct.cells.expr.sj",
                  "n.cells.expr.gene", "pct.cells.expr.gene",
                  "sj.count.total", "gene.count.total", "psi"
                  )
        
        results <- results[, cols]
        
        names(results)[which(names(results)=="gene_short_name.start")] <- "gene_short_name"
        
        # Indicate cell group
        names(results)[-which(names(results) %in% c("coord.intron", "gene_short_name"))] <- paste(names(results)[-which(names(results) %in% c("coord.intron", "gene_short_name"))], ".g2", sep="")
        
        # Save as new object
        results.g2 <- results
        
    # Merge group 1, 2
    table(results.g1$coord.intron==results.g2$coord.intron)
    table(results.g1$gene_short_name==results.g2$gene_short_name)
    
    
    index.l <- table(results.g1$coord.intron==results.g2$coord.intron)
    index.true <- length(which(names(index.l)==TRUE))
    index.false <- length(which(names(index.l)==FALSE))
     
    if(index.true==1 & index.false==0) {
    
        results.g2$coord.intron <- NULL
        results.g2$gene_short_name <- NULL
        results <- cbind.data.frame(results.g1, results.g2)
        
    } else {
        
        message("Error in merging tables from Group 1 and 2")

    }
    
    # Compute log2fc, delta
    results$log2fc <- log2( (results$psi.g2 + 1) / (results$psi.g1 + 1) )
    results$delta <- results$psi.g2 - results$psi.g1
    
    # Save as new object
    results.obs <- results
    
    #################################################################
    ##################### COMPUTE P-VALUES ##########################
    #################################################################
    
    # Report progress
    message("Creating null distributions...")
     
    # Set random num. generator
    set.seed(seed)
    
    if(show.progress==TRUE) {
        
        pb <- txtProgressBar(1, n.iterations, style=3)
        
    }
    
    .list.results.perm <- list()
    
    for(i in 1:n.iterations) {
        
        # Create null distribution
            # Shuffle cell ids
            cell.ids.shuffled <- sample(colnames(df.sj.count), size=ncol(df.sj.count), replace=FALSE)
            
            # Shuffle SJ matrix
            df.sj.count.shuffled <- df.sj.count
            colnames(df.sj.count.shuffled) <- cell.ids.shuffled
            
            # Match column in gene matrix
            df.gene.count.shuffled <- df.gene.count
            colnames(df.gene.count.shuffled) <- cell.ids.shuffled
            
            # Check alignment
            table(colnames(df.sj.count.shuffled)==colnames(df.gene.count.shuffled))
            
            # Report progress
            index.l <- table(colnames(df.sj.count.shuffled)==colnames(df.gene.count.shuffled))
            index.true <- length(which(names(index.l)==TRUE))
            index.false <- length(which(names(index.l)==FALSE))
             
            if(index.true==1 & index.false==0) {
            
                #message(paste("Iteration ", 1, " ...", sep=""))
                
            } else {
                
                return(message(paste("Error in iteration ", 1, " ...", sep="")))

            }
        
        # Subset cells, expressed genes, SJs
            # Group 1
            df.sj.count.g1 <- df.sj.count.shuffled[results.obs$coord.intron, cell.group.g1]
            df.gene.count.g1 <- df.gene.count.shuffled[unique(results.obs$gene_short_name), cell.group.g1]
            
            # Group 2
            df.sj.count.g2 <- df.sj.count.shuffled[results.obs$coord.intron, cell.group.g2]
            df.gene.count.g2 <- df.gene.count.shuffled[unique(results.obs$gene_short_name), cell.group.g2]
            
            # Check alignment
            table(row.names(df.sj.count.g1)==row.names(df.sj.count.g2))
            table(row.names(df.gene.count.g1)==row.names(df.gene.count.g2))
            
        # Compute PSI: Group 1
            # Report progress
            #message("Computing PSI for cell group 1...")
            
            # Compute cell group size
            #n.cells.total <- ncol(df.sj.count.g1)
            
            # Tabulate SJ metrices
                # Compute % expressed SJ
                #n.cells.expr.sj <- apply(df.sj.count.g1, 1, function(x) {sum(x != 0)})
                #pct.cells.expr.sj <- round(n.cells.expr.sj/n.cells.total * 100, digits=2)
                
                # Compute total SJ counts
                sj.count.total <- apply(df.sj.count.g1, 1, function(x) {sum(x)})
                
                # Save into data frame
                results.sj <- data.frame("coord.intron"=names(sj.count.total),
                                         #"n.cells.total"=n.cells.total,
                                         #"n.cells.expr.sj"=n.cells.expr.sj,
                                         #"pct.cells.expr.sj"=pct.cells.expr.sj,
                                         "sj.count.total"=sj.count.total,
                                         stringsAsFactors=FALSE
                                         )
                row.names(results.sj) <- NULL
                
                # Annotate gene
                results.sj <- join(results.sj, sj.metadata[, c("coord.intron", "gene_short_name.start")], by="coord.intron", type="left")
                
            # Tabulate gene metrices
                # Compute % expressed genes
                #n.cells.expr.gene <- apply(df.gene.count.g1, 1, function(x) {sum(x != 0)})
                #pct.cells.expr.gene <- round(n.cells.expr.gene/n.cells.total * 100, digits=2)
                
                # Compute total gene counts
                gene.count.total <- apply(df.gene.count.g1, 1, function(x) {sum(x)})
            
                # Tabulate results
                results.gene <- data.frame("gene_short_name.start"=names(gene.count.total),
                                           #"n.cells.expr.gene"=n.cells.expr.gene,
                                           #"pct.cells.expr.gene"=pct.cells.expr.gene,
                                           "gene.count.total"=gene.count.total,
                                           stringsAsFactors=FALSE
                                           )
                                           
            # Merge
            results <- join(results.sj, results.gene, by="gene_short_name.start", type="left")
            
            # Compute PSI
            results$psi <- round(results$sj.count.total / results$gene.count.total * 100, digits=2)
            
            # Reorder columns
            #cols <- c("coord.intron", "gene_short_name.start", "n.cells.total",
                      #"n.cells.expr.sj", "pct.cells.expr.sj",
                      #"n.cells.expr.gene", "pct.cells.expr.gene",
                      #"sj.count.total", "gene.count.total", "psi"
                      #)
                                  
            #results <- results[, cols]
            
            names(results)[which(names(results)=="gene_short_name.start")] <- "gene_short_name"
            
            # Indicate cell group
            names(results)[-which(names(results) %in% c("coord.intron", "gene_short_name"))] <- paste(names(results)[-which(names(results) %in% c("coord.intron", "gene_short_name"))], ".g1", sep="")
         
            # Save as new object
            results.g1 <- results
            #results.g1 <- results[, c("coord.intron", "psi.g1")]

        # Compute PSI: Group 2
            # Report progress
            #message("Computing PSI for cell group 2...")
            
            # Compute cell group size
            #n.cells.total <- ncol(df.sj.count.g2)
            
            # Tabulate SJ metrices
                # Compute % expressed SJ
                #n.cells.expr.sj <- apply(df.sj.count.g2, 1, function(x) {sum(x != 0)})
                #pct.cells.expr.sj <- round(n.cells.expr.sj/n.cells.total * 100, digits=2)
                
                # Compute total SJ counts
                sj.count.total <- apply(df.sj.count.g2, 1, function(x) {sum(x)})
                
                # Save into data frame
                results.sj <- data.frame("coord.intron"=names(sj.count.total),
                                         #"n.cells.total"=n.cells.total,
                                         #"n.cells.expr.sj"=n.cells.expr.sj,
                                         #"pct.cells.expr.sj"=pct.cells.expr.sj,
                                         "sj.count.total"=sj.count.total,
                                         stringsAsFactors=FALSE
                                         )
                row.names(results.sj) <- NULL
                
                # Annotate gene
                results.sj <- join(results.sj, sj.metadata[, c("coord.intron", "gene_short_name.start")], by="coord.intron", type="left")
                
            # Tabulate gene metrices
                # Compute % expressed genes
                #n.cells.expr.gene <- apply(df.gene.count.g2, 1, function(x) {sum(x != 0)})
                #pct.cells.expr.gene <- round(n.cells.expr.gene/n.cells.total * 100, digits=2)
                
                # Compute total gene counts
                gene.count.total <- apply(df.gene.count.g2, 1, function(x) {sum(x)})
            
                # Tabulate results
                results.gene <- data.frame("gene_short_name.start"=names(gene.count.total),
                                           #"n.cells.expr.gene"=n.cells.expr.gene,
                                           #"pct.cells.expr.gene"=pct.cells.expr.gene,
                                           "gene.count.total"=gene.count.total,
                                           stringsAsFactors=FALSE
                                           )
                                           
            # Merge
            results <- join(results.sj, results.gene, by="gene_short_name.start", type="left")
            
            # Compute PSI
            results$psi <- round(results$sj.count.total / results$gene.count.total * 100, digits=2)
            
            # Reorder columns
            #cols <- c("coord.intron", "gene_short_name.start", "n.cells.total",
                      #"n.cells.expr.sj", "pct.cells.expr.sj",
                      #"n.cells.expr.gene", "pct.cells.expr.gene",
                      #"sj.count.total", "gene.count.total", "psi"
                      #)
            
            #results <- results[, cols]
            
            names(results)[which(names(results)=="gene_short_name.start")] <- "gene_short_name"
            
            # Indicate cell group
            names(results)[-which(names(results) %in% c("coord.intron", "gene_short_name"))] <- paste(names(results)[-which(names(results) %in% c("coord.intron", "gene_short_name"))], ".g2", sep="")
            
            # Save as new object
            results.g2 <- results
            #results.g2 <- results[, c("coord.intron", "psi.g2")]
            
        # Merge
        results <- join(results.g1, results.g2, by="coord.intron", type="left")
        
        # Compute permutated delta
        results$delta.perm <- results$psi.g2 - results$psi.g1
        
        # Save into list
        row.names(results) <- results$coord.intron
        results <- results[, "delta.perm", drop=FALSE]
        
        .list.results.perm[[i]] <- results
             
        if(show.progress==TRUE) {
        
            # Track progress
            setTxtProgressBar(pb, i)
            
        }
        
    }
        
    results.perm <- do.call(cbind.data.frame, .list.results.perm)
    
    # Compute pval
        # Report progress
        message("Computing P values...")
        
        # Annotate delta observed
        index.l <- table(results.obs$coord.intron==row.names(results.perm))
        index.true <- length(which(names(index.l)==TRUE))
        index.false <- length(which(names(index.l)==FALSE))
         
        if(index.true==1 & index.false==0) {
        
            results.perm <- cbind.data.frame(results.perm, results.obs[,"delta",drop=FALSE])
                        
        } else {
            
            return(message("Error in consolidating observed and permutated results"))

        }
        
        # Compute indices
        index.delta.perm <- grep("perm", names(results.perm))
        index.delta.obs <- which(names(results.perm)=="delta")

        # Compute pval
        pval <- apply(results.perm[, ], 1, function(x) {
                
                
                    #x <- as.numeric(results.perm[2, ]) # Test
                                
                    # Retrieve delta observed, perm
                    delta.obs <- x[index.delta.obs]
                    delta.perm <- x[index.delta.perm]
                    
                    # Compute pval (1-sided)
                    #if(delta.obs < 0) {
                        
                        #pval <- sum(delta.perm < delta.obs) / n.iterations
                        
                    #} else if(delta.obs >= 0) {
                        
                        #pval <- sum(delta.perm > delta.obs) / n.iterations
                        
                    #}
                    
                    # Compute pval (2-sided)
                    pval <- sum(abs(delta.perm) > abs(delta.obs)) / n.iterations
                        
                    return(pval)
                    
                })
    
        results.obs$pval <- pval
        
        # Check results
        #results.small <- results.obs[which(results.obs$pval < 0.05 & results.obs$delta > 0), ]
        #results.small <- results.obs[which(results.obs$pval < 0.05 & results.obs$delta < 0), ]
        #table(results.obs$pval < 0.05)
        
    # Re-order by pval
    results.obs <- results.obs[order(results.obs$pval), ]
    
    ################################################################

    # Compute mean norm gene expression
    #. <- apply(df.gene.norm, 1, function(x) {mean(log2(x + 1))})
    #. <- data.frame("gene_short_name"=names(.),
                    #"mean.expr.gene.norm.g1.g2"=as.numeric(.),
                    #stringsAsFactors=FALSE
                    #)
    
    # Annotate
    results.obs <- join(results.obs, mean.combined.df, by="gene_short_name", type="left")
    
    ################################################################
 
    # Save into new slot
    MarvelObject$DE$SJ$Table <- results.obs
    MarvelObject$DE$SJ$cell.group.g1 <- cell.group.g1
    MarvelObject$DE$SJ$cell.group.g2 <- cell.group.g2
    
    # Return final object
    return(MarvelObject)
            
}


