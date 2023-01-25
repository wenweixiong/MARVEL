#' @title Differential splice junction analysis at donor level
#'
#' @description Performs differential splice junction analysis between two groups of cells by aggregating cells by their corresponding donor IDs. Statistical test used here is the permutation test.
#'
#' @param MarvelObject Marvel object. S3 object generated from \code{CheckAlignment.10x} function.
#' @param MarvelObject cell.group.list List of character strings. List defining the two cell groups and their respective donor IDs and their corresponding cells IDs.
#' @param min.pct.cells.gene.expr Numeric value. Mininum percentage of cells across a donor expressing a given gene for that gene to be considered expressed. Genes below this threshold will be considered not expressed. Default is \code{10}. To be used in conjunction with \code{min.n.cells.gene.expr} and \code{min.gene.counts.total}.
#' @param min.n.cells.gene.expr Numeric value. Mininum number of cells across a donor expressing a given gene for that gene to be considered expressed. Genes below this threshold will be considered not expressed. Default is \code{10}. To be used in conjunction with \code{min.pct.cells.gene.expr} and \code{min.gene.counts.total}.
#' @param min.gene.counts.total Numeric value. Mininum number of UMI counts across a donor for that gene to be considered expressed. Genes below this threshold will be considered not expressed. Default is \code{3}. To be used in conjunction with \code{min.n.cells.gene.expr} and \code{min.pct.cells.gene.expr}.
#' @param min.sj.count Numeric value. Mininum number of UMI counts across a donor for that splice junction to be considered expressed. Splice junctions below this threshold will be considered not expressed. Default is \code{1}.
#' @param min.donor.gene.expr Numeric value. Minimum number of donors expressing a given gene based on \code{min.pct.cells.gene.expr}, \code{min.n.cells.gene.expr}, and \code{min.gene.counts.total} options for that gene to be included for analysis. Default is \code{3}.
#' @param min.donor.sj.expr Numeric value. Minimum number of donors expressing a given splice junction based on \code{min.sj.count} option for that splice junction to be included for analysis. Default is \code{3}.
#'
#' @return An object of class S3 with a new slots \code{MarvelObject$DE$Pseudobulk$SJ$Table} and \code{MarvelObject$DE$Pseudobulk$SJ$cell.group.list}.
#'
#' @importFrom plyr join
#' @import Matrix
#' @importFrom utils txtProgressBar setTxtProgressBar
#' 
#' @export

CompareValues.SJ.DonorLevel.10x <- function(MarvelObject, cell.group.list, min.pct.cells.gene.expr=10, min.n.cells.gene.expr=10, min.gene.counts.total=3, min.sj.count=1, min.donor.gene.expr=3, min.donor.sj.expr=3) {
        
    # Define arguments
    MarvelObject <- MarvelObject
    sample.metadata <- MarvelObject$sample.metadata
    sj.metadata <- MarvelObject$sj.metadata
    df.gene.norm <- MarvelObject$gene.norm.matrix
    df.gene.count <- MarvelObject$gene.count.matrix
    df.sj.count <- MarvelObject$sj.count.matrix
    cell.group.list <- cell.group.list
    min.pct.cells.gene.expr <- min.pct.cells.gene.expr
    min.n.cells.gene.expr <- min.n.cells.gene.expr
    min.gene.counts.total <- min.gene.counts.total
    min.sj.count <- min.sj.count
    min.donor.gene.expr <- min.donor.gene.expr
    min.donor.sj.expr <- min.donor.sj.expr

    # Example arguments
    #MarvelObject <- marvel
    #sample.metadata <- MarvelObject$sample.metadata
    #sj.metadata <- MarvelObject$sj.metadata
    #df.gene.norm <- MarvelObject$gene.norm.matrix
    #df.gene.count <- MarvelObject$gene.count.matrix
    #df.sj.count <- MarvelObject$sj.count.matrix
    #cell.group.list <- cell.group.list
    #min.pct.cells.gene.expr <- 10
    #min.n.cells.gene.expr <- 10
    #min.gene.counts.total <- 3
    #min.sj.count <- 1
    #min.donor.gene.expr <- 5
    #min.donor.sj.expr <- 5

    #################################################################
    ################# CREATE SAMPLE METADATA ########################
    #################################################################
    
    .list.2 <- list()
    
    for(j in 1:length(cell.group.list)) {
        
        cell.group.list.small <- cell.group.list[[j]]
        
        .list <- list()
        
        for(i in 1:length(cell.group.list.small)) {
            
            .list[[i]] <- data.frame("cell.id"=unlist(cell.group.list.small[[i]]),
                                     "cell.group"=names(cell.group.list)[j],
                                     "donor.id"=names(cell.group.list.small)[i]
                                     )
                                            
        }
        
        .list.2[[j]] <- do.call(rbind.data.frame, .list)
        
    }
    
    sample.metadata <- do.call(rbind.data.frame, .list.2)
        
    #################################################################
    ############### SUBSET EXPRESSED GENES: GROUP 1 #################
    #################################################################
    
    # Define donor.ids
    index <- which(sample.metadata$cell.group==names(cell.group.list)[1])
    donor.ids <- unique(sample.metadata[index, "donor.id"])
    
    # Tabulate metrices
    .list <- list()
    
    message("Retrieving expressed genes for cell group 1...")
    
    pb <- txtProgressBar(1, length(donor.ids), style=3)
    
    for(i in 1:length(donor.ids)) {
        
        # Define donor.id
        donor.id <- donor.ids[i]
        
        # Definie cells
        index <- which(sample.metadata$donor.id==donor.id)
        cell.ids <- sample.metadata[index, "cell.id"]
                                
        # Subset matrix
        df.gene.count.small <- df.gene.count[, cell.ids]
        
        # Compute metrices: Cells
        . <- apply(df.gene.count.small, 1, function(x) { sum(x != 0)})
        results <- data.frame("cell.group"=names(cell.group.list)[1],
                              "donor.id"=donor.id,
                              "gene_short_name"=names(.),
                              "n.cells.total"=length(cell.ids),
                              "n.cells.gene.expr"=as.numeric(.),
                              "pct.cells.gene.expr"=round(as.numeric(.)/length(cell.ids) * 100, digits=2),
                              stringsAsFactors=FALSE
                              )
                                    
        # Compute metrices: Counts
        . <- apply(df.gene.count.small, 1, function(x) { sum(x)})
        . <- data.frame("gene_short_name"=names(.),
                        "gene.counts.total"=.,
                        stringsAsFactors=FALSE
                        )

        results <- join(results, ., by="gene_short_name", type="left")
            
        # Save into list
        .list[[i]] <- results
        
        # Track progress
        setTxtProgressBar(pb, i)
                
    }
    
    results <- do.call(rbind.data.frame, .list)
    
    # Subset expressed genes
    index <- which(results$pct.cells.gene.expr >= min.pct.cells.gene.expr &
                   results$n.cells.gene.expr >= min.n.cells.gene.expr &
                   results$gene.counts.total >= min.gene.counts.total
                   )
    
    results.small <- results[index, ]
    
    # Tabulate n donor expressing genes
    . <- as.data.frame(table(results.small$gene_short_name))
    names(.) <- c("gene_short_name", "freq")
    
    # Subest expressed genes
    index <- which(.$freq >= min.donor.gene.expr)
    gene_short_names <- .[index, "gene_short_name"]
    gene_short_names <- as.character(gene_short_names)
   
    # Report progress
    message(paste(length(gene_short_names), " genes expressed in cell group 1", sep=""))
    
    # Save as new object
    gene_short_names.1 <- gene_short_names
    
    #################################################################
    ############### SUBSET EXPRESSED GENES: GROUP 2 #################
    #################################################################
    
    # Define donor.ids
    index <- which(sample.metadata$cell.group==names(cell.group.list)[2])
    donor.ids <- unique(sample.metadata[index, "donor.id"])
    
    # Tabulate metrices
    .list <- list()
    
    message("Retrieving expressed genes for cell group 2...")
    
    pb <- txtProgressBar(1, length(donor.ids), style=3)
    
    for(i in 1:length(donor.ids)) {
        
        # Define donor.id
        donor.id <- donor.ids[i]
        
        # Definie cells
        index <- which(sample.metadata$donor.id==donor.id)
        cell.ids <- sample.metadata[index, "cell.id"]
                                
        # Subset matrix
        df.gene.count.small <- df.gene.count[, cell.ids]
        
        # Compute metrices: Cells
        . <- apply(df.gene.count.small, 1, function(x) { sum(x != 0)})
        results <- data.frame("cell.group"=names(cell.group.list)[2],
                              "donor.id"=donor.id,
                              "gene_short_name"=names(.),
                              "n.cells.total"=length(cell.ids),
                              "n.cells.gene.expr"=as.numeric(.),
                              "pct.cells.gene.expr"=round(as.numeric(.)/length(cell.ids) * 100, digits=2),
                              stringsAsFactors=FALSE
                              )
                                    
        # Compute metrices: Counts
        . <- apply(df.gene.count.small, 1, function(x) { sum(x)})
        . <- data.frame("gene_short_name"=names(.),
                        "gene.counts.total"=.,
                        stringsAsFactors=FALSE
                        )

        results <- join(results, ., by="gene_short_name", type="left")
            
        # Save into list
        .list[[i]] <- results
        
        # Track progress
        setTxtProgressBar(pb, i)
                
    }
    
    results <- do.call(rbind.data.frame, .list)
    
    # Subset expressed genes
    index <- which(results$pct.cells.gene.expr >= min.pct.cells.gene.expr &
                   results$n.cells.gene.expr >= min.n.cells.gene.expr &
                   results$gene.counts.total >= min.gene.counts.total
                   )
    
    results.small <- results[index, ]
    
    # Tabulate n donor expressing genes
    . <- as.data.frame(table(results.small$gene_short_name))
    names(.) <- c("gene_short_name", "freq")
    
    # Subest expressed genes
    index <- which(.$freq >= min.donor.gene.expr)
    gene_short_names <- .[index, "gene_short_name"]
    gene_short_names <- as.character(gene_short_names)
   
    # Report progress
    message(paste(length(gene_short_names), " genes expressed in >= ", min.donor.gene.expr, " donors in cell group 2", sep=""))
    
    # Save as new object
    gene_short_names.2 <- gene_short_names
    
    #################################################################
    ############ SUBSET EXPRESSED GENES: GROUP 1 & 2 ################
    #################################################################
    
    # Find overlaps
    gene_short_names <- intersect(gene_short_names.1, gene_short_names.2)
    
    # Track progress
    message(paste(length(gene_short_names), " genes expressed in BOTH cell groups", sep=""))
    
    #################################################################
    ################ SUBSET EXPRESSED SJ: GROUP 1 ###################
    #################################################################
    
    # Define donor.ids
    index <- which(sample.metadata$cell.group==names(cell.group.list)[1])
    donor.ids <- unique(sample.metadata[index, "donor.id"])
    
    # Tabulate metrices
    .list <- list()
    
    message("Retrieving expressed SJs for cell group 1...")
    
    pb <- txtProgressBar(1, length(donor.ids), style=3)

    for(i in 1:length(donor.ids)) {
        
        # Define donor.id
        donor.id <- donor.ids[i]
        
        # Definie cells
        index <- which(sample.metadata$donor.id==donor.id)
        cell.ids <- sample.metadata[index, "cell.id"]
        
        # Subset cells
        df.sj.count.small <- df.sj.count[, cell.ids]
        
        # Subset expressed genes
        coord.introns <- sj.metadata[which(sj.metadata$gene_short_name.start %in% gene_short_names), "coord.intron"]
        length(coord.introns)
        df.sj.count.small <- df.sj.count.small[coord.introns, ]
        
        # Compute n cells express
        . <- apply(df.sj.count.small, 1, function(x) { sum(x) })
        results <- data.frame("cell.group"=names(cell.group.list)[1],
                              "donor.id"=donor.id,
                              "coord.intron"=names(.),
                              "sj.count"=as.numeric(.),
                              stringsAsFactors=FALSE
                              )
        
        # Save into list
        .list[[i]] <- results
        
        # Track progress
        setTxtProgressBar(pb, i)
                
    }

    results <- do.call(rbind.data.frame, .list)
    
    # Subset expressed SJs
    index <- which(results$sj.count >= min.sj.count)
    results.small <- results[index, ]
    
    # Tabulate n donor expressing SJs
    . <- as.data.frame(table(results.small$coord.intron))
    names(.) <- c("coord.intron", "freq")
    
    # Subset expressed SJs
    index <- which(.$freq >= min.donor.sj.expr)
    coord.introns <- .[index, "coord.intron"]
    coord.introns <- as.character(coord.introns)
   
    # Report progress
    message(paste(length(coord.introns), " SJs expressed in cell group 1", sep=""))
    
    # Save as new object
    coord.introns.1 <- coord.introns
    
    #################################################################
    ################ SUBSET EXPRESSED SJ: GROUP 2 ###################
    #################################################################
    
    # Define donor.ids
    index <- which(sample.metadata$cell.group==names(cell.group.list)[2])
    donor.ids <- unique(sample.metadata[index, "donor.id"])
    
    # Tabulate metrices
    .list <- list()
    
    message("Retrieving expressed SJs for cell group 1...")
    
    pb <- txtProgressBar(1, length(donor.ids), style=3)

    for(i in 1:length(donor.ids)) {
        
        # Define donor.id
        donor.id <- donor.ids[i]
        
        # Definie cells
        index <- which(sample.metadata$donor.id==donor.id)
        cell.ids <- sample.metadata[index, "cell.id"]
        
        # Subset cells
        df.sj.count.small <- df.sj.count[, cell.ids]
        
        # Subset expressed genes
        coord.introns <- sj.metadata[which(sj.metadata$gene_short_name.start %in% gene_short_names), "coord.intron"]
        length(coord.introns)
        df.sj.count.small <- df.sj.count.small[coord.introns, ]
        
        # Compute n cells express
        . <- apply(df.sj.count.small, 1, function(x) { sum(x) })
        results <- data.frame("cell.group"=names(cell.group.list)[2],
                              "donor.id"=donor.id,
                              "coord.intron"=names(.),
                              "sj.count"=as.numeric(.),
                              stringsAsFactors=FALSE
                              )
        
        # Save into list
        .list[[i]] <- results
        
        # Track progress
        setTxtProgressBar(pb, i)
                
    }

    results <- do.call(rbind.data.frame, .list)
    
    # Subset expressed SJs
    index <- which(results$sj.count >= min.sj.count)
    results.small <- results[index, ]
    
    # Tabulate n donor expressing SJs
    . <- as.data.frame(table(results.small$coord.intron))
    names(.) <- c("coord.intron", "freq")
    
    # Subset expressed SJs
    index <- which(.$freq >= min.donor.sj.expr)
    coord.introns <- .[index, "coord.intron"]
    coord.introns <- as.character(coord.introns)
   
    # Report progress
    message(paste(length(coord.introns), " SJs expressed in cell group 2", sep=""))
    
    # Save as new object
    coord.introns.2 <- coord.introns
        
    #################################################################
    ############# SUBSET EXPRESSED SJ: GROUP 1 & 2 ##################
    #################################################################
        
    # Find overlaps
    coord.introns <- unique(c(coord.introns.1, coord.introns.2))
        
    # Report progress
    message(paste(length(coord.introns), " SJs expressed in EITHER cell groups", sep=""))
    
    # Report final numbers
    n.sj <- length(coord.introns)
    n.genes <- length(unique(sj.metadata[which(sj.metadata$coord.intron %in% coord.introns), "gene_short_name.start"]))
    message(paste("Total of ", n.sj, " SJs from ", n.genes, " genes included for DE analysis", sep=""))
        
    #################################################################
    ################ AGGREGATE COUNTS BY DONOR ID ###################
    #################################################################
    
    # Subset SJs for analysis
        # SJ metadata
        sj.metadata <- sj.metadata[which(sj.metadata$coord.intron %in% coord.introns), ]
        
        # SJ count matrix
        df.sj.count <- df.sj.count[sj.metadata$coord.intron, ]
        
        # Gene count matrix
        df.gene.count <- df.gene.count[unique(sj.metadata$gene_short_name.start), ]
        
    # Create pseudobulk
    donor.ids <- unique(sample.metadata$donor.id)
    
    df.sj.count.merged.list <- list()
    df.gene.count.merged.list <- list()
    
    message("Creating pseudobulks...")
    
    pb <- txtProgressBar(1, length(donor.ids), style=3)
    
    for(i in 1:length(donor.ids)) {
        
        # Define donor id
        donor.id <- donor.ids[i]
        
        # Retrieve cell ids
        index <- which(sample.metadata$donor.id==donor.id)
        cell.ids <- sample.metadata[index, "cell.id"]
        
        # Subset, merge matrix
            # SJ
            . <- df.sj.count[,cell.ids]
            . <- rowSums(.)
            . <- data.frame(.)
            names(.) <- donor.id
            df.sj.count.merged.list[[i]] <- .
            
            # Gene
            . <- df.gene.count[,cell.ids]
            . <- rowSums(.)
            . <- data.frame(.)
            names(.) <- donor.id
            df.gene.count.merged.list[[i]] <- .
            
        # Track progress
        setTxtProgressBar(pb, i)
        
    }
    
    
    df.sj.count.merged <- do.call(cbind.data.frame, df.sj.count.merged.list)
    df.gene.count.merged <- do.call(cbind.data.frame, df.gene.count.merged.list)
    
    # Check donor id alignment
    index.l <- table(names(df.sj.count.merged)==names(df.gene.count.merged))
    index.true <- length(which(names(index.l)==TRUE))
    index.false <- length(which(names(index.l)==FALSE))
     
    if(index.true==1 & index.false==0) {
        
        message("Merged SJ and gene count matrix columns MATCHED")
    
    } else {
        
        message("Merged SJ and gene count matrix columns DO NOT MATCHED")

    }
    
    #################################################################
    ######################### COMPUTE PSI ###########################
    #################################################################
    
    # Report progress
    message("Computing PSI...")
    
    # Define SJs
    coord.introns <- sj.metadata$coord.intron
    
    .list <- list()
    
    pb <- txtProgressBar(1, length(coord.introns), style=3)
    
    for(i in 1:length(coord.introns)) {
        
        # Define SJ
        coord.intron <- coord.introns[i]
        
        # Define gene
        index <- which(sj.metadata$coord.intron==coord.intron)
        gene_short_name <- sj.metadata[index, "gene_short_name.start"]
        
        # Retrieve counts
        df.sj.count.merged.small <- df.sj.count.merged[coord.intron, ]
        df.gene.count.merged.small <- df.gene.count.merged[gene_short_name, ]
        
        # Compute PSI
        df.psi <- round(df.sj.count.merged.small/df.gene.count.merged.small * 100, digits=2)
        
        # Save into list
        .list[[i]] <- df.psi
        
        # Track progress
        setTxtProgressBar(pb, i)
        
    }
        
    df.psi <- do.call(rbind.data.frame, .list)
    
    #################################################################
    ##################### STATISTICAL TEST ##########################
    #################################################################
    
    # Report progress
    message("Performing statistical test...")
    
    # Define cell groups
    donor.ids.1 <- names(cell.group.list[[1]])
    donor.ids.2 <- names(cell.group.list[[2]])
        
    # Define SJs
    coord.introns <- sj.metadata$coord.intron
    
    gene_short_names <- NULL
    mean.psi.g1 <- NULL
    mean.psi.g2 <- NULL
    delta <- NULL
    log2fc <- NULL
    pval <- NULL
    
    pb <- txtProgressBar(1, length(coord.introns), style=3)
    
    for(i in 1:length(coord.introns)) {
        
        # Define SJ
        coord.intron <- coord.introns[i]
        
        # Subset SJ
        . <- df.psi[coord.intron, ]
        . <- as.data.frame(t(.))
        names(.) <- "psi"
        .$donor.id <- row.names(.)
        row.names(.) <- NULL
        
        # Annotate cell group
        . <- join(., unique(sample.metadata[,c("donor.id", "cell.group")]), by="donor.id", type="left")
        
        # Retrieve values
        x <- .[which(.$donor.id %in% donor.ids.1), "psi"]
        y <- .[which(.$donor.id %in% donor.ids.2), "psi"]
        
        # Compute statistics
        mean.psi.g1[i] <- mean(x)
        mean.psi.g2[i] <- mean(y)
        delta[i] <- mean(y) - mean(x)
        log2fc[i] <- log2(mean(y)) - log2(mean(x))
        
        # Permutation test
        if(sd(x) != 0 | sd(y) != 0) {
            
            pval[i] <- coin::pvalue(coin::independence_test(psi ~ as.factor(cell.group), ., alternative="two.sided"))
            
        } else {
            
            pval[i] <- 1.0
            
        }
        
        # Retrieve gene name
        index <- which(sj.metadata$coord.intron==coord.intron)
        gene_short_names[i] <- sj.metadata[index, "gene_short_name.start"]
        
        # Track progress
        setTxtProgressBar(pb, i)
        
    }
        
    results <- data.frame("coord.intron"=coord.introns,
                          "gene_short_name"=gene_short_names,
                          "mean.psi.g1"=mean.psi.g1,
                          "mean.psi.g2"=mean.psi.g2,
                          "log2fc"=log2fc,
                          "delta"=delta,
                          "pval"=pval,
                          stringsAsFactors=FALSE
                          )
    
    results <- results[order(results$pval), ]
        
    ################################################################
 
    # Save into new slot
    MarvelObject$DE$Pseudobulk$SJ$Table <- results
    MarvelObject$DE$Pseudobulk$SJ$cell.group.list <- cell.group.list
    
    # Return final object
    return(MarvelObject)
            
}


