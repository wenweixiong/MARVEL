#' @title Check Alignment for both splicing and gene data
#'
#' @description \code{CheckAlignment} subset overlapping samples between splicing and gene data.
#'
#' @details This function subset overlapping samples between splicing and gene data and aligns the order of sampleIDs between the two data types.
#'
#' @param MarvelObject S3 object generated from \code{CheckAlignment} function.
#'
#' @export
#'
#' @return An object of class S3. The original \code{MarvelObject$SplicePheno}, \code{MarvelObject$PSI}, \code{MarvelObject$GenePheno}, and \code{MarvelObject$Exp} are updated.
#'
#' @author Sean Wen <sean.wenwx@gmail.com>
#'
#' @import methods
#'
#' @examples
#' # Load input
#' marvel <- readRDS(system.file("extdata/Data", "MarvelObject.rds", package="MARVEL"))
#'
#' # Run example
#' marvel <- CheckAlignment.PSI.Exp(MarvelObject=marvel)

CheckAlignment.PSI.Exp <- function(MarvelObject) {
        
    # Define arguments
    MarvelObject <- MarvelObject
    
    # Example argument
    #MarvelObject <- marvel
    
    # Retrieve non-empty matrix event types
    ncol <- c(ncol(MarvelObject$PSI[["SE"]]),
              ncol(MarvelObject$PSI[["MXE"]]),
              ncol(MarvelObject$PSI[["RI"]]),
              ncol(MarvelObject$PSI[["A5SS"]]),
              ncol(MarvelObject$PSI[["A3SS"]])
              )
    event.types <- c("SE", "MXE", "RI", "A5SS", "A3SS")
    event.types <- event.types[which(ncol >= 2)]
    
    #########################################################################
    ############################### PHENODATA ###############################
    #########################################################################
    
    # Retrieve sample IDs from splicing phenoData
    df.pheno.psi <- MarvelObject$SplicePheno
    sample.ids.psi <- df.pheno.psi$sample.id
    
    print(paste(length(sample.ids.psi), " samples identified in splicing phenoData", sep=""))
    
    # Retrieve sample IDs from splicing phenoData
    df.pheno.exp <- MarvelObject$GenePheno
    sample.ids.exp  <- df.pheno.exp $sample.id
    
    print(paste(length(sample.ids.exp), " samples identified in gene phenoData", sep=""))
    
    # Find overlap
    overlap.psi.exp <- intersect(sample.ids.psi, sample.ids.exp)
    print(paste(length(overlap.psi.exp), " overlapping samples identified and retained", sep=""))
    
    # Subset overlaps
    df.pheno.psi <- df.pheno.psi[which(df.pheno.psi$sample.id %in% overlap.psi.exp),]
    df.pheno.exp <- df.pheno.exp[which(df.pheno.exp$sample.id %in% overlap.psi.exp),]
    
    # Ensure both splicing and gene phenoData has same row order
    row.names(df.pheno.psi) <- df.pheno.psi$sample.id
    row.names(df.pheno.exp) <- df.pheno.exp$sample.id
    
    df.pheno.psi <- df.pheno.psi[overlap.psi.exp, ]
    df.pheno.exp <- df.pheno.exp[overlap.psi.exp, ]
    
    row.names(df.pheno.psi) <- NULL
    row.names(df.pheno.exp) <- NULL
    
    # Check alignment
    index.l <- table(df.pheno.psi$sample.id==df.pheno.exp$sample.id)
    index.true <- length(which(names(index.l)==TRUE))
    index.false <- length(which(names(index.l)==FALSE))
     
    if(index.true==1 & index.false==0) {
    
        print("splicing and gene phenoData sample IDs MATCHED")
        
    } else {
        
        
        print("splicing and gene phenoData sample IDs NOT MATCHED")

    }

    # Save as into MARVEL slot
    MarvelObject$SplicePheno <- df.pheno.psi
    MarvelObject$GenePheno <- df.pheno.exp
    
    # Retrieve overlapping samples IDs: PSI matrix
        # Retrieve overlapping IDs
        .list <- list()
        
        for(i in 1:length(event.types)) {
             
            .list[[i]] <- names(MarvelObject$PSI[[event.types[i]]])
            
            
        }
        
        sample.ids.matrix <- unique(unlist(.list))
        
        overlap <- intersect(overlap.psi.exp, sample.ids.matrix)
        
        # Subset and return MARVEL object
        for(i in 1:length(event.types)) {
             
             df <- MarvelObject$PSI[[event.types[i]]]
            
             df.small <- df[, c("tran_id", overlap)]
             
             MarvelObject$PSI[[event.types[i]]] <- df.small
            
        }
        
    # Retrieve overlapping samples IDs: Expression matrix
        # Retrieve overlapping IDs
        df.exp <- MarvelObject$Exp
        
        overlap <- intersect(overlap.psi.exp, names(df.exp)[-1])
        
        # Subset and return MARVEL object
        df.exp <- df.exp[,c("gene_id", overlap)]
        MarvelObject$Exp <- df.exp
        
    # Check alignment
    for(i in 1:length(event.types)) {
        
        index.l <- table(names(MarvelObject$PSI[[event.types[i]]])[-1]==names(MarvelObject$Exp)[-1])
        index.true <- length(which(names(index.l)==TRUE))
        index.false <- length(which(names(index.l)==FALSE))
         
        if(index.true==1 & index.false==0) {
        
            print(paste(event.types[i], " splicing and gene matrix column names MATCHED", sep=""))
            
        } else {
            
            
            print(paste(event.types[i], " splicing and gene matrix column names MATCHED", sep=""))
            
        }
        
    }

    #########################################################################
    
    # Return final object
    return(MarvelObject)
    
}


