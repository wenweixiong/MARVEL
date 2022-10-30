#' @title Create Marvel object for plate-based RNA-sequencing data
#'
#' @description Creates an S3 object named \code{Marvel} for downstream analysis, specifically for plate-based RNA-sequencing data.
#'
#' @param SplicePheno Data frame. Sample metadata.
#' @param SpliceJunction Data frame. Splice junction counts matrix.
#' @param IntronCounts Data frame. Intron coverage matrix.
#' @param SpliceFeature List of data frames. Each data frame is the exon-level alternative splicing event metadata.
#' @param SpliceFeatureValidated List of data frames. Each data frame is the validated (high-quality) exon-level alternative splicing event metadata.
#' @param PSI Data frame. PSI matrix.
#' @param GeneFeature Data frame. Gene metadata.
#' @param Exp Data frame. Normalised, non-log2-transformed gene expression matrix.
#' @param GTF Data frame. GTF used for generating the exon-level alternative splicing event metadata.
#'
#' @return An object of class S3.
#'
#' @importFrom plyr join
#' @import methods
#'
#' @export
#'
#' @examples
#' marvel.demo <- readRDS(system.file("extdata/data", "marvel.demo.rds", package="MARVEL"))
#'
# Retrieve, observe format of pre-saved input files
#' SpliceJunction <- marvel.demo$SpliceJunction
#' SpliceJunction[1:5,1:5]
#'
#' SplicePheno <- marvel.demo$SplicePheno
#' SplicePheno[1:5,]
#'
#' SpliceFeature <- marvel.demo$SpliceFeature
#' SpliceFeature[["SE"]][1:5, ]
#'
#' IntronCounts <- marvel.demo$IntronCounts
#' IntronCounts[1:5,1:5]
#'
#' GeneFeature <- marvel.demo$GeneFeature
#' GeneFeature[1:5, ]
#'
#' Exp <- marvel.demo$Exp
#' Exp[1:5,1:5]
#'
#' marvel <- CreateMarvelObject(SpliceJunction=SpliceJunction,
#'                              SplicePheno=SplicePheno,
#'                              SpliceFeature=SpliceFeature,
#'                              IntronCounts=IntronCounts,
#'                              GeneFeature=GeneFeature,
#'                              Exp=Exp
#'                              )
#' class(marvel)

CreateMarvelObject <- function(SplicePheno=NULL,
                               SpliceJunction=NULL,
                               IntronCounts=NULL,
                               SpliceFeature=NULL,
                               SpliceFeatureValidated=NULL,
                               PSI=NULL,
                               GeneFeature=NULL,
                               Exp=NULL,
                               GTF=NULL
                               ) {
        
    # Create s3 object
    s3 <- list()
    class(s3) <- "Marvel"
    
    # Fill slots
    s3$SplicePheno <- SplicePheno
    s3$SpliceJunction <- SpliceJunction
    s3$IntronCounts <- IntronCounts
    
    s3$SpliceFeature$SE <- SpliceFeature[["SE"]]
    s3$SpliceFeature$MXE <- SpliceFeature[["MXE"]]
    s3$SpliceFeature$RI <- SpliceFeature[["RI"]]
    s3$SpliceFeature$A5SS <- SpliceFeature[["A5SS"]]
    s3$SpliceFeature$A3SS <- SpliceFeature[["A3SS"]]
    s3$SpliceFeature$ALE <- SpliceFeature[["ALE"]]
    s3$SpliceFeature$AFE <- SpliceFeature[["AFE"]]
    
    s3$SpliceFeatureValidated$SE <- SpliceFeatureValidated[["SE"]]
    s3$SpliceFeatureValidated$MXE <- SpliceFeatureValidated[["MXE"]]
    s3$SpliceFeatureValidated$RI <- SpliceFeatureValidated[["RI"]]
    s3$SpliceFeatureValidated$A5SS <- SpliceFeatureValidated[["A5SS"]]
    s3$SpliceFeatureValidated$A3SS <- SpliceFeatureValidated[["A3SS"]]
    s3$SpliceFeatureValidated$ALE <- SpliceFeatureValidated[["ALE"]]
    s3$SpliceFeatureValidated$AFE <- SpliceFeatureValidated[["AFE"]]
    
    s3$PSI$SE <- PSI[["SE"]]
    s3$PSI$MXE <- PSI[["MXE"]]
    s3$PSI$RI <- PSI[["RI"]]
    s3$PSI$A5SS <- PSI[["A5SS"]]
    s3$PSI$A3SS <- PSI[["A3SS"]]
    s3$PSI$ALE <- PSI[["ALE"]]
    s3$PSI$AFE <- PSI[["AFE"]]
    
    #s3$GenePheno <- GenePheno
    s3$GeneFeature <- GeneFeature
    s3$Exp <- Exp
    
    s3$GTF <- GTF
    
    # Returen final object
    return(s3)
            
}


