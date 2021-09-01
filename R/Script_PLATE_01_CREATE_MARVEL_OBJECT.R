#' @title Create Marvel Object
#'
#' @description \code{CreateMarvelObject} creates an S3 object named \code{Marvel} for downstream analysis.
#'
#' @details This function creates an S3 object named \code{Marvel} for downstream analysis. It can take both splicing and gene expression data. Gene expression data is highly encouraged so that users can compare and contrast splicing and gene expression profiles using other functionalities by \code{MARVEL}.
#'
#' @param SplicePheno Data frame containing sample metadata. This object should consist of at least 2 columns. Mandatory columns are \code{sample.id} and \code{cell.type}. \code{sample.id} indicates all unique sample IDs. \code{cell.type} indicates the group names of each sample. Additional columns, if present, may contain additional details of each sample such sequencing QC details etc..
#' @param SpliceJunction Data frame containing splice junction counts. First column should be named \code{coord.intron} and indicate the splice junction position in the form of chr:start:end. Subsequent columns should contain the splice junction counts for each sample. These junctions can be detected using external softwares such as STAR, featureCounts etc..
#' @param IntronCounts Data frame. Columns indicate sample IDs, rows indicate intron coordinates, and values indicate total intron coverage. The first column needs to be named \code{coord.intron}. These values will be combined with splice junction counts in the MARVEL object to compute PSI values.
#' @param SpliceFeature List containing splicing event metadata. Each element in the list is a data frame for each splicing event type, i.e. SE, MXE, RI, A5SS, and A3SS. Names of each element should reflect the splicing event type. Mandatory columns in each data frame are \code{tran_id} and \code{gene_id}. \code{tran_id} indicates all the splicing event coordinates. These events can be detected using external softwares such as rMATS, MISO, BRIE etc.. Other columns, if present, may contain additional details of each splicing event such as gene name, gene type etc..
#' @param SpliceFeatureValidated List containing validated splicing event metadata. Required when \code{SpliceFeature} and \code{SpliceJunction} not specified. Each element in the list is a data frame for each splicing event type, i.e. SE, MXE, RI, A5SS, and A3SS. Names of each element should reflect the splicing event type. Mandatory columns in each data frame are \code{tran_id} and \code{gene_id}. \code{tran_id} indicates all the splicing event coordinates. These events can be detected using external softwares such as rMATS, MISO, BRIE etc.. Additional columns, if present, may contain additional details of each splicing event such as gene name, gene type etc..
#' @param PSI Data frame containing pre-computed PSI values. Required when \code{SpliceFeature} and \code{SpliceJunction} not specified. The first column should be named \code{tran_id} and second column onwards should be sample names containing the PSI values.
#' @param GenePheno Data frame containing sample metadata. Optional but highly recommended. This object should consist of at least 2 columns. Mandatory columns are \code{sample.id} and \code{cell.type}. \code{sample.id} indicates all unique sample IDs. \code{cell.type} indicates the group names of each sample. Additional columns, if present, may contain additional details of each sample such sequencing QC details etc..
#' @param GeneFeature Data frame containing the gene metadata. Optional but highly recommended. Mandatory column is \code{gene_id}. Additional columns, if present, may contain details of each gene, e.g. gene name, gene type etc..
#' @param Exp Data frame containing normalised and log-transformed gene expression values. Optional but highly recommended. The first column should be named \code{gene_id} and second column onwards should be sample names containing the gene expression values.
#'
#' @export
#'
#' @return An object of class S3. Each slot in the object is named after the \code{CreateMarvelObject} arguments.
#'
#' @author Sean Wen <sean.wenwx@gmail.com>
#'
#' @import methods
#' 
#' @examples
#' # Read sample metadata file
#' path_to_file <- system.file("extdata/Data", "SJ_phenoData.txt", package="MARVEL")
#' df.pheno <- read.table(path_to_file, sep="\t", header=TRUE, stringsAsFactors=FALSE,
#'              na.strings="NA")
#'
#' # Read splice junction file
#' path_to_file <- system.file("extdata/Data", "SJ.txt", package="MARVEL")
#' sj <- read.table(paste(path_to_file), sep="\t", header=TRUE, stringsAsFactors=FALSE,
#'                   na.strings="NA")
#'
#' # Read splicing event metadata file
#' df.feature.list <- list()
#'
#' path_to_file <- system.file("extdata/Data", "SE_featureData.txt", package="MARVEL")
#' df.feature.list[[1]] <- read.table(paste(path_to_file), sep="\t", header=TRUE,
#'                          stringsAsFactors=FALSE, na.strings="NA")
#'
#' names(df.feature.list) <- "SE"
#'
#' # Create MARVEL object
#' marvel <- CreateMarvelObject(SplicePheno=df.pheno,
#'                               SpliceJunction=sj,
#'                               SpliceFeature=df.feature.list
#'                               )
#'
#' class(marvel)


CreateMarvelObject <- function(SplicePheno=NULL,
                              SpliceJunction=NULL,
                              IntronCounts=NULL,
                              SpliceFeature=NULL,
                              SpliceFeatureValidated=NULL,
                              PSI=NULL,
                              GenePheno=NULL,
                              GeneFeature=NULL,
                              Exp=NULL) {
        
        
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
    
    s3$SpliceFeatureValidated$SE <- SpliceFeatureValidated[["SE"]]
    s3$SpliceFeatureValidated$MXE <- SpliceFeatureValidated[["MXE"]]
    s3$SpliceFeatureValidated$RI <- SpliceFeatureValidated[["RI"]]
    s3$SpliceFeatureValidated$A5SS <- SpliceFeatureValidated[["A5SS"]]
    s3$SpliceFeatureValidated$A3SS <- SpliceFeatureValidated[["A3SS"]]
    
    s3$PSI$SE <- PSI[["SE"]]
    s3$PSI$MXE <- PSI[["MXE"]]
    s3$PSI$RI <- PSI[["RI"]]
    s3$PSI$A5SS <- PSI[["A5SS"]]
    s3$PSI$A3SS <- PSI[["A3SS"]]

    s3$GenePheno <- GenePheno
    s3$GeneFeature <- GeneFeature
    s3$Exp <- Exp
    
    # Returen final object
    return(s3)
            
}


