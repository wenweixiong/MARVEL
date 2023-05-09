#' @title Preprocess rMATS files for MARVEL
#'
#' @description Preprocess rMATS coordinates for MARVEL input and annotate the gene type for each gene using the GTF provided.
#'
#' @param file Data frame. Tab-delimited file output from rMATS. Typically the file is named \code{fromGTF.*.txt}.
#' @param GTF Data frame. The same GTF file used in the rMATS step.
#' @param event_type. Indicate the type of splicing event. Options are \code{"SE"}, \code{"MXE"}, \code{"RI"}, \code{"A5SS"}, or \code{"A3SS"}.
#'
#' @return An object of class data frame that may be used as input for MARVEL
#'
#' @importFrom plyr join
#' @import methods
#'
#' @export

Proproces_rMATS <- function(file, GTF, EventType) {
    
    ###############################################################
    
    if(EventType=="SE"){
        
        Proproces_rMATS.SE(file=file, GTF=GTF)
        
    } else if(EventType=="MXE"){
        
        Proproces_rMATS.MXE(file=file, GTF=GTF)
        
    } else if(EventType=="RI"){
        
        Proproces_rMATS.RI(file=file, GTF=GTF)
        
    } else if(EventType=="A5SS"){
        
        Proproces_rMATS.A5SS(file=file, GTF=GTF)
        
    } else if(EventType=="A3SS"){
        
        Proproces_rMATS.A3SS(file=file, GTF=GTF)
        
    }
    
    ###############################################################
   
}
