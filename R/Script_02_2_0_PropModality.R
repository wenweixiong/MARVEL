#' @title Modality Proportion
#'
#' @description
#' \code{PropModality} tabulates and plots the proportion of each modality.
#'
#' @details
#' This function tabulates and plots the proportion of each modality for a specified splicing event type(s) or compares proportion of each modality across specified splicing event types.
#'
#' @param MarvelObject S3 object generated from \code{AssignModality} function.
#' @param modality.column Character string. Can take the value \code{"modality"}, \code{"modality.var"} or \code{"modality.bimodal.adj"}. Please refer to \code{AssignModality} function help page for more details.
#' @param modality.type Character string. \code{basic} indicates that only the main modalities (included, excluded, bimodal, middle, multimodal) are analysed. Sub-modalities (primary and dispersed) will be merged. \code{extended} indicates that both main and sub-modalities are analysed. Sub-modalities will not be merged.
#' @param event.type Character string. To indicate which event type to analyse. Can take the value \code{"SE"}, \code{"MXE"}, \code{"RI"}, \code{"A5SS"} or \code{"A3SS"}. Specify \code{"all"} to include all event types.
#' @param across.event.type Logical. If set to \code{TRUE}, the proportion of modality will be compared across the specified event types
#' @param prop.test Character string. Only applicable when \code{across.event.type} set to \code{TRUE}. \code{chisq} Chi-squared test used to compare the proportion of modalities across the different event splicing type. \code{fisher} Fisher test used to compare the proportion of modalities across the different splicing event type.
#' @param prop.adj Character string. Only applicable when \code{across.event.type} set to \code{TRUE}. Adjust p-values generated from \code{prop.test} for multiple testing. Options available as per \code{p.adjust} function.
#' @export
#' @return An object of class S3 containing all the original slots as inputted by the user in addition to one new slot. If \code{across.event.type} set to \code{FALSE}, results returned to \code{$Modality$Prop$DoughnutChart} slot. If \code{across.event.type} set to \code{TRUE}, results returned to \code{$Modality$Prop$BarChart} slot.
#' @author Sean Wen <sean.wenwx@gmail.com>
#' @importFrom fitdistrplus fitdist
#' @import methods
#' @import stats
#' @import ggplot2
#' @importFrom plyr join

PropModality <- function(MarvelObject, modality.column, modality.type, event.type, across.event.type, prop.test=NULL, prop.adj=NULL) {
    
    if(across.event.type==FALSE) {
        
        PropModality.Doughnut(MarvelObject=MarvelObject,
                              modality.column=modality.column,
                              modality.type=modality.type,
                              event.type=event.type
                              )
        
    } else if(across.event.type==TRUE) {
        
        PropModality.Bar(MarvelObject=MarvelObject,
                              modality.column=modality.column,
                              modality.type=modality.type,
                              event.type=event.type,
                              prop.test=prop.test,
                              prop.adj=prop.adj
                              )
        
    }
            
}
