#' @title Tabulate modality proportion
#'
#' @description Tabulates and plots the proportion of each modality. This is a wrapper function for \code{PropModality.Doughnut} and \code{PropModality.Bar} functions.
#'
#' @param MarvelObject Marvel object. S3 object generated from \code{AssignModality} function.
#' @param modality.column Character string. Can take the value \code{"modality"}, \code{"modality.var"} or \code{"modality.bimodal.adj"}. Please refer to \code{AssignModality} function help page for more details.
#' @param modality.type Character string. \code{basic} indicates that only the main modalities (included, excluded, bimodal, middle, multimodal) are analysed. Sub-modalities (primary and dispersed) will be merged. \code{complete} indicates that both main and sub-modalities are analysed. Sub-modalities will not be merged.
#' @param event.type Character string. To indicate which event type to analyse. Can take the value \code{"SE"}, \code{"MXE"}, \code{"RI"}, \code{"A5SS"} or \code{"A3SS"}. Specify \code{"all"} to include all event types.
#' @param across.event.type Logical. If set to \code{TRUE}, the proportion of modality will be compared across the specified event types
#' @param prop.test Character string. Only applicable when \code{across.event.type} set to \code{TRUE}. \code{chisq} Chi-squared test used to compare the proportion of modalities across the different event splicing type. \code{fisher} Fisher test used to compare the proportion of modalities across the different splicing event type.
#' @param prop.adj Character string. Only applicable when \code{across.event.type} set to \code{TRUE}. Adjust p-values generated from \code{prop.test} for multiple testing. Options available as per \code{p.adjust} function.
#' @param xlabels.size Numeric value. Only applicable when \code{across.event.type} set to \code{TRUE}. Size of x-axis labels as per \code{ggplot2} function. Default is 8.
#' @param zoom Logical value. Only applicable if \code{across.event.type} set to \code{TRUE}. If set to \code{TRUE}, users can specify the range of the y-axis using \code{yinterval} argument. Useful when scrutinasing low-frequency event types, e.g. middle and multimodal.
#' @param yinterval Logical value. Only applicable if \code{across.event.type} set to \code{TRUE} and \code{zoom} set to \code{TRUE}.
#'
#' @return An object of class S3 containing with new slot \code{$Modality$Prop$DoughnutChart} or \code{$Modality$Prop$BarChart}.
#'
#' @importFrom plyr join
#' @importFrom stats p.adjust p.adjust.methods
#' @import methods
#' @import ggplot2
#'
#' @export
#'
#' @examples
#' marvel.demo <- readRDS(system.file("extdata/data", "marvel.demo.rds", package="MARVEL"))
#'
#' marvel.demo <- PropModality(MarvelObject=marvel.demo,
#'                             modality.column="modality.bimodal.adj",
#'                             modality.type="extended",
#'                             event.type=c("SE", "MXE", "RI", "A5SS", "A3SS", "AFE", "ALE"),
#'                             across.event.type=FALSE
#'                             )
#'
#' # Check outputs
#' marvel.demo$Modality$Prop$DoughnutChart$Table
#' marvel.demo$Modality$Prop$DoughnutChart$Plot

PropModality <- function(MarvelObject, modality.column, modality.type, event.type, across.event.type, prop.test=NULL, prop.adj=NULL, xlabels.size=8, zoom=FALSE, yinterval=NULL) {
    
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
                              prop.adj=prop.adj,
                              xlabels.size=xlabels.size,
                              zoom=zoom,
                              yinterval=yinterval
                              )
        
    }
            
}
