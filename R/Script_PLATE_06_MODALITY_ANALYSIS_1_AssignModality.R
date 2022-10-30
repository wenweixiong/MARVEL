#' @title Assign modalities
#'
#' @description Assigns modalities to each splicing event for a specified group of cells.
#'
#' @param MarvelObject Marvel object. S3 object generated from \code{TransformExpValues} function.
#' @param sample.ids Vector of character strings. Sample IDs that constitute the cell group.
#' @param min.cells Numeric value. The minimum no. of cells expressing the splicing event for the event to be included for modality assignment.
#' @param sigma.sq Numeric value. The variance threshold below which the included/excluded modality will be defined as primary sub-modality, and above which it will be defined as dispersed sub-modality.
#' @param bimodal.adjust Logical. When set to \code{TRUE}, MARVEL will identify false bimodal modalities and reassign them as included/excluded modality.
#' @param bimodal.adjust.fc Numeric value. The ratio between the proportion of cells with >0.75 PSI vs <0.25 PSI (and vice versa) below which the splicing event will be classified as bimodal. Only applicable when \code{bimodal.adjust} set to \code{TRUE}. To be used in conjunction with \code{bimodal.adjust.diff}.
#' @param bimodal.adjust.diff Numeric value. The difference between the percentage of cells with >0.75 PSI vs <0.25 PSI (and vice versa) below which the splicing event will be classified as bimodal. Only applicable when \code{bimodal.adjust} set to \code{TRUE}. To be used in conjunction with \code{bimodal.adjust.fc}.
#' @param seed Numeric value. Ensure the \code{fitdist} function returns the same values for alpha and beta paramters each time this function is executed using the same random number generator.
#' @param tran_ids Character strings. Specific vector of transcript IDs for modality assignment. This will be a subset of all transcripts expressed in sufficient number of cells as defined in \code{min.cells} option.
#'
#' @return An object of class S3 containing with new slot \code{MarvelObject$Modality$Results}.
#'
#' @author Sean Wen <sean.wenwx@gmail.com>
#'
#' @importFrom plyr join
#' @importFrom stats runif
#' @import methods
#'
#' @export
#'
#' @examples
#' marvel.demo <- readRDS(system.file("extdata/data", "marvel.demo.rds", package="MARVEL"))
#'
# Define cell group for analysis
#' df.pheno <- marvel.demo$SplicePheno
#' sample.ids <- df.pheno[which(df.pheno$cell.type=="iPSC"), "sample.id"]
#'
#' # Assign modality
#' marvel.demo <- AssignModality(MarvelObject=marvel.demo,
#'                               sample.ids=sample.ids,
#'                               min.cells=5
#'                               )
#'
#' # Check output
#' head(marvel.demo$Modality$Results)

AssignModality <- function(MarvelObject, sample.ids, min.cells=25, sigma.sq=0.001, bimodal.adjust=TRUE, bimodal.adjust.fc=3, bimodal.adjust.diff=50, seed=1, tran_ids=NULL) {

    # Define arguments
    psi <- do.call(rbind.data.frame, MarvelObject$PSI)
    psi.feature <- do.call(rbind.data.frame, MarvelObject$SpliceFeatureValidated)
    psi.pheno <- MarvelObject$SplicePheno
    sample.ids <- sample.ids
    min.cells <- min.cells
    sigma.sq <- sigma.sq
    bimodal.adjust <- bimodal.adjust
    bimodal.adjust.fc <- bimodal.adjust.fc
    bimodal.adjust.diff <- bimodal.adjust.diff
    seed <- seed
    
    # Example arguments
    #MarvelObject <- marvel.demo
    #psi <- do.call(rbind.data.frame, MarvelObject$PSI)
    #psi.feature <- do.call(rbind.data.frame, MarvelObject$SpliceFeatureValidated)
    #psi.pheno <- MarvelObject$SplicePheno
    #sample.ids <- sample.ids
    #min.cells <- 5
    #sigma.sq <- 0.001
    #bimodal.adjust <- TRUE
    #bimodal.adjust.fc <- 3.0
    #bimodal.adjust.diff <- 50.0
    #seed <- 1
    #tran_ids <- results$tran_id
    
    ######################################################################
    
    # Generate row names
    row.names(psi) <- psi$tran_id
    psi$tran_id <- NULL
    
    # Subset relevant cells
    psi.pheno <- psi.pheno[which(psi.pheno$sample.id %in% sample.ids), ]
    psi <- psi[, which(names(psi) %in% psi.pheno$sample.id)]

    # Subset events with sufficient cells
    . <- apply(psi, 1, function(x) {sum(!is.na(x))})
    index.keep <- . >= min.cells
    
    if(sum(index.keep) == 0) {
        
        message("No expressed events found")
        
        MarvelObject$Modality$Results <- NULL
        
        return(MarvelObject)
    
    }
    
    psi <- psi[index.keep, , drop=FALSE]
    psi.feature <- psi.feature[index.keep, , drop=FALSE]
    
    # Subset user-defined events
    if(!is.null(tran_ids[1])) {
        
        index <- which(psi.feature$tran_id %in% tran_ids)
        psi.feature <- psi.feature[index, , drop=FALSE]
        psi <- psi[psi.feature$tran_id, , drop=FALSE]
        
    }
    
    # Compute num. of cells analysed
    n.cells <- apply(psi, 1, function(x) {sum(!is.na(x))})
    psi.feature$n.cells <- n.cells
    
    # Retrieve parameters
        # Define function
        estbetaParams <- function(x) {

            # Convert to numeric
            values <- as.numeric(x)
            
            # Remove missing values
            values <- values[!is.na(values)]
            
            # Round off values
            #values <- round(values, digits=4)
            
            # Jitter exact values
            set.seed(seed)
            values[values==1 & !is.na(values)] <- runif(sum(values==1, na.rm=TRUE), min=0.98, max=0.9999)
            values[values==0 & !is.na(values)] <- runif(sum(values==0, na.rm=TRUE), min=0.0001, max=0.02)
                                        
            # Build model
            #model <- fitdist(data=values, distr="beta", method="mle")
            model <- tryCatch(fitdistrplus::fitdist(data=values, distr="beta", method="mle"), error=function(err) "Error")
            
            # Retrieve parameters and log-likelihood
            if(inherits(model, "fitdist", TRUE)==1) {
            
                alpha <- model$estimate[1]
                beta <- model$estimate[2]
                log.likelihood <- summary(model)$loglik
                variance <- var(values)
                params <- list(alpha=alpha, beta=beta, log.likelihood=log.likelihood, variance=variance)
                return(params)
                
            } else  {
            
                alpha <- NA
                beta <- NA
                log.likelihood <- NA
                variance <- NA
                params <- list(alpha=alpha, beta=beta, log.likelihood=log.likelihood, variance=variance)
                return(params)
            
            }
            
        }

        # Retrieve parameters
        param <- apply(psi, 1, estbetaParams)
            
        # For debugging
        #param <- NULL
        
        #for(i in 1:nrow(psi)) {
                   
           #param[[i]] <- estbetaParams(na.omit(as.numeric(psi[i,])))
           
           #message(i)
           
           
        #}
        
        #values <- as.numeric(psi[4397,])

        # Annotate parameters
        psi.feature$alpha <- as.numeric(sapply(param, function(x) {x[1]}))
        psi.feature$beta <- as.numeric(sapply(param, function(x) {x[2]}))
        psi.feature$log.likelihood <- as.numeric(sapply(param, function(x) {x[3]}))
        psi.feature$variance <- as.numeric(sapply(param, function(x) {x[4]}))
        
    # Re-assign modality for missing modalities
    if(sum(is.na(psi.feature$alpha)) >= 1) {
        
        # Annotate row no. for re-ordering later
        psi.feature$row.num. <- c(1:nrow(psi))
        
        # Subset events with missing modality
        tran_ids <- psi.feature[is.na(psi.feature$alpha), "tran_id"]
        psi.na <- psi[tran_ids, ]
        psi.feature.na <- psi.feature[which(psi.feature$tran_id %in% tran_ids), ]
        
        # Retrieve parameters
            # Define function
            estbetaParams <- function(x) {

                # Convert to numeric
                values <- as.numeric(x)
                
                # Remove missing values
                values <- values[!is.na(values)]
                
                # Round off values
                #values <- round(values, digits=4)
                
                # Jitter values
                set.seed(seed)
                values <- values + runif(n=length(values), min=0.0001, max=0.01)
                values[which(values >=1)] <- values[which(values >=1)] - runif(n=values[which(values >=1)], min=0.0001, max=0.01)
                                            
                # Build model
                #model <- fitdist(data=values, distr="beta", method="mle")
                model <- tryCatch(fitdistrplus::fitdist(data=values, distr="beta", method="mle"), error=function(err) "Error")
                
                # Retrieve parameters and log-likelihood
                if(inherits(model, "fitdist", TRUE)==1) {
                
                    alpha <- model$estimate[1]
                    beta <- model$estimate[2]
                    log.likelihood <- summary(model)$loglik
                    variance <- var(values)
                    params <- list(alpha=alpha, beta=beta, log.likelihood=log.likelihood, variance=variance)
                    return(params)
                    
                } else  {
                
                    alpha <- NA
                    beta <- NA
                    log.likelihood <- NA
                    variance <- NA
                    params <- list(alpha=alpha, beta=beta, log.likelihood=log.likelihood, variance=variance)
                    return(params)
                
                }
                
            }
            
            # Retrieve parameters
            param.na <- apply(psi.na, 1, estbetaParams)
            
            # Annotate parameters
            psi.feature.na$alpha <- as.numeric(sapply(param.na, function(x) {x[1]}))
            psi.feature.na$beta <- as.numeric(sapply(param.na, function(x) {x[2]}))
            psi.feature.na$log.likelihood <- as.numeric(sapply(param.na, function(x) {x[3]}))
            psi.feature.na$variance <- as.numeric(sapply(param.na, function(x) {x[4]}))
            
            # Reorder results table
            psi.feature <- psi.feature[!is.na(psi.feature$alpha), ]
            psi.feature <- rbind.data.frame(psi.feature, psi.feature.na)
            psi.feature <- psi.feature[order(psi.feature$row.num.), ]
            psi.feature$row.num. <- NULL
            
        }
    
    # Report persistent missing modalities
    modality.na <- sum(is.na(psi.feature$alpha))
    
    if(modality.na >= 1) {
        
        message(paste("Modality couldn't be computer for ", modality.na, " events. Please re-run using a different seed number.", sep=""))
        
    }

    # Assign modalities
        # Create new column
        psi.feature$modality <- NA
        
        # Indicate missing values
        psi.feature$modality[which(is.na(psi.feature$alpha))] <- "Missing"
        
        # Bimodal
        psi.feature$modality[which(is.na(psi.feature$modality) &
                                   (psi.feature$alpha <= 0.4 | psi.feature$beta <= 0.4)
                                   )] <- "Bimodal"
        
        # Included (alpha > 2, beta < 1)
        psi.feature$modality[which(is.na(psi.feature$modality) &
                                   psi.feature$alpha >= 2.0 &
                                   psi.feature$beta <= 1
                                   )] <- "Included"
        
        # Included (by FC)
        psi.feature$modality[which(is.na(psi.feature$modality) &
                                   (psi.feature$alpha/psi.feature$beta) > 2.0
                                   )] <- "Included"
        
        # Included (beta > 2, alpha < 1)
        psi.feature$modality[which(is.na(psi.feature$modality) &
                                   psi.feature$beta >= 2.0 &
                                   psi.feature$alpha <= 1
                                   )] <- "Excluded"
                                   
        # Excluded (by FC)
        psi.feature$modality[which(is.na(psi.feature$modality) &
                                   (psi.feature$beta/psi.feature$alpha) > 2.0
                                   )] <- "Excluded"

            
        # Middle
        psi.feature$modality[which(is.na(psi.feature$modality) &
                                   psi.feature$alpha >= 1.6 &
                                   psi.feature$beta >= 1.6
                                   )] <- "Middle"
                                   
        # Multimodal
        psi.feature$modality[which(is.na(psi.feature$modality))] <- "Multimodal"
        

    # Further stratify included and excluded modalities
        # Create new column for new modality
        psi.feature$modality.var <- NA
        
        # Indicate missing values
        psi.feature$modality.var[which(is.na(psi.feature$alpha))] <- "Missing"
        
        # Included
        psi.feature$modality.var[which(is.na(psi.feature$modality.var) &
                                           psi.feature$variance <= sigma.sq &
                                           psi.feature$modality=="Included"
                                           )] <- "Included.Primary"
        
        
        psi.feature$modality.var[which(is.na(psi.feature$modality.var) &
                                           psi.feature$modality=="Included"
                                           )] <- "Included.Dispersed"
        
        # Excluded
        psi.feature$modality.var[which(is.na(psi.feature$modality.var) &
                                           psi.feature$variance <= sigma.sq &
                                           psi.feature$modality=="Excluded"
                                           )] <- "Excluded.Primary"
        
        
        psi.feature$modality.var[which(is.na(psi.feature$modality.var) &
                                           psi.feature$modality=="Excluded"
                                           )] <- "Excluded.Dispersed"

        # Non-included/excluded
        psi.feature$modality.var[which(is.na(psi.feature$modality.var))] <-
        psi.feature$modality[which(is.na(psi.feature$modality.var))]
        
        ########################################################################
        ############################ BIMODAL ADJUST ############################
        ########################################################################
        
        if(length(which(psi.feature$modality.var=="Bimodal")) != 0) {
        
            if(bimodal.adjust==TRUE) {
            
                # Compute feature
                pct.lower <- apply(psi, 1, function(x) {sum(x[which(!is.na(x))] < 0.25) / length(x[which(!is.na(x) & (x < 0.25 | x > 0.75))]) * 100})
                pct.higher <- apply(psi, 1, function(x) {sum(x[which(!is.na(x))] > 0.75) / length(x[which(!is.na(x) & (x < 0.25 | x > 0.75))]) * 100})
                psi.feature$pct.fc <- ifelse(pct.higher > pct.lower, pct.higher/pct.lower, pct.lower/pct.higher)
                psi.feature$pct.diff <- ifelse(pct.higher > pct.lower, pct.higher - pct.lower, pct.lower - pct.higher)
                psi.feature$psi.average <- apply(psi, 1, function(x) {mean(x[which(!is.na(x))])})
                
                # Split into bimodal/non-bimodal
                bi <- psi.feature[which(psi.feature$modality.var=="Bimodal"), ]
                non.bi <- psi.feature[which(psi.feature$modality.var!="Bimodal"), ]
                
                # Annotate true/false bimodal (MOST IMPORTANT STEP)
                bi$bimodal.class <- ifelse( bi$alpha <= 0.4 & bi$beta <= 0.4 &
                                            bi$pct.fc <= bimodal.adjust.fc &
                                            bi$pct.diff <= bimodal.adjust.diff
                                            , "pass", "fail"
                                            )
                
                # Reclassify false bimodals
                if(length(which(bi$bimodal.class=="fail")) != 0) {
                
                    # Subset false bimodals
                    bi.fail <- bi[which(bi$bimodal.class=="fail"), ]
                    
                    # Assign modalities
                    bi.fail$modality.bimodal.adj <- NA
                
                    # Included
                    bi.fail$modality.bimodal.adj[which(bi.fail$psi.average >= 0.5 &
                                                       bi.fail$variance <= sigma.sq
                                                       )] <- "Included.Primary"
                    
                    bi.fail$modality.bimodal.adj[which(is.na(bi.fail$modality.bimodal.adj) &
                                                     bi.fail$psi.average >= 0.5
                                                     )] <- "Included.Dispersed"
                                                                        
                    # Excluded
                    bi.fail$modality.bimodal.adj[which(is.na(bi.fail$modality.bimodal.adj) &
                                                       bi.fail$psi.average < 0.5 &
                                                       bi.fail$variance <= sigma.sq
                                                       )] <- "Excluded.Primary"
                    
                    bi.fail$modality.bimodal.adj[which(is.na(bi.fail$modality.bimodal.adj) &
                                                     bi.fail$psi.average < 0.5
                                                     )] <- "Excluded.Dispersed"
                                                                                                        
                # Merge bi-pass/fail
                    # Format bi-pass columns to match bi-fail
                    bi.pass <- bi[which(bi$bimodal.class=="pass"), ]
                    bi.pass$modality.bimodal.adj <- bi.pass$modality.var
                    
                    # Merge
                    bi <- rbind.data.frame(bi.pass, bi.fail)
                
            }
            
            # Merge bi/non-bi
                if(nrow(non.bi) != 0) {
                    
                    # Format non-bi columns to match bi
                    non.bi$bimodal.class <- NA
                    non.bi$modality.bimodal.adj <- non.bi$modality.var
                    
                    # Merge
                    psi.feature <- rbind.data.frame(bi, non.bi)
                    
                    # Reorder as per psi data frame
                    row.names(psi.feature) <- psi.feature$tran_id
                    psi.feature <- psi.feature[row.names(psi), ]
                    row.names(psi.feature) <- NULL
                    
                } else {
                    
                    # Merge
                    psi.feature <- bi
                    
                    # Reorder as per psi data frame
                    row.names(psi.feature) <- psi.feature$tran_id
                    psi.feature <- psi.feature[row.names(psi), ]
                    row.names(psi.feature) <- NULL
                    
                }
            
            }
                
        }
    
    # Recode missing values
    if(length(names(psi.feature)[which(names(psi.feature)=="modality.bimodal.adj")]) == 0) {
        
            # Special case: 1 row only, bimodal=pass
            
            psi.feature$modality.bimodal.adj <- psi.feature$modality.var
        
        } else {
            
            # Normal cases
            psi.feature$modality[which(psi.feature$modality=="Missing")] <- NA
            psi.feature$modality.var[which(psi.feature$modality.var=="Missing")] <- NA
            
            if(bimodal.adjust==TRUE & length(which(psi.feature$modality.var=="Bimodal")) == 0) {
            
                psi.feature$modality.bimodal.adj <- psi.feature$modality.var
            
            } else {
                
                psi.feature$modality.bimodal.adj[which(psi.feature$modality.bimodal.adj=="Missing")] <- NA
                
            }
        
    }
    
    ########################################################################
    
    # Remove intermediate columns
    #psi.feature$alpha <- NULL
    #psi.feature$beta <- NULL
    #psi.feature$log.likelihood <- NULL
    #psi.feature$variance <- NULL
    #psi.feature$pct.fc <- NULL
    #psi.feature$pct.diff <- NULL
    #psi.feature$psi.average <- NULL
    #psi.feature$bimodal.class <- NULL
    
    # Save to new slots
    MarvelObject$Modality$Results <- psi.feature
    
    return(MarvelObject)
            
}
