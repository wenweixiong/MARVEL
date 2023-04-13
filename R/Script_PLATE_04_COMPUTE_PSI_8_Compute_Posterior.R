#' @title Compute PSI with Bayesian approach
#'
#' @description Combines a prior with likelihood (observed SJ counts) to obtain a posterior value of PSI.
#'
#' @details This function computes posterior value of PSI based on a prior and likelihood (observed SJ counts). The resulting PSI values may be subsequently used for downstream principal component analysis (PCA). Because a prior is used here, there are no arguments to set minimum coverage threshold for computing the PSI values.
#'
#' @param MarvelObject S3 object generated from \code{ComputePSI} function.
#'
#' @return An object of class S3 with new slots \code{$Counts}.
#'
#' @import methods
#'
#' @export

ComputePSI.Posterior <- function(MarvelObject) {

    # Define arguments
    # None
    
    # Example arguments
    #MarvelObject <- marvel
    
    ######################################################################
    ######################### RETRIEVE COUNTS ############################
    ######################################################################
    
    # SE
    if("SE" %in% names(MarvelObject$Counts)) {
        
        # Track progress
        print("Computing posterior PSI for SE...")
        
        # Subset SJ counts
        df.included.1 <- MarvelObject$Counts$SE$sj.included.1
        df.included.2 <- MarvelObject$Counts$SE$sj.included.2
        df.excluded <- MarvelObject$Counts$SE$sj.excluded
        
        # Add row names
        row.names(df.included.1) <- df.included.1$tran_id
        df.included.1$tran_id <- NULL
        
        row.names(df.included.2) <- df.included.2$tran_id
        df.included.2$tran_id <- NULL
        
        row.names(df.excluded) <- df.excluded$tran_id
        df.excluded$tran_id <- NULL
        
        # Compute observed counts
            # Numerator
            df.num.obs <- df.included.1 + df.included.2

            # Denominator
            df.den.obs <- (df.included.1 + df.included.2) + (2 * df.excluded)
            
        # Compute total counts
            # Numerator
            num.total <- apply(df.num.obs, 1, function(x) {sum(x, na.rm=TRUE)})

            # Denominator
            den.total <- apply(df.den.obs, 1, function(x) {sum(x, na.rm=TRUE)})

        # Compute pseudo counts
            # Numerator
            num.prior <- num.total/den.total
            
            # Denominator
            den.prior <- den.total/num.total

        # Compute PSI
            # Add pseudo counts to observed counts: Numerator
            df.num.final <- sweep(df.num.obs, 1, num.prior, "+")
         
            # Add pseudo counts to observed counts: Denominator
            df.den.final <- sweep(df.den.obs, 1, den.prior, "+")

            # Compute
            df.psi <- df.num.final/df.den.final
            
        # Add tran_id column
        . <- data.frame("tran_id"=row.names(df.psi))
        df.psi <- cbind.data.frame(., df.psi)
        row.names(df.psi) <- NULL
            
        # Save into list
        MarvelObject$PSI.Posterior$SE <- df.psi
                        
    }
    
    # MXE
    if("MXE" %in% names(MarvelObject$Counts)) {
        
        # Track progress
        print("Computing posterior PSI for MXE...")
        
        # Subset SJ counts
        df.included.1 <- MarvelObject$Counts$MXE$sj.included.1
        df.included.2 <- MarvelObject$Counts$MXE$sj.included.2
        df.excluded.1 <- MarvelObject$Counts$MXE$sj.excluded.1
        df.excluded.2 <- MarvelObject$Counts$MXE$sj.excluded.2
        
        # Add row names
        row.names(df.included.1) <- df.included.1$tran_id
        df.included.1$tran_id <- NULL
        
        row.names(df.included.2) <- df.included.2$tran_id
        df.included.2$tran_id <- NULL
        
        row.names(df.excluded.1) <- df.excluded.1$tran_id
        df.excluded.1$tran_id <- NULL
        
        row.names(df.excluded.2) <- df.excluded.2$tran_id
        df.excluded.2$tran_id <- NULL
        
        # Compute observed counts
            # Numerator
            df.num.obs <- df.included.1 + df.included.2

            # Denominator
            df.den.obs <- (df.included.1 + df.included.2) + (df.excluded.1 + df.excluded.2)
 
     # Compute total counts
         # Numerator
         num.total <- apply(df.num.obs, 1, function(x) {sum(x, na.rm=TRUE)})

         # Denominator
         den.total <- apply(df.den.obs, 1, function(x) {sum(x, na.rm=TRUE)})

     # Compute pseudo counts
         # Numerator
         num.prior <- num.total/den.total
         
         # Denominator
         den.prior <- den.total/num.total

     # Compute PSI
         # Add pseudo counts to observed counts: Numerator
         df.num.final <- sweep(df.num.obs, 1, num.prior, "+")
      
         # Add pseudo counts to observed counts: Denominator
         df.den.final <- sweep(df.den.obs, 1, den.prior, "+")

         # Compute
         df.psi <- df.num.final/df.den.final

     # Add tran_id column
     . <- data.frame("tran_id"=row.names(df.psi))
     df.psi <- cbind.data.frame(., df.psi)
     row.names(df.psi) <- NULL
         
     # Save into list
     MarvelObject$PSI.Posterior$MXE <- df.psi
 
    }
    
    # RI
    if("RI" %in% names(MarvelObject$Counts)) {
        
        # Track progress
        print("Computing posterior PSI for RI...")
        
        # Subset SJ counts
        df.included <- MarvelObject$Counts$RI$counts.included
        df.excluded <- MarvelObject$Counts$RI$counts.excluded
        
        # Add row names
        row.names(df.included) <- df.included$tran_id
        df.included$tran_id <- NULL

        row.names(df.excluded) <- df.excluded$tran_id
        df.excluded$tran_id <- NULL
        
        # Compute observed counts
            # Numerator
            df.num.obs <- df.included

            # Denominator
            df.den.obs <- df.included + df.excluded
            
        # Compute total counts
            # Numerator
            num.total <- apply(df.num.obs, 1, function(x) {sum(x, na.rm=TRUE)})

            # Denominator
            den.total <- apply(df.den.obs, 1, function(x) {sum(x, na.rm=TRUE)})

        # Compute pseudo counts
            # Numerator
            num.prior <- num.total/den.total
            
            # Denominator
            den.prior <- den.total/num.total

        # Compute PSI
            # Add pseudo counts to observed counts: Numerator
            df.num.final <- sweep(df.num.obs, 1, num.prior, "+")
         
            # Add pseudo counts to observed counts: Denominator
            df.den.final <- sweep(df.den.obs, 1, den.prior, "+")

            # Compute
            df.psi <- df.num.final/df.den.final
  
        # Add tran_id column
        . <- data.frame("tran_id"=row.names(df.psi))
        df.psi <- cbind.data.frame(., df.psi)
        row.names(df.psi) <- NULL
  
        # Save into list
        MarvelObject$PSI.Posterior$RI <- df.psi
        
    }
    
    # A5SS
    if("A5SS" %in% names(MarvelObject$Counts)) {
        
        # Track progress
        print("Computing posterior PSI for A5SS...")
        
        # Subset SJ counts
        df.included <- MarvelObject$Counts$A5SS$sj.included
        df.excluded <- MarvelObject$Counts$A5SS$sj.excluded
        
        # Add row names
        row.names(df.included) <- df.included$tran_id
        df.included$tran_id <- NULL

        row.names(df.excluded) <- df.excluded$tran_id
        df.excluded$tran_id <- NULL
        
        # Compute observed counts
            # Numerator
            df.num.obs <- df.included

            # Denominator
            df.den.obs <- df.included + df.excluded

        # Compute total counts
            # Numerator
            num.total <- apply(df.num.obs, 1, function(x) {sum(x, na.rm=TRUE)})

            # Denominator
            den.total <- apply(df.den.obs, 1, function(x) {sum(x, na.rm=TRUE)})

        # Compute pseudo counts
            # Numerator
            num.prior <- num.total/den.total
            
            # Denominator
            den.prior <- den.total/num.total

        # Compute PSI
            # Add pseudo counts to observed counts: Numerator
            df.num.final <- sweep(df.num.obs, 1, num.prior, "+")
         
            # Add pseudo counts to observed counts: Denominator
            df.den.final <- sweep(df.den.obs, 1, den.prior, "+")

            # Compute
            df.psi <- df.num.final/df.den.final

        # Add tran_id column
        . <- data.frame("tran_id"=row.names(df.psi))
        df.psi <- cbind.data.frame(., df.psi)
        row.names(df.psi) <- NULL
            
        # Save into list
        MarvelObject$PSI.Posterior$A5SS <- df.psi
        
    }
    
    # A3SS
    if("A3SS" %in% names(MarvelObject$Counts)) {
        
        # Track progress
        print("Computing posterior PSI for A3SS...")
        
        # Subset SJ counts
        df.included <- MarvelObject$Counts$A3SS$sj.included
        df.excluded <- MarvelObject$Counts$A3SS$sj.excluded
        
        # Add row names
        row.names(df.included) <- df.included$tran_id
        df.included$tran_id <- NULL

        row.names(df.excluded) <- df.excluded$tran_id
        df.excluded$tran_id <- NULL
        
        # Compute observed counts
            # Numerator
            df.num.obs <- df.included

            # Denominator
            df.den.obs <- df.included + df.excluded

        # Compute total counts
            # Numerator
            num.total <- apply(df.num.obs, 1, function(x) {sum(x, na.rm=TRUE)})

            # Denominator
            den.total <- apply(df.den.obs, 1, function(x) {sum(x, na.rm=TRUE)})

        # Compute pseudo counts
            # Numerator
            num.prior <- num.total/den.total
            
            # Denominator
            den.prior <- den.total/num.total

        # Compute PSI
            # Add pseudo counts to observed counts: Numerator
            df.num.final <- sweep(df.num.obs, 1, num.prior, "+")
         
            # Add pseudo counts to observed counts: Denominator
            df.den.final <- sweep(df.den.obs, 1, den.prior, "+")

            # Compute
            df.psi <- df.num.final/df.den.final

        # Add tran_id column
        . <- data.frame("tran_id"=row.names(df.psi))
        df.psi <- cbind.data.frame(., df.psi)
        row.names(df.psi) <- NULL
            
        # Save into list
        MarvelObject$PSI.Posterior$A3SS <- df.psi
    
    }

    # AFE
    if("AFE" %in% names(MarvelObject$Counts)) {
        
        # Track progress
        print("Computing posterior PSI for AFE...")
        
        # Subset SJ counts
        df.included <- MarvelObject$Counts$AFE$sj.included
        df.excluded <- MarvelObject$Counts$AFE$sj.excluded
        
        # Add row names
        row.names(df.included) <- df.included$tran_id
        df.included$tran_id <- NULL

        row.names(df.excluded) <- df.excluded$tran_id
        df.excluded$tran_id <- NULL
        
        # Compute observed counts
            # Numerator
            df.num.obs <- df.included

            # Denominator
            df.den.obs <- df.included + df.excluded
 
         # Compute total counts
             # Numerator
             num.total <- apply(df.num.obs, 1, function(x) {sum(x, na.rm=TRUE)})

             # Denominator
             den.total <- apply(df.den.obs, 1, function(x) {sum(x, na.rm=TRUE)})

         # Compute pseudo counts
             # Numerator
             num.prior <- num.total/den.total
             
             # Denominator
             den.prior <- den.total/num.total

         # Compute PSI
             # Add pseudo counts to observed counts: Numerator
             df.num.final <- sweep(df.num.obs, 1, num.prior, "+")
          
             # Add pseudo counts to observed counts: Denominator
             df.den.final <- sweep(df.den.obs, 1, den.prior, "+")

             # Compute
             df.psi <- df.num.final/df.den.final
 
         # Add tran_id column
         . <- data.frame("tran_id"=row.names(df.psi))
         df.psi <- cbind.data.frame(., df.psi)
         row.names(df.psi) <- NULL
 
         # Save into list
         MarvelObject$PSI.Posterior$AFE <- df.psi
                
    }
    
    # ALE
    if("ALE" %in% names(MarvelObject$Counts)) {
        
        # Track progress
        print("Computing posterior PSI for ALE...")
        
        # Subset SJ counts
        df.included <- MarvelObject$Counts$ALE$sj.included
        df.excluded <- MarvelObject$Counts$ALE$sj.excluded
        
        # Add row names
        row.names(df.included) <- df.included$tran_id
        df.included$tran_id <- NULL

        row.names(df.excluded) <- df.excluded$tran_id
        df.excluded$tran_id <- NULL
        
        # Compute observed counts
            # Numerator
            df.num.obs <- df.included

            # Denominator
            df.den.obs <- df.included + df.excluded
        
        # Compute total counts
            # Numerator
            num.total <- apply(df.num.obs, 1, function(x) {sum(x, na.rm=TRUE)})

            # Denominator
            den.total <- apply(df.den.obs, 1, function(x) {sum(x, na.rm=TRUE)})

        # Compute pseudo counts
            # Numerator
            num.prior <- num.total/den.total
            
            # Denominator
            den.prior <- den.total/num.total

        # Compute PSI
            # Add pseudo counts to observed counts: Numerator
            df.num.final <- sweep(df.num.obs, 1, num.prior, "+")
         
            # Add pseudo counts to observed counts: Denominator
            df.den.final <- sweep(df.den.obs, 1, den.prior, "+")

            # Compute
            df.psi <- df.num.final/df.den.final
 
        # Add tran_id column
        . <- data.frame("tran_id"=row.names(df.psi))
        df.psi <- cbind.data.frame(., df.psi)
        row.names(df.psi) <- NULL
         
        # Save into list
        MarvelObject$PSI.Posterior$ALE <- df.psi
        
    }
    
    
    ######################################################################
    ###################### RETURN FINAL OBJECTS ##########################
    ######################################################################
    
    # New slots already created above
    return(MarvelObject)
        
}
