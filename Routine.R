log.likelihood <- function(y, ly, omega, ispsx, ciphi) {
    log.like1 <- colSums((ispsx * (y - ly %*% omega))^2)
    
    log.like3 <- colSums((ciphi %*% omega)^2)
    
    return((log.like1 + log.like3))
}

MCMCSummary <- function(BFAResult) {
    
    MCMCEstimate <- list()
    
    MCMCEstimate$LoadingMatrix$PosteriorMean <- apply(BFAResult$LoadingMatrix, 
        FUN = mean, MAR = c(2, 3))
    NK <- dim(MCMCEstimate$LoadingMatrix$PosteriorMean)[2]
    MCMCEstimate$Intercept$PosteriorMean <- apply(BFAResult$Intercept, 
        FUN = mean, MAR = c(2))
    MCMCEstimate$UniqueVariance$PosteriorMean <- apply(BFAResult$UniqueVariance, 
        FUN = mean, MAR = c(2))
    MCMCEstimate$FactorCovariance$PosteriorMean <- matrix(apply(BFAResult$FactorCovariance, 
        FUN = mean, MAR = c(2)), NK, NK)
    
    
    MCMCEstimate$LoadingMatrix$SE <- apply(BFAResult$LoadingMatrix, 
        FUN = sd, MAR = c(2, 3))
    MCMCEstimate$Intercept$SE <- apply(BFAResult$Intercept, FUN = sd, 
        MAR = c(2))
    MCMCEstimate$UniqueVariance$SE <- apply(BFAResult$UniqueVariance, 
        FUN = sd, MAR = c(2))
    MCMCEstimate$FactorCovariance$SE <- matrix(apply(BFAResult$FactorCovariance, 
        FUN = sd, MAR = c(2)), NK, NK)
    
    MCMCEstimate$LoadingMatrix$PIP <- apply(BFAResult$PIP, FUN = mean, 
        MAR = c(2, 3))
    
    return(MCMCEstimate)
}

