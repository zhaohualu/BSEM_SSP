BFA <- function(DataMat, InitialValues, Identification, Hyperpara, 
    NMCMC = 10000, NBurn = 3000, Nthin = 1, StoreLatentVariables = F) {
    
    pb <- txtProgressBar(style = 3)
    
    # Some prelimary data management for MCMC sampling
    
    # Data
    Y <- t(DataMat)
    
    # Parameters
    LY <- InitialValues$Para$LoadingMatrix
    PSX <- InitialValues$Para$UniqueVariance
    PHI <- InitialValues$Para$FactorCovariance
    MU <- InitialValues$Para$Intercept
    
    # Latent variables
    Omega <- t(InitialValues$LatentVariable)
    
    # Dimensions
    N <- dim(DataMat)[1]
    NY <- dim(DataMat)[2]
    NK <- dim(LY)[2]
    NZ <- NK
    
    # Hyperparameters
    PLY <- Hyperpara$LoadingMatrix$Mean
    PPLY <- Hyperpara$LoadingMatrix$precision
    
    PMU <- Hyperpara$Intercept$Mean
    PPMU <- Hyperpara$Intercept$precision
    
    alphax <- Hyperpara$UniqueVariance$alpha
    betax <- Hyperpara$UniqueVariance$beta
    
    RouZero <- Hyperpara$FactorCovariance$Scalar
    IWMat <- Hyperpara$FactorCovariance$Matrix
    
    prLY <- Hyperpara$LoadingMatrix$PIP
    prMU <- Hyperpara$Intercept$PIP
    
    # Identification
    IDY <- Identification$LoadingMatrix
    IDMU <- Identification$Intercept
    
    # Convenient variables
    iPSX <- 1/PSX
    inv.sqrt.PSX <- 1/sqrt(PSX)
    
    inv.PHI <- chol2inv(chol(PHI))
    c.inv.PHI <- chol(inv.PHI)
    
    AC.XI <- 0
    
    # Null interface
    NM <- 0
    
    NANA <- 0
    AD <- array(0, dim = c(NY, NANA))
    AZ <- array(0, dim = c(NANA, N))
    
    PAD <- array(0, dim = c(NY, NANA))
    PPAD <- array(0, dim = c(NY, NANA))
    
    IDA <- array(0, dim = c(NY, NANA))
    
    # Joint matrix
    COE <- cbind(MU, AD, LY)
    PCOE <- cbind(matrix(PMU, ncol = 1), PAD, PLY)
    PPCOE <- cbind(matrix(PPMU, ncol = 1), PPAD, PPLY)
    ICE <- cbind(IDMU, IDA, IDY)
    IDE <- array(0, dim = dim(ICE))
    IDE[ICE == 1] <- 1
    COV <- rbind(rep(1, N), AZ, Omega)
    pr <- cbind(prMU, prLY)
    
    # Array storing MCMC samples
    Nrec <- NMCMC - NBurn
    EILY <- array(0, dim = c(Nrec, NY, NK + NANA + 1))
    ELY2 <- array(0, dim = c(Nrec, NY, NK + NANA + 1))
    EPSX <- array(0, dim = c(Nrec, NY))
    EPHI <- array(0, dim = c(Nrec, (NZ * NZ)))
    if (StoreLatentVariables) {
        EOmega <- array(0, dim = c(Nrec, NK, N))
    }
    
    # Start MCMC sampling
    
    for (g in 1:NMCMC) {
        for (gthin in 1:Nthin) {
            gm <- g - NBurn
            
            # Samleing the latent factor
            if (NK > 0) {
                ISG <- crossprod(inv.sqrt.PSX * LY)
                if (NZ > 0) 
                  ISG[(NM + 1):NK, (NM + 1):NK] <- ISG[(NM + 1):NK, 
                    (NM + 1):NK] + inv.PHI
                ISG <- ISG  #*VAR1
                
                SIG <- chol2inv(chol(ISG))
                cSIG <- chol(SIG)
                OMEN <- Omega + crossprod(cSIG, matrix(rnorm(N * 
                  NK, 0, 1), nrow = NK, ncol = N))
                
                Ycen <- MU + AD %*% AZ
                
                r1 <- log.likelihood(Y - Ycen, LY, Omega, inv.sqrt.PSX, 
                  c.inv.PHI)
                r2 <- log.likelihood(Y - Ycen, LY, OMEN, inv.sqrt.PSX, 
                  c.inv.PHI)
                r <- exp(0.5 * (r1 - r2))
                comr <- runif(N)
                crit <- (comr < r)
                COV[(2 + NANA):(1 + NANA + NK), crit] <- Omega[, 
                  crit] <- OMEN[, crit]
                AC.XI <- AC.XI + crit
            }
            
            
            # Sampleing the indicator in SSP, intercepts, loading matrix and
            # unique variance for every item
            
            XX <- tcrossprod(COV)
            
            for (j in 1:NY) {
                
                # Sampling the R for SSP in case some elements assigned SSP
                
                subs <- ICE[j, ] == 0  # fixed item
                Ycen <- as.vector(Y[j, , drop = FALSE] - COE[j, 
                  (subs), drop = F] %*% COV[(subs), , drop = F])
                XY <- COV %*% Ycen
                Ycens2 <- sum(Ycen^2)
                
                iPSigb0 <- diag(PPCOE[j, ])
                Pmean <- PCOE[j, ]
                
                idx0 <- IDE[j, ] > 0  # current !=0 components
                
                idxtmp <- which(ICE[j, ] == 2)
                for (k in idxtmp) {
                  if (IDE[j, k] > 0) {
                    idx1 <- idx0
                    idx2 <- idx0
                    idx2[k] <- F
                  } else {
                    idx1 <- idx0
                    idx2 <- idx0
                    idx1[k] <- T
                  }
                  
                  iPSigb1 <- iPSigb0[idx1, idx1, drop = F]
                  Pb1 <- Pmean[idx1]
                  iSigb1 <- XX[idx1, idx1] + iPSigb1
                  Sigb1 <- solve(iSigb1)
                  bets1 <- Sigb1 %*% (XY[idx1] + iPSigb1 %*% Pb1)
                  
                  # Calcualte the posterior proability of the indicators
                  
                  logR <- log(pr[j, k]/(1 - pr[j, k]))
                  
                  if (sum(idx2) > 0) {
                    iPSigb2 <- iPSigb0[idx2, idx2, drop = F]
                    Pb2 <- Pmean[idx2]
                    iSigb2 <- XX[idx2, idx2] + iPSigb2
                    Sigb2 <- solve(iSigb2)
                    bets2 <- Sigb2 %*% (XY[idx2] + iPSigb2 %*% Pb2)
                    
                    logR <- logR - log(sqrt(det(Sigb2)))  #*det(iPSigb2)
                    logR <- logR + (N/2 + alphax[j]) * log(Ycens2 + 
                      sum(Pb2 * (iPSigb2 %*% Pb2)) - sum(bets2 * 
                      (iSigb2 %*% bets2)))
                  }
                  
                  
                  logR <- logR - (N/2 + alphax[j]) * log(Ycens2 + 
                    sum(Pb1 * (iPSigb1 %*% Pb1)) - sum(bets1 * (iSigb1 %*% 
                    bets1)))
                  logR <- logR + log(sqrt(det(Sigb1)))  #det(iPSigb1)
                  if (PPCOE[j, k] > 0) 
                    logR <- logR + log(sqrt(PPCOE[j, k]))
                  R <- exp(logR)
                  p <- 1 - 1/(1 + R)
                  itmp <- rbinom(1, 1, p)
                  IDE[j, k] <- itmp
                  idx0 <- IDE[j, ] > 0  # current !=0 components    
                }
                idx1 <- IDE[j, ] == 0 & ICE[j, ] == 2
                COE[j, idx1] <- 0
                
                # Sampling the intercepts, loadings and unique variances
                
                subs <- IDE[j, ] > 0  # current active components to be optimized.    
                Ycen <- COE[j, (!subs), drop = F] %*% COV[(!subs), 
                  , drop = F]
                
                omesub <- COV[subs, , drop = FALSE]
                ssubs <- sum(subs)
                iPSiginv <- diag(PPCOE[j, subs], ssubs)
                Pmean1 <- PCOE[j, subs]
                
                # Sampling unique variances
                Ycen <- as.vector(Y[j, , drop = FALSE] - Ycen)
                alphastar <- alphax[j] + N/2
                betastar <- betax[j] + 1/2 * sum(Ycen^2)
                if (ssubs > 0) {
                  calsmnpsx <- chol2inv(chol((XX[subs, subs, drop = F] + 
                    iPSiginv)))
                  temp <- (omesub %*% Ycen + iPSiginv %*% Pmean1)
                  LYnpsx <- calsmnpsx %*% temp
                  betastar <- betastar + 1/2 * (sum(Pmean1 * iPSiginv * 
                    Pmean1) - sum(temp * LYnpsx))
                }
                iPSX[j] <- rgamma(1, shape = alphastar, rate = betastar)
                PSX[j] <- 1/iPSX[j]
                inv.sqrt.PSX[j] <- sqrt(iPSX[j])
                
                # Sampling the intercepts and loadings
                if (ssubs > 0) {
                  COE[j, idx0] <- mvrnorm(1, LYnpsx, PSX[j] * calsmnpsx)
                }
                
            }
            
            MU[] <- COE[, 1, drop = F]
            if (NANA > 0) 
                AD[, ] <- COE[, 2:(1 + NANA), drop = F]
            if (NK > 0) 
                LY[, ] <- COE[, (2 + NANA):(1 + NANA + NK), drop = F]
            
            
            # Sampling the factor covariance matrix
            if (NZ > 0) {
                
                inv.PHI[, ] <- rWishart(1, RouZero + N, solve(tcrossprod(Omega[(NM + 
                  1):NK, , drop = F]) + IWMat))
                c.inv.PHI <- chol(inv.PHI)
                PHI <- chol2inv(c.inv.PHI)
                
                # if(FLAG){ for(j in 1:NZ){
                # inv.PHI[j,j]<-rgamma(1,N/2+5,4+0.5*sum(Omega[j+NM,]^2))
                # PHI[j,j] = 1/inv.PHI[j,j] } c.inv.PHI<-chol(inv.PHI) }
                
                
            }
            
            # After burn-in, store the MCMC samples
            if (gm > 0) {
                EILY[gm, , ] <- IDE
                ELY2[gm, , ] <- COE
                
                if (NZ > 0) 
                  EPHI[gm, ] <- as.vector(PHI)
                EPSX[gm, ] <- PSX
                if (StoreLatentVariables) {
                  EOmega[gm, , ] <- Omega
                }
            }
            
            
            setTxtProgressBar(pb, g/NMCMC)
        }
    }  #end of MCMC sampling
    
    # Put the MCMC samples in a list and return
    
    MCMCResult <- list(LoadingMatrix = ELY2[, , (2 + NANA):(1 + 
        NANA + NK), drop = F], PIP = EILY, Intercept = ELY2[, , 
        1], UniqueVariance = EPSX, FactorCovariance = EPHI)
    if (StoreLatentVariables) {
        MCMCResult$LatentVariable <- EOmega
    }
    return(MCMCResult)
    
}

