#' Log Likelihood function of factor model
LFYTK <- function(Y, mu, Lambda, PSX, PHI) {
    DY <- dim(Y)
    N <- DY[1]
    P <- DY[2]
    
    SIG <- tcrossprod(Lambda %*% PHI, Lambda) + PSX
    
    return(sum(dmvnorm(Y, mu, SIG, log = T)))
}

#' Log density function of inverse gamma distribution
digamma <- function(x, shape, rate, log = T) {
    return(shape * log(rate) - lgamma(shape) - (shape + 1) * log(x) - 
        rate/x)
}

# ` Log density of multivariate gamma function
mgf <- function(p, a, logflag = T) {
    lp <- 0.25 * p * (p - 1) * log(pi) + sum(lgamma(a + (1 - (1:p))/2))
    if (logflag) {
        return(lp)
    } else {
        return(exp(lp))
    }
}

#' Log density function of inverse Wishart distribution without the constants
diwish0 <- function(PHI, v, PHI0, logflag = T) {
    p <- dim(PHI0)[1]
    
    dtmp <- -0.5 * (v + p + 1) * log(det(PHI)) - 0.5 * (sum(PHI0 * 
        solve(PHI)))
    if (logflag) {
        return(dtmp)
    } else {
        return(exp(dtmp))
    }
}

#' Log density function of inverse Wishart distribution
diwish <- function(PHI, v, PHI0, FUN = diwish0, MAR = c(1, 2), logflag = T) {
    p <- dim(PHI0)[1]
    const <- 0.5 * v * log(det(PHI0)) - (0.5 * v * p * log(2) + 
        mgf(p, v/2, TRUE))
    # lp0 = apply(PHI,FUN=FUN,MAR=MAR,v=v,PHI0=PHI0,logflag=logflag)
    lp0 <- diwish0(PHI, v, PHI0, logflag = logflag)
    lp <- const + lp0
    if (logflag) {
        return(lp)
    } else {
        return(exp(lp))
    }
}

#' Log density function of loading matrix (multivariate normal distribution) 
gf1 <- function(Lambda, MLam, VLam) {
    return((dmvnorm(Lambda, MLam, VLam, log = T)))
}

#' Log density function of unique variance (inverse gamma distribution) 
gf2 <- function(PSX, shape, rate) {
    return(sum(digamma(PSX, shape, rate, log = T)))
}

# ` log density function of loading matrix (proposal
# distribution)
LPT1 <- function(Lambda, MLam, VLam) {
    return(sum(dnorm(as.vector(Lambda), as.vector(MLam), sqrt(as.vector(VLam)), 
        log = T)))
}

# ` log density function of unique variance (inverse gamma
# distribution)
LPT2 <- function(PSX, a, b) {
    return(sum(digamma(as.vector(PSX), as.vector(a), rate = (as.vector(b)), 
        log = T)))
}


BayesianMCC <- function(BFAResult, DataMat, Hyperpara, Identification, 
    BridgeSamplerMargLik = T, BIC = T, DIC = T, LOO = T) {
    
    MCCResults <- c()
    MCCNames <- c()
    
    # Identification
    
    lamfreeidx <- which(t(Identification$LoadingMatrix == 1))
    
    # Hyperparameters
    PLY <- Hyperpara$LoadingMatrix$Mean
    PPLY <- Hyperpara$LoadingMatrix$precision
    
    PMU <- Hyperpara$Intercept$Mean
    PPMU <- Hyperpara$Intercept$precision
    
    alphax <- Hyperpara$UniqueVariance$alpha
    betax <- Hyperpara$UniqueVariance$beta
    
    RouZero <- Hyperpara$FactorCovariance$Scalar
    IWMat <- Hyperpara$FactorCovariance$Matrix
    
    ## CALCULATE MEAN AND COVARIANCE
    
    Nrec <- dim(BFAResult$LoadingMatrix)[1]
    NY <- dim(BFAResult$LoadingMatrix)[2]
    NK <- dim(BFAResult$LoadingMatrix)[3]
    NZ <- NK
    NM <- 0
    NANA <- 0
    
    mcmc.omega <- array(0, dim = c(Nrec, NK * N))
    mcmc.mu <- array(0, dim = c(Nrec, NY))
    mcmc.lam <- array(0, dim = c(Nrec, NY * NK))
    mcmc.sigma <- array(0, dim = c(Nrec, NY))
    mcmc.phx <- array(0, dim = c(Nrec, NZ * NZ))
    
    # ASSIGN
    mcmc.omega[, ] <- BFAResult$LatentVariable[, (NM + 1):(NM + 
        NZ), ]
    
    for (j in 1:Nrec) {
        mcmc.lam[j, ] <- t(BFAResult$LoadingMatrix[j, , ])
    }
    mcmc.mu[, ] <- BFAResult$Intercept
    mcmc.sigma[, ] <- BFAResult$UniqueVariance
    mcmc.phx[, ] <- BFAResult$FactorCovariance
    
    
    # Calcualte Posterior mean of parameters
    Emu <- colMeans(BFAResult$Intercept)
    ELam <- apply(BFAResult$LoadingMatrix, MAR = c(2, 3), FUN = mean)
    Esigma <- colMeans(BFAResult$UniqueVariance)
    Ephx <- matrix(colMeans(BFAResult$FactorCovariance), nrow = NZ, 
        ncol = NZ, byrow = T)
    
    phxidx <- matrix(1:(NZ * NZ), nrow = NZ)
    phxidx1 <- phxidx[lower.tri(phxidx, diag = TRUE)]
    mcmc.theta <- cbind(mcmc.mu, mcmc.lam[, lamfreeidx], mcmc.sigma, 
        mcmc.phx[, phxidx1])
    Etheta <- colMeans(mcmc.theta)
    Vtheta <- cov(mcmc.theta)
    
    
    if (BridgeSamplerMargLik) {
        
        #### LOG LIKELIHOOD BASED ON MCMC SAMPLES (mlfytk)
        Y1 <- DataMat
        mlfytk <- numeric(nrow(mcmc.lam))
        for (j in 1:nrow(mcmc.lam)) {
            Lambj <- matrix(mcmc.lam[j, ], nrow = NY, ncol = NK, 
                byrow = T)
            PHIj <- matrix(mcmc.phx[j, ], nrow = NZ, ncol = NZ, 
                byrow = T)
            Sigmaj <- diag(mcmc.sigma[j, ])
            
            mlfytk[j] <- LFYTK(Y1, mcmc.mu[j, ], Lambj, Sigmaj, 
                PHIj)
        }
        
        #### LOG PRIOR DENSITY FOR MCMC SAMPLES (SLPT)
        tPLY <- t(PLY)
        tPPLY <- t(PPLY)
        # Prior density of Intercept
        SLPT0 <- apply(mcmc.mu, FUN = LPT1, MAR = 1, MLam = PMU, 
            VLam = 1/(PPMU))
        # Prior density of free elements in the loading matrix
        SLPT1 <- apply(mcmc.lam[, lamfreeidx], FUN = LPT1, MAR = 1, 
            MLam = tPLY[lamfreeidx], VLam = 1/(tPPLY[lamfreeidx]))
        
        SLPT2 <- numeric(Nrec)
        SLPT3 <- numeric(Nrec)
        for (j in 1:Nrec) {
            # Prior density of unique variance
            SLPT2[j] <- LPT2(mcmc.sigma[j, ], a = alphax, b = betax)
            # Prior density of factor covariance matrix
            PHItmp <- matrix(mcmc.phx[j, ], NZ, NZ, byrow = F)
            SLPT3[j] <- diwish(PHItmp, RouZero, IWMat, logflag = T)
        }
        SLPT <- SLPT0 + SLPT1 + SLPT2 + SLPT3
        
        
        
        ### Calculate the parameters of the g distibutions
        
        # Mean and Covariance of multivariate normal distribution for
        # the Intercepts
        mmu <- colMeans(mcmc.mu)
        vmu <- cov(mcmc.mu)
        
        # Mean and Covariance of multivariate normal distribution for
        # the free elements in the loading matrix
        mLam <- colMeans(mcmc.lam[, lamfreeidx])
        vLam <- cov(mcmc.lam[, lamfreeidx])
        
        # Shape and rate parameters of inverse gamma distribution for
        # the unique variances
        pm <- colMeans(mcmc.sigma)
        pv <- apply(mcmc.sigma, FUN = var, MAR = 2)
        pa <- as.numeric(pm^2/pv + 2)
        pb <- as.numeric((pa - 1) * pm)
        
        # Parameters of the inverse Wishart distribution for the factor
        # covariance matrix
        Exiv <- colMeans(mcmc.omega)
        Exi <- matrix(Exiv, nrow = NK, ncol = N)
        Exi2 <- tcrossprod(Exi) + IWMat
        Erho <- mean((Exi2)/Ephx) + NZ + 1
        if (Erho <= 0) {
            Erho <- mean(diag(Exi2)/diag(Ephx)) + NZ + 1
        }
        
        ## Denominator
        
        ### LOG DENSITY FUNCTION of proposal g distributions OF MCMC
        ### SAMPLES (Sgf)
        
        Sgf0 <- gf1(mcmc.mu, mmu, vmu)
        Sgf1 <- gf1(mcmc.lam[, lamfreeidx], mLam, vLam)
        Sgf2 <- apply(mcmc.sigma, FUN = gf2, MAR = 1, pa, pb)
        Sgf3 <- numeric(Nrec)
        for (j in 1:Nrec) {
            PHItmp <- matrix(mcmc.phx[j, ], NZ, NZ, byrow = F)
            Sgf3[j] <- diwish(PHItmp, Erho, Exi2, logflag = T)
        }
        Sgf <- Sgf0 + Sgf1 + Sgf2 + Sgf3
        
        
        
        
        #### GENERATE ADDITIONAL SAMPLES from the g density FOR BRIDGE
        #### SAMPLER
        
        gpsx <- matrix(1/rgamma(Nrec * length(pa), rep(pa, each = Nrec), 
            rate = rep(pb, each = Nrec)), nrow = Nrec, ncol = length(pm))
        gmu <- mvrnorm(Nrec, mmu, vmu)
        gLam <- mvrnorm(Nrec, mLam, vLam)
        giPHI <- rWishart(Nrec, Erho, solve(Exi2))
        gPHI <- array(0, dim = dim(giPHI))
        for (j in 1:Nrec) {
            gPHI[, , j] <- solve(giPHI[, , j])
        }
        
        #### Log LIKELIHOOD BASED ON ADDITIONAL (g) SAMPLES (gmlfytk)
        gmlfytk <- numeric(Nrec)
        for (j in 1:Nrec) {
            tLambj <- matrix(mcmc.lam[1, ], nrow = NK, ncol = NY)
            # tLambj[lamfixidx] = lamfixval
            tLambj[lamfreeidx] <- gLam[j, ]
            Lambj <- t(tLambj)
            PHIj <- gPHI[, , j]
            Sigmaj <- diag(gpsx[j, ])
            
            gmlfytk[j] <- LFYTK(Y1, gmu[j, ], Lambj, Sigmaj, PHIj)
        }
        
        #### Log PRIOR FOR ADDITIONAL (g) SAMPLES (SLPT)
        tPLY <- t(PLY)
        tPPLY <- t(PPLY)
        # Prior density of Intercept
        gSLPT0 <- apply(gmu, FUN = LPT1, MAR = 1, MLam = PMU, VLam = 1/(PPMU))
        # Prior density of free elements in the loading matrix
        gSLPT1 <- apply(gLam, FUN = LPT1, MAR = 1, MLam = tPLY[lamfreeidx], 
            VLam = 1/(tPPLY[lamfreeidx]))
        gSLPT2 <- numeric(Nrec)
        gSLPT3 <- numeric(Nrec)
        
        for (j in 1:Nrec) {
            # Prior density of unique variance
            gSLPT2[j] <- LPT2(gpsx[j, ], a = alphax, b = betax)
            # Prior density of factor covariance matrix
            PHItmp <- gPHI[, , j]
            gSLPT3[j] <- diwish(PHItmp, RouZero, IWMat, logflag = T)
        }
        gSLPT <- gSLPT0 + gSLPT1 + gSLPT2 + gSLPT3
        
        
        # g LOG DENSITY OF g SAMPLES (Sgf) Intercept
        gSgf0 <- gf1(gmu, mmu, vmu)
        # Free loading elements
        gSgf1 <- gf1(gLam, mLam, vLam)
        # Unique Variance
        gSgf2 <- apply(gpsx, FUN = gf2, MAR = 1, pa, pb)
        # Factor Covariance
        gSgf3 <- numeric(Nrec)
        for (j in 1:Nrec) {
            PHItmp <- gPHI[, , j]
            gSgf3[j] <- diwish(PHItmp, Erho, Exi2, logflag = T)
        }
        # 
        gSgf <- gSgf0 + gSgf1 + gSgf2 + gSgf3
        
        # Bridge Sampler
        lN <- 0.5 * (gmlfytk + gSLPT - gSgf)
        lD <- -0.5 * (mlfytk + SLPT - Sgf)
        mlN <- mean(lN)
        mlD <- mean(lD)
        Marglik_BS <- log(mean(exp(lN - mlN))/mean(exp(lD - mlD))) + 
            (mlN - mlD)
        
        MCCResults <- c(MCCResults, Marglik_BS)
        MCCNames <- c(MCCNames, "BayesFactor_BridgeSampler")
    }
    
    if (BIC) {
        Emu <- colMeans(mcmc.mu)
        ELam <- matrix(colMeans(mcmc.lam), nrow = NY, ncol = NK, 
            byrow = T)
        Esigma <- colMeans(mcmc.sigma)
        Ephx <- matrix(colMeans(mcmc.phx), nrow = NZ, ncol = NZ, 
            byrow = T)
        
        library(mvtnorm)
        # BIC, log likelihood and number of parameters
        BIC <- -2 * sum(dmvnorm(DataMat, Emu, ELam %*% Ephx %*% 
            t(ELam) + diag(Esigma), log = T)) + (length(lamfreeidx) + 
            NY + NY + NZ * (NZ + 1)/2) * log(N)
        MCCResults <- c(MCCResults, BIC)
        MCCNames <- c(MCCNames, "BIC")
    }
    
    if (DIC) {
        # Log likelihood at posterior mean
        dbt <- -2 * sum(dmvnorm(DataMat, Emu, ELam %*% Ephx %*% 
            t(ELam) + diag(Esigma), log = T))
        dtb <- 0
        # Posterior mean of log likelihood
        for (j in 1:Nrec) {
            Lambj <- matrix(mcmc.lam[j, ], nrow = NY, ncol = NK, 
                byrow = T)
            PHIj <- matrix(mcmc.phx[j, ], nrow = NZ, ncol = NZ, 
                byrow = T)
            Sigmaj <- diag(mcmc.sigma[j, ])
            
            dtb <- dtb + (-2 * sum(dmvnorm(DataMat, mcmc.mu[j, ], 
                Lambj %*% PHIj %*% t(Lambj) + (Sigmaj), log = T)))
        }
        dtb <- dtb/Nrec
        pd <- dtb - dbt
        DIC <- dtb + 2 * pd
        
        MCCResults <- c(MCCResults, DIC)
        MCCNames <- c(MCCNames, "DIC")
    }
    
    
    
    if (LOO) {
        ## Calculate LOO with the the loo function in the loo package
        ## Create a matrix recording the log likelihood of each subject
        ## givenparameters in each MCMC iteration , # of row is number of
        ## MCMC sample , # of column if # of subject
        LogLikMcMat <- matrix(nrow = Nrec, ncol = nrow(DataMat))
        for (j in 1:Nrec) {
            Lambj <- matrix(mcmc.lam[j, ], nrow = NY, ncol = NK, 
                byrow = T)
            PHIj <- matrix(mcmc.phx[j, ], nrow = NZ, ncol = NZ, 
                byrow = T)
            Sigmaj <- diag(mcmc.sigma[j, ])
            
            LogLikMcMat[j, ] <- dmvnorm(DataMat, mcmc.mu[j, ], Lambj %*% 
                PHIj %*% t(Lambj) + (Sigmaj), log = T)
        }
        lootmp <- loo(LogLikMcMat)
        LooEst <- lootmp$looic
        
        
        MCCResults <- c(MCCResults, LooEst)
        MCCNames <- c(MCCNames, "LOO-PSIS")
        
    }
    names(MCCResults) <- MCCNames
    print(MCCResults)
    return(MCCResults)
}

# Calculate the estimated marginal model probability by SSP
MMP_SSP <- function(PIP, CandidateModels) {
    
    # A vector record the estimated marginal model probability by
    # SSP
    TrueModPostProb <- numeric(length(CandidateModels))
    names(TrueModPostProb) <- names(CandidateModels)
    for (mcrep in 1:dim(PIP)[1]) {
        EILYtmp <- PIP[mcrep, , ]
        
        # Count the number of iteration where the candidate model is the
        # selected model
        for (j in 1:length(CandidateModels)) {
            # Index of free parameters to be compared
            Tidx <- CandidateModels[[j]] == 3
            # Index of fixed parameters to be compared
            Fidx <- CandidateModels[[j]] == 2
            # Count
            TrueModPostProb[j] <- TrueModPostProb[j] + (all(EILYtmp[Tidx] == 
                1) & all(EILYtmp[Fidx] == 0))
        }
    }
    print("Marginal model probability calculated by SSP")
    print(TrueModPostProb)
    return(TrueModPostProb)
}



