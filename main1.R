# Load some libraries and R functions
install.packages("mvtnorm")
install.packages("loo")


library(stats)
library(MASS)
library(mvtnorm)
library(loo)


source("Routine.R",local=T)
source("BFA.R")
source("MCCRoutine.R")

# Dimension Setting
NMCMC <- 8000  # Number of MCMC samples generated
NBurn <- 4000  # Number of MCMC samples discarded as burn-in samples
Nthin <- 1     # thinning of MCMC sample, the number of MCMC iterations is NMCMC * Nthin 
N <- 100 # Number of subjects
NY <- 9  # Number of items
NK <- 3  # Number of factors
NZ <- NK

# Simulate data

set.seed(1000)

source("SimulateData.R")


#set.seed(1000)

# Input data
DataMat<-read.table("SimulatedData.txt",header=F)
# Define the hyperparameters of  prior distributions
source("Prior.R",local=T)
# Define the identification conditions
source("ind.R",local=T)
# Define the initial values of parameters and latent factors
source("init1.R",local=T)

#set.seed(10)
# Generate MCMC samples for the factor analysis model
# The results is a list containing MCMC samples of parameters and latent factors
BFAResult <- BFA(DataMat,InitialValues,Identification,Hyperpara,NMCMC=NMCMC,NBurn=NBurn,Nthin=Nthin,StoreLatentVariables = T)

#dump(c("BFAResult"),file="MCMCSamples.R")
save("BFAResult",file="MCMCSamples.RData")

# Calculate posterior mean, SE and posterior inclusion probabilities of the parameters stored in a list
MCMCSummary(BFAResult)

#set.seed(10)
# Calculate the model comparison statistics
MCCResults <- BayesianMCC(BFAResult,DataMat,Hyperpara,Identification,BridgeSamplerMargLik=T,BIC=T,DIC=T,LOO=T)


# Only use if SSP is used
# Store the structure of the candidate models in a list CandidateModels
source("CandidateModelsIdentification.R")
# Calculate the estimated marginal model probability by SSP
MMP <- MMP_SSP(BFAResult$PIP,CandidateModels)







