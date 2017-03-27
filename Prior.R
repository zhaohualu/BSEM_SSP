###############attention must be given to structure.

Hyperpara <- list()

# Set the mean of the normal distribution for the elements in the loading matrix
Hyperpara$LoadingMatrix$Mean <-matrix(c(
1,0,0.1,
1,0,0,
1,0,0,
0.1,1,0,
0,1,0,
0,1,0,
0,0.1,1,
0,0,1,
0,0,1
),nrow=NY,ncol=NK,byr=T)

# Set the 1/variance of the normal distribution for the elements in the loading matrix
Hyperpara$LoadingMatrix$precision <- matrix(c(
  1,0,1,
  1,0,0,
  1,0,0,
  1,1,0,
  0,1,0,
  0,1,0,
  0,1,1,
  0,0,1,
  0,0,1
),nrow=NY,ncol=NK,byr=T)

# Prior probability of the nonzero mixture component of the loading matrix, 0 if ssp is not used
pr = array(0,dim=c(NY,NK))
pr[1,3]=0.5
pr[4,1]=0.5
pr[7,2]=0.5
Hyperpara$LoadingMatrix$PIP <- pr

# Set the mean and 1/variance of the normal distribution for the intercepts
Hyperpara$Intercept$Mean<-rep(0,NY)
Hyperpara$Intercept$precision<-rep(1,NY)
# Prior probability of the nonzero mixture component of the intercept, 0 if ssp is not used
Hyperpara$Intercept$PIP <- rep(0,NY)

# Hyperparameters of Unique Variance
# rate and shape parameters of the inverse gamma distribution
Hyperpara$UniqueVariance$alpha <- rep(11,NY)
Hyperpara$UniqueVariance$beta <- rep(3,NY)


# Hyperparameters of factor covaraince
# Set the scale and matrix of the inverse Wishart distribution
Hyperpara$FactorCovariance$Scalar <- 6
TmpMat <- matrix(0.5,NZ,NZ)
diag(TmpMat) <- 1
Hyperpara$FactorCovariance$Matrix <- TmpMat*(Hyperpara$FactorCovariance$Scalar-NK-1)


