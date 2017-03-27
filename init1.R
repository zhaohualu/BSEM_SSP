###################need careful revision when model changed

InitialValues <- list()

# Initial value for the loading matrix
InitialValues$Para$LoadingMatrix<-matrix(c(
1,0,0.1,
1,0,0,
1,0,0,
0.1,1,0,
0,1,0,
0,1,0,
0,0.1,1,
0,0,1,
0,0,1),
ncol=NK,nrow=NY,byr=T)

# Initial value for the intercepts
InitialValues$Para$Intercept<-rep(0,NY)

# Initial value for the unique variance
InitialValues$Para$UniqueVariance<-rep(0.1,NY)

# Initial value for the factor covariance matrix
PHI<-matrix(0.5,nrow=NZ,ncol=NZ)
diag(PHI)<-1.0
InitialValues$Para$FactorCovariance <- PHI

# Initial value for the latent variables
InitialValues$LatentVariable<-mvrnorm(N,rep(0,NZ),Sigma=PHI)


