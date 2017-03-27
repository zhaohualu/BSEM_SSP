###########pay attention to structure



LY<-matrix(c(
  1,0,0.0,
  1,0,0,
  1,0,0,
  0.0,1,0,
  0,1,0,
  0,1,0,
  0,0.0,1,
  0,0,1,
  0,0,1
),ncol=NK,nrow=NY,byr=T)

MU<-rep(0,NY)

PSX<-rep(0.3,NY)

PHI<-matrix(0.5,nrow=NZ,ncol=NZ)
diag(PHI)<-1

#gen True Latent Variable
Omega<-t(mvrnorm(N,rep(0,NZ),Sigma=PHI))

Y <- array(0,dim=c(NY,N))

theta<-LY%*%Omega+MU
for(j in 1:NY){
    Y[j,]<-rnorm(N,theta[j,],sqrt(PSX[j]))
}

write(Y,"SimulatedData.txt",ncol=NY,append=F)

