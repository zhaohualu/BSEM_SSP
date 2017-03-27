# Put the identification of intercepts and loading matrix together in a big matrix
# Put the matrices of all candidate models in a list
CandidateModels <- list()

# Identification for the elements in the loading matrix
# 0 - fixed and known common in all candidate models
# 1 - free parameter common in all candidate models
# 2 - parameters fixed (usually 0) in a specific candidate models
# 3 - free parameters in a specific candidate models
# SSP were assigned to all elements with 2 or 3 


CandidateModels$M0<-matrix(c(  
  1,1,0,2,
  1,0,0,0,
  1,1,0,0,  
  1,2,1,0,
  1,0,0,0,
  1,0,1,0,  
  1,0,2,1,
  1,0,0,0,
  1,0,0,1
),nrow=NY,ncol=NK+1,byr=T)

CandidateModels$M1<-matrix(c(
  1,1,0,3,
  1,0,0,0,
  1,1,0,0,  
  1,3,1,0,
  1,0,0,0,
  1,0,1,0,  
  1,0,3,1,
  1,0,0,0,
  1,0,0,1  
),nrow=NY,ncol=NK+1,byr=T)

