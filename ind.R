###################need careful revised when changing model

Identification <- list()

# Identification for the elements in the loading matrix
# 0 - the element is fixed and known
# 1 - the element is a free parameter to be estimated with usually Bayesian CFA approach
# 2 - the element is a free parameter and a SSP is assigned to determine the inclusion/exclusion of the parameter and estimate the model probabilities of the related models
Identification$LoadingMatrix<-matrix(c(
0,0,1,
1,0,0,
1,0,0,
1,0,0,
0,1,0,
0,1,0,
0,1,0,
0,0,1,
0,0,1
),ncol=NK,nrow=NY,byr=T)

# Identification of the intercepts
# Set to 1 (free) or 0 (no intercept)
Identification$Intercept<-rep(1,NY)


