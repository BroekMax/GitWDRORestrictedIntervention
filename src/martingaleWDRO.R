####################################################################
# Code that has the functions for the Martingale estimator
####################################################################
pkgs=c("quadprog")
for(package in pkgs){
  if(!require(package,character.only = TRUE)) {install.packages(package)}
  library(package, character.only = TRUE)
}

####################################################################
# Function computing the Martingale estimator
# Inputs:
# X: Design matrix
# y: Response vector
# A: Matrix t(I-B)%*%(I-B)
# epsilon: Radius ambiguity set, functions as regularisation parameter
# Output:
# beta_estimated: The estimated coefficient 
####################################################################
#Martinale WDRO estimator
tikhonov_constrained=function(X,y,A,epsilon){
  n=nrow(X)
  #Concatenate covariates and response
  xy=cbind(X,y)
  k=ncol(xy)
  #Constraint matrix
  Amat=matrix(0,nrow=k,ncol=1)
  Amat[k,1]=1
  #Last coefficient corresponding to the response should be -1
  bvec=-1
  #Matrix for objective function
  Dmat=t(xy)%*%xy+epsilon*solve(A)
  dvec=rep(0,k)
  #Solve the QCQP, constraint is equality
  solution=solve.QP(Dmat=Dmat,dvec=dvec,Amat=Amat,bvec=bvec,meq=1)
  #Retrieve value
  beta_estimated=c(solution$solution)
  
  #Return beta without the constrained value
  return(beta_estimated[-k])
}