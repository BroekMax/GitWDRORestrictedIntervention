####################################################################
# This code contains the functions related to anchor transformations
####################################################################


pkgs=c("MASS", "glmnet", "CVXR", "RPtests")
for(package in pkgs){
  if(!require(package,character.only = TRUE)) {install.packages(package)}
  library(package, character.only = TRUE)
}
suppressWarnings(library(CVXR, warn.conflicts=FALSE))

####################################################################
# Function computing the PH-Anchor-RWP estimator
# Inputs:
# yanch: Anchor transformed response data
# Xanch: Anchor transformed design matrix
# delta: PH parameter
# Outputs:
# coef: Coefficients
# a0: The intercept
####################################################################
#PH-Anchor-RWP
pseudo_huber_anchor_RWP=function(yanch, Xanch, delta){
  #Draw standard normals to compute the variance of Z
  e=rnorm(100000)
  #Intialise coefficients and intercept
  betaHat=Variable(ncol(Xanch))
  a0=Variable(1)
  #Approximation of covariance scaler
  scaler=mean(e^2/(e^2+delta^2))
  Sigma=as.matrix(Xanch)%*%t(as.matrix(Xanch)) 
  #Draw samples from Z
  Z=mvrnorm(n=100,mu=rep(0,dim(Sigma)[1]),delta^2*Sigma*scaler)
  #Dual norm of cost function is infinity-norm
  Z2=apply(abs(Z),1,max)
  #Ambiguity radius according to RWP
  epsilon=quantile(Z2,0.95)/sqrt(nrow(Xanch))
  #Objective function
  objective=Minimize(delta^2*(mean((cvxr_norm((cbind((yanch-Xanch%*%betaHat-a0)/delta,rep(1,nrow(Xanch)))),axis=1))-1))
                     +delta*epsilon*cvxr_norm(betaHat,1))
  problem=Problem(objective)
  #Default solver ECOS can fail to converge
  result=solve(problem,solver="SCS")
  #Retrieve coefficients
  coef=result$getValue(betaHat)
  a0=result$getValue(a0)
  return(list(coef=coef,a0=a0))
}

####################################################################
# Function computing the PH-Anchor estimator
# Inputs:
# yanch: Anchor transformed response data
# Xanch: Anchor transformed design matrix
# delta: PH parameter
# epsilon: size of the ambiguity set used in the regularisation parameter
# Outputs:
# coef: Coefficients
# a0: The intercept
####################################################################
#PH-Anchor, epsilon is assumed to be given
pseudo_huber_anchor=function(yanch, Xanch, delta,epsilon){
  #Intialise coefficients and intercept
  betaHat=Variable(ncol(Xanch))
  a0=Variable(1)
  #Objective function
  objective=Minimize(delta^2*(mean((cvxr_norm((cbind((yanch-Xanch%*%betaHat-a0)/delta,rep(1,nrow(Xanch)))),axis=1))-1))
                     +delta*epsilon*cvxr_norm(betaHat,1))
  problem=Problem(objective)
  #Default solver ECOS can fail to converge
  result=solve(problem,solver="SCS")
  #Retrieve coefficients
  coef=result$getValue(betaHat)
  a0=result$getValue(a0)
  return(list(coef=coef,a0=a0))
}

####################################################################
# Function computing the anchor regression estimator
# Inputs:
# yanch: Anchor transformed response data
# Xanch: Anchor transformed design matrix
# lambdas: Grid for CV
# Outputs:
# anchor_out: glmnet model
####################################################################
#Anchor regression
anchor_regression = function(yanch, Xanch,lambdas) {
  #CV to find lambda
  lambda=cv.glmnet(as.matrix(Xanch),as.matrix(yanch),lambda=lambdas)$lambda.min
  #Fit anchor; Lasso on transformed data
  anchor_out = glmnet(as.matrix(Xanch), as.matrix(yanch), lambda = lambda)
  return(anchor_out)
}

####################################################################
# Function computing Huber from a WDRO point of view
# Inputs:
# yanch: Anchor transformed response data
# Xanch: Anchor transformed design matrix
# delta_h: PH parameter
# Outputs:
# coef_h: Coefficients
# a0_h: The intercept
####################################################################
#Huber loss with Anchor transformation
huber_anchor=function(yanch, Xanch, delta_h){
  #Confidence
  alpha=0.05
  delta=delta_h
  #Intialising
  betaHat=Variable(ncol(Xanch))
  a0=Variable(1)
  zi=Variable(nrow(Xanch))
  #Regularisation parameter via RWP
  epsilon=pi/(pi-2)*qnorm(1-alpha/(2*ncol(Xanch)))/sqrt(nrow(Xanch))
  #Objective function
  objective=Minimize(mean(1/2*zi^2)+
                       delta*mean(abs(-yanch+Xanch%*%betaHat-zi+a0))+
                       delta*epsilon^2*cvxr_norm(betaHat,1))
  problem=Problem(objective)
  #Default solver ECOS can fail to converge
  result=solve(problem,solver="SCS")
  #Retrieve coefficients
  coef_h=result$getValue(betaHat)
  a0_h=result$getValue(a0)
  return(list(coef_h=coef_h,a0_h=a0_h))
}

####################################################################
# Function computing the Sqrt-Lasso-RWP estimator
# Inputs:
# yanch: Anchor transformed response data
# Xanch: Anchor transformed design matrix
# Outputs:
# lasso: Sqrt-lasso model
####################################################################
#Sqrt-lasso with RWP regularisation parameter
sqrt_lasso_RWP=function(yanch, Xanch){
  #Confidence level
  alpha=0.05
  #Regularisation parameter according to RWP, here named epsilon
  epsilon=pi/(pi-2)*qnorm(1-alpha/(2*ncol(Xanch)))/sqrt(nrow(Xanch))
  suppressWarnings(lasso <- sqrt_lasso(as.matrix(Xanch), yanch, lambda=epsilon, output_all = TRUE))
  return(lasso)
}