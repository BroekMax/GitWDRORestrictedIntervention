####################################################################
# Function to compute the Wasserstein causal regulariser
# Inputs: 
# X: Design matrix of observational data
# y: Response vector of observational data
# X_e: List of design matrices of environments
# y_e: List of response vectors of all environments
# beta_init: Initialisation of parameter
# gamma: Value of regularisation parameter gamma
# n_iter: Number of iterations used for updating
# Output:
# beta: Coefficient vector
####################################################################
wasserstein_causal_regulariser=function(X,y,X_e,y_e,beta_init=NULL,gamma=1,n_iter=50){
  #Number of dimensions d
  d=ncol(X)
  
  #Initialising list that stores the ordered design and response according to the size of the residual 
  resid_e=list()
  X_e_order=list()
  y_e_order=list()
  
  #Initialising beta
  if (is.null(beta_init)) {
    beta=numeric(d)
  } else {
    beta=beta_init
  }
  
  #Weight in each environment, here uniform
  weights=rep(1/(length(X_e)),length(X_e))
  #Update beta and residuals in turns for n_iter iterations
  for(j in 1:n_iter){
    residual=y-X%*%beta
    #Order residual
    resid=sort(residual)
    #Store design and response according to the order
    X_order=X[order(residual),]
    y_order=y[order(residual)]
    #Initialise helper matrices needed for matrix product computation
    X_helper=matrix(0,d,d)
    cross_helper=rep(0,d)
    Xe_helper=matrix(0,d,d)
    crosse_helper=rep(0,d)
    for(i in 1:length(X_e)){
      #Order residuals in each environment
      resid_env=y_e[[i]]-X_e[[i]]%*%beta
      resid_e[[i]]=sort(resid_env)
      #Store X,y according to residuals for each environment
      X_e_order[[i]]=X_e[[i]][order(resid_env),]
      y_e_order[[i]]=y_e[[i]][order(resid_env)]
      #Compute cross terms needed for the estimator
      X_helper=X_helper+weights[i]*t(X_e_order[[i]]-X_order)%*%(X_e_order[[i]]-X_order)
      cross_helper=cross_helper+weights[i]*t(X_e_order[[i]]-X_order)%*%(y_e_order[[i]]-y_order)
      Xe_helper=Xe_helper+t(X_e[[i]])%*%X_e[[i]]
      crosse_helper=crosse_helper+t(X_e[[i]])%*%y_e[[i]]
    }
    
    #Obtain solution for the given residuals
    beta=solve(t(X)%*%X+Xe_helper+gamma*X_helper)%*%(t(X)%*%y+crosse_helper+gamma*cross_helper)
  }
  return(beta)
}

####################################################################
# Function to compute benchmark causal regulariser
# Inputs: 
# X: Design matrix of observational data
# y: Response vector of observational data
# X_e: List of design matrices of environments
# y_e: List of response vectors of all environments
# beta_init: Initialisation of parameter
# gamma: Value of regularisation parameter gamma
# Output:
# beta: Coefficient vector
#################################################################### 
causal_regulariser=function(X,y,X_e,y_e,gamma=0.99){
  d=ncol(X)
  Gplus=matrix(0,d,d)
  GDelta=matrix(0,d,d)
  Zplus=rep(0,d)
  ZDelta=rep(0,d)
  
  weights=rep(1/(length(X_e)),length(X_e))
  #Compute the matrices by adding the products in each environment
  for(i in 1:length(X_e)){
    Gplus=Gplus+weights[i]^2*t(X_e[[i]])%*%X_e[[i]]
    GDelta=GDelta+weights[i]^2*t(X_e[[i]])%*%X_e[[i]]
    Zplus=Zplus+weights[i]^2*t(X_e[[i]])%*%y_e[[i]]
    ZDelta=ZDelta+weights[i]^2*t(X_e[[i]])%*%y_e[[i]]
  }
  #Add/subtract the observational terms
  Gplus=Gplus+t(X)%*%X
  GDelta=GDelta-t(X)%*%X
  Zplus=Zplus+t(X)%*%y
  ZDelta=ZDelta-t(X)%*%y
  #Compute closed form solution
  beta=solve(Gplus+gamma*GDelta)%*%(Zplus+gamma*ZDelta)
  return(beta)
}