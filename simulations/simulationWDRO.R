####################################################################
# Working directory should be the root of the code folder
####################################################################
source("./src/martingaleWDRO.R") #Martingale Estimator
source("./src/iterative_scheme.R") #Iterative scheme WDRORI Estimator

####################################################################
# This is a simulation from Section 6.2 with the Martingale and Iterative scheme
####################################################################

#Sample Size
n=100
#Test size
ntest=1000
#10 covariates
d=10
k=d+1
B = matrix(0, k, k)
#DGP initialisation
beta_true=rep(0,d)
beta_true[9]=1
beta_true[10]=-2
pa=which(beta_true!=0)
B[k,-k]=beta_true
#Identity matrix
I=diag(k)

#Matrix for estimators, essentially the cost function squared 
A=t(I-B)%*%(I-B)

#Number of simulations
nsim=100
#Vectors to store the MSE
MSE_OLS=rep(0,nsim)
MSE_causal=rep(0,nsim)
MSE_WDRO=rep(0,nsim)
MSE_M1=rep(0,nsim)
MSE_M2=rep(0,nsim)
MSE_W1=rep(0,nsim)
MSE_W2=rep(0,nsim)
MSE_matrix=matrix(0,4,7)
#Half of the sample is the interventional distribution
frac=0.5

set.seed(200008)
#For each setting do 100 repetitions
for(setting in 1:4){
  #Choice of epsilon discussed in the paper
  epsilon=2000*(setting==1)+50*(setting>=2)
  epsilon_2=2000*(setting==1)+50*(setting>=2)
  epsilon_3=2
  epsilon_sgd=0.5*(setting==1)+0.25*(setting>1)
  for(i in 1:nsim){
    #Standard normal covariates
    X=matrix(rnorm(n * d,0,1), n, d)

    #Shifts only in setting 4
    Xtest=matrix((rnorm(ntest * d,0,1)+(setting>=4)*rnorm(ntest*d,sqrt(0.5),1)), ntest, d) 
    X[ceiling(frac*n):n,9]=X[ceiling(frac*n):n,9]+rnorm(n-ceiling(frac*n)+1,0,5)*(setting>=3)
    
    #Latent confounder from setting 2 onwards
    H=rnorm(n,0,1)
    Htest=rnorm(ntest,0,1)
    X=X+H*(setting>=2)
    Xtest=Xtest+Htest*(setting>=2)
    
    #In case propagation is inserted
    for(j in 1:d){
      X[,j]=B[j,1:d]%*%t(X)+X[,j]     
    }
    for(j in 1:d){
      Xtest[,j]=B[j,1:d]%*%t(Xtest)+Xtest[,j]    
    }
    
    #Response variable
    y=as.numeric(B[k,-k]%*%t(X))+rnorm(n,0,1)+H*(setting>=2)
    ytest=as.numeric(B[k,-k]%*%t(Xtest))+rnorm(ntest,0,1)+Htest*(setting>=2)
    Xtest=sweep(Xtest,2,colMeans(X),FUN="-")
    
    #Standardising
    X=sweep(X,2,colMeans(X),FUN="-")
    ytest=ytest-mean(y)
    y=y-mean(y)
    
    #Parent oracle matrix 
    B_2= matrix(0, k, k)
    causal=lm(y~X[,pa])
    B_2[k,pa]=causal$coefficients[-1]
    ols=lm(y~X)
    #OLS
    beta_ols=ols$coefficients[-1]
    A_2=t(I-B_2)%*%(I-B_2)

    #Causal estimator
    causal=lm(y~X[,pa])
    beta_estimated_OLS=ols$coefficients[-1]
    ypred_OLS=Xtest%*%beta_estimated_OLS
    MSE_OLS[i]=(mean((ypred_OLS-ytest)^2))
    beta_estimated_causal=causal$coefficients[-1]
    ypred_causal=Xtest[,pa]%*%beta_estimated_causal
    MSE_causal[i]=(mean((ypred_causal-ytest)^2))
    
    #Martingale estimator with graph oracle
    beta_estimated_M1=tikhonov_constrained(X,y,A,epsilon)
    #Martingale estimator with parent oracle
    beta_estimated_M2=tikhonov_constrained(X,y,A_2,epsilon_2)

    #Classic WDRO estimator
    beta_estimated_WDRO=tikhonov_constrained(X,y,I,epsilon_3)
    #Predictions
    ypred_M1=Xtest%*%beta_estimated_M1
    ypred_M2=Xtest%*%beta_estimated_M2
    ypred_WDRO=Xtest%*%beta_estimated_WDRO
    
    #Store prediction errors
    MSE_M1[i]=(mean((ypred_M1-ytest)^2))
    MSE_M2[i]=(mean((ypred_M2-ytest)^2))
    MSE_WDRO[i]=(mean((ypred_WDRO-ytest)^2))
    
    #Iterative scheme solver for graph- and parent oracle
    beta_estimated_W1=iterative_scheme(X,y,A,150000,burn_in=0,beta_init=beta_ols,epsilon=epsilon_sgd)$beta[-k] 
    beta_estimated_W2=iterative_scheme(X,y,A_2,150000,burn_in=0,beta_init=beta_ols,epsilon=epsilon_sgd)$beta[-k]    
    
    ypred_W1=Xtest%*%beta_estimated_W1
    MSE_W1[i]=(mean((ypred_W1-ytest)^2))
    ypred_W2=Xtest%*%beta_estimated_W2
    MSE_W2[i]=(mean((ypred_W2-ytest)^2))

  }
  #Record MSE
  MSE_matrix[setting,1]=mean(MSE_OLS)
  MSE_matrix[setting,2]=mean(MSE_causal)
  MSE_matrix[setting,3]=mean(MSE_WDRO)
  MSE_matrix[setting,4]=mean(MSE_M1)
  MSE_matrix[setting,5]=mean(MSE_M2)
  MSE_matrix[setting,6]=mean(MSE_W1)
  MSE_matrix[setting,7]=mean(MSE_W2)
}

#Print results
print(MSE_matrix)
