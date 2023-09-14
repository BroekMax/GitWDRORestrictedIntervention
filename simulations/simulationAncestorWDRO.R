####################################################################
# Working directory should be the root of the code folder
####################################################################
source("./src/martingaleWDRO.R") #Martingale Estimator
path="./simulations/results/"

pkgs=c("MASS","ggplot2","glmnet")
for(package in pkgs){
  if(!require(package,character.only = TRUE)) {install.packages(package)}
  library(package, character.only = TRUE)
}

####################################################################
# This is a simulation from Section 6.2 with Ancestor Regression
####################################################################

#vector of shifts
shifts=seq(0,10,by=0.5)
nshifts=length(shifts)
#Sample size
n=250
#Test size
ntest=1000
#Dimension
d=10
k=d+1
#DGP matrix with true coefficients
B = matrix(0, k, k)
B[k,1:10]=0
B[k,9]=-2
B[k,10]=1
#Causal parent interaction, X_9 is a causal parent of X_10
#Set to 0 to get the other simulation result
B[10,9]=1.5
I=diag(k)
A=t(I-B)%*%(I-B)
pa=which(B[k,-k]!=0)
nsim=100
#Vectors to store the MSE of each run
MSE_causal=rep(0,nsim)
MSE_OLS=rep(0,nsim)
MSE_M2=rep(0,nsim)
MSE_M3=rep(0,nsim)
MSE_WDRO=rep(0,nsim)
MSE_Lasso=rep(0,nsim)
MSE_matrix=matrix(0,nshifts,6)

set.seed(66+B[10,9]*4)
for(r in 1:nshifts){
  shift=shifts[r]
  for(i in 1:nsim){
    X=matrix(runif(n * d,-1,1), n, d)
    #Variance shift at test time
    Xtest=matrix((runif(ntest * d,-0.5,0.5)+runif(ntest*d,-shift,shift)), ntest, d) 
    #4 environments with variance interventions
    for(hh in 1:4){
      X[(hh*50+1):(hh*50+50),9]=X[(hh*50+1):(hh*50+50),9]+runif(50,-hh,hh)
    }
    
    #Propagation of the network
    for(j in 1:d){
      X[,j]=B[j,1:d]%*%t(X)+X[,j]     
    }
    for(j in 1:d){
      Xtest[,j]=B[j,1:d]%*%t(Xtest)+Xtest[,j]    
    }
    
    #Response variable
    y=as.numeric(B[k,-k]%*%t(X))+rgamma(n,1,1)
    ytest=as.numeric(B[k,-k]%*%t(Xtest))+rgamma(ntest,1,1)
    Xtest=sweep(Xtest,2,colMeans(X),FUN="-")
    
    #Standardise
    X=sweep(X,2,colMeans(X),FUN="-")
    ytest=ytest-mean(y)
    y=y-mean(y)
    #Leave two environments as validation sets
    samp=151:250
    Xval=X[samp,]
    Xtrain=X[-samp,]
    yval=y[samp]
    ytrain=y[-samp]
    
    #Ancestor regression
    ols=lm(y^3~X+y)
    #Bonferroni correction, find significant variables, these are the anctors
    parents=which(summary(ols)$coefficients[c(-1,-(k+1)),4]*10<0.05)
    B_2= matrix(0, k, k)
    #Matrix for Ancestor
    B_3= matrix(0, k, k)
    #If Ancestor Regression finds ancestors then fit on these ancestors to obtain the matrix B (here B_3)
    if(length(parents)>=1){
      epsilon_3=20
      causal_ancestor=lm(y~X[,parents])
      B_3[k,parents]=causal_ancestor$coefficients[-1]
      A_3=t(I-B_3)%*%(I-B_3)
      lowest=Inf
      #Cross validation for Martingale estimator if there are parents
      for(bb in seq(1,2000,by=10)){
        be=tikhonov_constrained(Xtrain,ytrain,A_3,bb)
        MSE_be=mean((yval-Xval%*%be)^2)
        if(MSE_be<lowest){
          lowest=MSE_be
          epsilon_3=bb
        }
      }
    }
    #Use OLS if no parents obtained
    else{
      B_3[k,-k]=rep(0,d)
      A_3=t(I-B_3)%*%(I-B_3)
      epsilon_3=0
    }
    #Causal estimator
    causal=lm(y~X[,pa])
    beta_estimated_causal=causal$coefficients[-1]
    #Set last row of B_2 to be the estimated causal coefficients
    B_2[k,pa]=causal$coefficients[-1]
    A_2=t(I-B_2)%*%(I-B_2)
    #CV for lambda used in Lasso
    lambda=cv.glmnet(as.matrix(X),as.matrix(y),intercept=FALSE)$lambda.min
    #Fit Lasso
    lasso_fit = glmnet(as.matrix(X), as.matrix(y),intercept=FALSE, lambda = lambda)
    #Lasso prediction
    ypred_Lasso=predict(lasso_fit,Xtest)
    
    
    #Pooled OLS
    ols=lm(y~X)
    beta_estimated_ols=ols$coefficients[-1]
    ypred_ols=Xtest%*%beta_estimated_ols
    MSE_OLS[i]=(mean((ypred_ols-ytest)^2))
    #Causal estimator uses the obtained parents
    if(length(parents>=1)){
      causal=lm(y~X[,parents])
      ypred_causal=as.matrix(Xtest[,parents])%*%causal$coefficients[-1]+causal$coefficients[1]
      MSE_causal[i]=(mean((ypred_causal-ytest)^2))
    }
    #If no parents then do standard OLS for causal estimator
    else{ypred_causal=Xtest%*%ols$coefficients[-1]+ols$coefficients[1]
    MSE_causal[i]=(mean((ypred_causal-ytest)^2))
    }
    #Cross-validation for Martingale parent oracle
    lowest=Inf
    for(bb in seq(1,2000,by=10)){
      be=tikhonov_constrained(Xtrain,ytrain,A_2,bb)
      MSE_be=mean((yval-Xval%*%be)^2)
      if(MSE_be<lowest){
        lowest=MSE_be
        epsilon_2=bb
      }
    }
    #Estimate and predict
    beta_estimated_M2=tikhonov_constrained(X,y,A_2,epsilon_2)
    ypred_M2=Xtest%*%beta_estimated_M2
    beta_estimated_M3=tikhonov_constrained(X,y,A_3,epsilon_3)
    ypred_M3=Xtest%*%beta_estimated_M3
    beta_estimated_WDRO=tikhonov_constrained(X,y,I,1)
    ypred_WDRO=Xtest%*%beta_estimated_WDRO
    #Compute MSE
    MSE_M2[i]=(mean((ypred_M2-ytest)^2))
    MSE_M3[i]=(mean((ypred_M3-ytest)^2))
    MSE_WDRO[i]=(mean((ypred_WDRO-ytest)^2))
    MSE_Lasso[i]=(mean((ypred_Lasso-ytest)^2))
  }
  #Average MSE for each shift
  MSE_matrix[r,1]=mean(MSE_OLS)
  MSE_matrix[r,2]=mean(MSE_causal)
  MSE_matrix[r,3]=mean(MSE_M2)
  MSE_matrix[r,4]=mean(MSE_M3)
  MSE_matrix[r,5]=mean(MSE_WDRO)
  MSE_matrix[r,6]=mean(MSE_Lasso)
  
}

MSE_matrix
#Smooth to make plot more visually informative
smoothedmatrix=MSE_matrix
shifts=seq(0,10,by=0.5)
for(b in 1:6){
  loess_fit=loess(smoothedmatrix[,b] ~ shifts,span=1)
  smoothedmatrix[,b]=predict(loess_fit)
}
result=as.data.frame(smoothedmatrix)
result$shifter=shifts

#Use ggplot to make the plot
p = ggplot(result, aes(x = shifter)) +
  geom_line(aes(y = V1, color = "Pooled OLS")) +
  geom_line(aes(y = V2, color = "Pooled Causal")) +
  geom_line(aes(y = V3, color = "Martingale")) +
  geom_line(aes(y = V5, color = "WDRO")) +
  geom_line(aes(y = V6, color = "Lasso")) +
  geom_line(aes(y = V4, color = "Ancestor Martingale")) +
  theme_minimal() +
  xlab("Perturbation strength") +
  ylab("Test MSE") +
  theme(axis.text.x = element_text(angle = 0, hjust = 1)) +
  scale_x_continuous(breaks = c(0,2.5,5,7.5,10))+
  scale_color_manual(
    values = c("Pooled OLS" = "green", 
               "Pooled Causal" = "magenta", 
               "Martingale" = "black", 
               "Ancestor Martingale" = "#0000FF",
               "WDRO" = "yellow",  
               "Lasso" = "#71690F"),
    breaks = c("Pooled Causal","WDRO","Pooled OLS","Lasso","Ancestor Martingale","Martingale"),
    labels = c("Pooled Causal","WDRO","Pooled OLS","Lasso","Ancestor Martingale","Martingale")
  ) +
  guides(
    color = guide_legend(override.aes = list(linewidth = rep(1, 6)))
  ) +
  theme(legend.position = c(0.175, 0.775),
        panel.background = element_blank(), 
        panel.grid = element_blank(), 
        legend.key = element_blank(),
        legend.background = element_rect(colour = 'black', fill = 'white', linetype = 'solid'),
        legend.title = element_blank(),
        legend.text = element_text(size = 6),
        panel.border = element_rect(colour = "black", fill=NA, size=1))
#Show the plot
print(p)
#Save the plot
ggsave(filename = "ancestor.pdf",path=path,plot = p, width = 8, height = 5)
