####################################################################
# Working directory should be the root of the code folder
####################################################################
pkgs=c("MASS","ggplot2")
for(package in pkgs){
  if(!require(package,character.only = TRUE)) {install.packages(package)}
  library(package, character.only = TRUE)
}
source("./src/martingaleWDRO.R") #Martingale estimator
path="./simulations/results/"


####################################################################
# Function that generates plots for the Martingale estimator
# Input: seed
# Output: ggplot2 graph
####################################################################
plot_generator=function(seed=302){
  set.seed(seed)
  n=100
  ntest=1000
  d=10
  k=d+1
  B = matrix(0, k, k)
  #The coefficients of the covariates into the response
  B[k,1:10]=0
  B[k,9]=1
  B[k,10]=2
  I=diag(k)
  A=t(I-B)%*%(I-B)
  #Finding the causal parents
  pa=which(B[k,-k]!=0)
  #Half of the data is from the interventional distribution
  frac=0.5
  #Use setting 4, where a shift is present in the test distribution
  X=matrix(rnorm(n * d,0,1), n, d)
  Xtest=matrix((rnorm(ntest * d,0,1)+rnorm(ntest*d,sqrt(0.5),1)), ntest, d) 
  X[ceiling(frac*n):n,9]=X[ceiling(frac*n):n,9]+rnorm(n-ceiling(frac*n)+1,0,5)
  
  #Normally distributed confounder
  H=rnorm(n,0,1)
  Htest=rnorm(ntest,0,1)
  #Add confounder
  X=X+H
  Xtest=Xtest+Htest
  
  #If propagation then can generate X via this loop:
  for(j in 1:d){
    X[,j]=B[j,1:d]%*%t(X)+X[,j]     
  }
  for(j in 1:d){
    Xtest[,j]=B[j,1:d]%*%t(Xtest)+Xtest[,j]    
  }
  
  #Generate response by 1*X9+2*X10 and a confounder
  y=as.numeric(B[k,-k]%*%t(X))+rnorm(n,0,1)+H
  ytest=as.numeric(B[k,-k]%*%t(Xtest))+rnorm(ntest,0,1)+Htest
  Xtest=sweep(Xtest,2,colMeans(X),FUN="-")
  
  #Standardise
  X=sweep(X,2,colMeans(X),FUN="-")
  ytest=ytest-mean(y)
  y=y-mean(y)
  
  B_2= matrix(0, k, k)
  B_22= matrix(0, k, k)
  causal=lm(y~X[,pa])
  #Construct B matrix with estimated causal coefficients
  B_2[k,pa]=causal$coefficients[-1]
  #Construct A that is used to solve the estimator
  A_2=t(I-B_2)%*%(I-B_2)
  
  #Values for the plot
  epsilons=seq(0.001,2500,by=0.4)
  
  #Save the MSE over different values of the radius (here named lambda)
  MSE_saver=rep(0,length(epsilons))
  coeffA=matrix(0,nrow=length(epsilons),ncol=d)
  i=1
  #MSE for each value of epsilon to construct the plot
  #coeffA can be used for the regularisation path
  for(epsilon in epsilons){
    beta_estimated_M2=tikhonov_constrained(X,y,A_2,epsilon)
    coeffA[i,]=beta_estimated_M2
    ypred_M2=Xtest%*%beta_estimated_M2
    MSE_saver[i]=mean((ypred_M2-ytest)^2)
    i=i+1
  }
  
  #Put results into dataframe to be able to use ggplot
  result=data.frame(Epsilon=epsilons,coeffA=coeffA)
  #Make the plot
  p=ggplot(result, aes(x = Epsilon)) +
    geom_line(aes(y=coeffA[,10], color = "A1")) +
    geom_line(aes(y=coeffA[,9], color = "A2")) +
    geom_line(aes(y=coeffA[,1], color = "A3")) +
    geom_line(aes(y=coeffA[,2], color = "A3")) +
    geom_line(aes(y=coeffA[,3], color = "A3")) +
    geom_line(aes(y=coeffA[,4], color = "A3")) +
    geom_line(aes(y=coeffA[,5], color = "A3")) +
    geom_line(aes(y=coeffA[,6], color = "A3")) +
    geom_line(aes(y=coeffA[,7], color = "A3")) +
    geom_line(aes(y=coeffA[,8], color = "A3")) +
    geom_vline(xintercept = which(MSE_saver==min(MSE_saver)) , 
               color = "blue")+
    geom_hline(yintercept = 2, linetype = "dashed", color = "red") +
    geom_hline(yintercept = 1, linetype = "dashed", color = "black") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "green") +
    xlab(expression(epsilon)) +
    ylab("Coefficient") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 0, hjust = 1)) +
    scale_color_manual(values = c("A1" = "red", "A2" = "black", "A3" = "green"))+
    theme(legend.position = "none")+
    theme(
      panel.grid = element_blank(), 
      panel.border = element_rect(colour = "black", fill=NA, linewidth=1))
  return(p)
}
#Regularisation path
plot1=plot_generator(2)
ggsave(filename ="regupath.pdf",path=path,plot=plot1,width=8,height=5)