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
  #Add a redundant covariate to the parents, here covariate 1
  pa2=c(1,pa)
  B_22[k,pa2]=lm(y~X[,pa2])$coefficients[-1]
  A_22=t(I-B_22)%*%(I-B_22)
  
  #Estimate OLS
  ols=lm(y~X)
  beta_estimated_ols=ols$coefficients[-1]
  ypred_ols=Xtest%*%beta_estimated_ols
  MSE_ols=(mean((ypred_ols-ytest)^2))
  
  #Estimate causal estimator
  beta_estimated_causal=causal$coefficients[-1]
  ypred_causal=Xtest[,pa]%*%beta_estimated_causal
  MSE_causal=(mean((ypred_causal-ytest)^2))
  
  #Values for the plot
  epsilons=seq(0.001,5000,by=0.4)
  
  #Save the MSE over different values of the radius (here named lambda)
  MSE_saver=rep(0,length(epsilons))
  MSE_saver2=rep(0,length(epsilons))
  coeffA=matrix(0,nrow=length(epsilons),ncol=d)
  i=1
  #MSE for each value of epsilon to construct the plot
  #coeffA can be used for the regularisation path
  for(epsilon in epsilons){
    beta_estimated_M2=tikhonov_constrained(X,y,A_2,epsilon)
    coeffA[i,]=beta_estimated_M2
    ypred_M2=Xtest%*%beta_estimated_M2
    beta_estimated_M22=tikhonov_constrained(X,y,A_22,epsilon)
    ypred_M22=Xtest%*%beta_estimated_M22
    MSE_saver[i]=mean((ypred_M2-ytest)^2)
    MSE_saver2[i]=mean((ypred_M22-ytest)^2)
    i=i+1
  }
  
  #Put results into dataframe to be able to use ggplot
  result=data.frame(Epsilon=epsilons,Test_MSE=MSE_saver,Test_MSE2=MSE_saver2)
  #Make the plot
  p=ggplot(result, aes(x = Epsilon)) +
    geom_line(aes(y=Test_MSE, color = "Martingale MSE "), linewidth = 1) +
    geom_line(aes(y=Test_MSE2, color = "Martingale MSE2 "), linewidth = 1) +
    geom_hline(aes(yintercept = MSE_causal, color = "Causal MSE "), linetype = "dashed", linewidth = 1) +
    geom_hline(aes(yintercept = MSE_ols, color = "OLS MSE "), linetype = "dashed", linewidth = 1) +
    xlab(expression(epsilon)) +
    ylab("Test MSE") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 0, hjust = 1)) +
    scale_color_manual(
      values = c("Causal MSE " = "red", "Martingale MSE " = "black", "Martingale MSE2 " = "blue", "OLS MSE "="green"),
      labels = c("Causal estimator","Martingale parent oracle", "Martingale redundant covariate","OLS")
    ) +
    guides(
      color=guide_legend(override.aes=list(linetype=c("dashed","solid","solid","dashed"), linewidth=c(0.5,1,1,0.5)))
    ) +
    theme(legend.position=c(0.865,0.165),panel.background = element_blank(), 
          panel.grid = element_blank(), 
          legend.key=element_blank(),
          legend.background =element_rect(colour = 'black', fill = 'white', linetype = 'solid'),
          legend.title = element_blank(),
          legend.text = element_text(size = 6),
          panel.border = element_rect(colour = "black", fill=NA, linewidth=1))
  
  
  return(p)
}
#Plots from the Thesis
plot1=plot_generator(32)
plot2=plot_generator(35013)
plot3=plot_generator(35001)
plot4=plot_generator(3409)
plot5=plot_generator(330)
plot6=plot_generator(37003)

ggsave(filename ="plot1.pdf",path=path,plot=plot1,width=8,height=5)
ggsave(filename ="plot2.pdf",path=path,plot=plot2,width=8,height=5)
ggsave(filename ="plot3.pdf",path=path,plot=plot3,width=8,height=5)
ggsave(filename ="plot4.pdf",path=path,plot=plot4,width=8,height=5)
ggsave(filename ="plot5.pdf",path=path,plot=plot5,width=8,height=5)
ggsave(filename ="plot6.pdf",path=path,plot=plot6,width=8,height=5)
