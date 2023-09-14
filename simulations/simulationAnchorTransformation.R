####################################################################
# Working directory should be the root of the code folder
####################################################################
source("./src/ph_anchor.R") #Anchor estimators
pkgs=c("ggplot2","MASS")
for(package in pkgs){
  if(!require(package,character.only = TRUE)) {install.packages(package)}
  library(package, character.only = TRUE)
}
####################################################################
#This is a simulation from Section 6.1 with the anchor transformation
####################################################################
#Do nsim iterations for each value of delta
simulation_delta=function(nsim,verbose=TRUE){
  set.seed(421)
  #Number of test points
  ntest=100
  #Sample size
  n=100
  #Dimension
  d=10
  #Identity matrix for anchor transformation
  Id=diag(1,n)
  gamma=5
  #Grid for values of delta
  deltas=seq(0.05,2.25,0.1)
  #Regularisation paremeter grid for Anchor regression
  lambdas=seq(0.001,0.2,0.005)
  MSE_test_anch_vec=rep(0,nsim)
  MSE_test_anch_ph_vec=matrix(0,nsim,length(deltas))
  #Fixed value for epsilon in PH-Anchor
  epsilon=0.1

  MSE_matrix=rep(0,length(deltas)+1)
  for(sim in 1:nsim){
    #Simulation generation, see paper for details about settings
    A1=rnorm(n,0,1)
    A2=rnorm(n,0,1)
    Atest1=rnorm(ntest,5,1)
    Atest2=rnorm(ntest,5,1)
    M1=rnorm(d)
    M1=(M1-mean(M1))/sd(M1)
    M2=rnorm(d)
    M2=(M2-mean(M2))/sd(M2)
    X=matrix(rnorm(n*d),nrow=n)
    Xtest=matrix(rnorm(ntest*d),nrow=ntest)
    for(i in 1:d){
      X[,i]=X[,i]+M1[i]*A1+M2[i]*A2
      Xtest[,i]=Xtest[,i]+M1[i]*Atest1+M2[i]*Atest2
      Xtest[,i]=Xtest[,i]-mean(X[,i])
      X[,i]=X[,i]-mean(X[,i])
    }
    y=X[,1]+2*X[,2]+rnorm(n)
    ytest=Xtest[,1]+2*Xtest[,2]+rnorm(ntest)
    
    
    #Standardising
    Xtest=sweep(Xtest,2,colMeans(X),FUN="-")
    X=sweep(X,2,colMeans(X),FUN="-")
    ytest=ytest-mean(y)
    y=y-mean(y)
    
    #Transforming the data via the Anchor Transformation
    A=cbind(A1,A2)
    projA=A%*%solve(t(A)%*%A)%*%t(A)
    yanch=(Id+(sqrt(gamma)-1)*projA)%*%y
    Xanch=(Id+(sqrt(gamma)-1)*projA)%*%X
    yanch=as.numeric(unlist(yanch))
    
    #Computing the estimators and the MSE on the test set
    #Anchor
    anchor = anchor_regression(yanch,Xanch,lambdas = lambdas)
    a0_anch = anchor$a0
    coef_anch = anchor$beta
    ypred_anch = (Xtest %*% coef_anch + a0_anch)
    MSE_test_anch = mean((ypred_anch - ytest)^2)
    for(i in 1:length(deltas)){
      #PH Anchor without RWP
      ph_anchor=pseudo_huber_anchor(yanch,Xanch,delta=deltas[i],epsilon)
      coef_anch_pha=ph_anchor$coef
      a0_anch_pha=ph_anchor$a0
      ypred_anch_pha=((Xtest%*%coef_anch_pha)+a0_anch_pha)
      MSE_test_anch_pha=mean((ypred_anch_pha-ytest)^2)
      
      
      #Recording the MSE
      MSE_test_anch_ph_vec[sim,i]=MSE_test_anch_pha
    }
    MSE_test_anch_vec[sim]=MSE_test_anch
    if(verbose==TRUE){
    cat("Iteration:", sim, "\n") 
    }
    
  }
  MSE_matrix[1:length(deltas)]=apply(MSE_test_anch_ph_vec,2,mean)
  MSE_matrix[24]=mean(MSE_test_anch_vec)
  
  return(MSE_matrix)
}

delta_sim_matrix=simulation_delta(50)

#Plotting:
deltas=seq(0.05,2.25,0.1)
result=data.frame(Epsilon=deltas,Test_MSE=delta_sim_matrix[1:23])
p=ggplot(result, aes(x = Epsilon, y = Test_MSE)) +
  geom_line(aes(color = "PH-Anchor"), linewidth = 1) +
  geom_hline(aes(yintercept = delta_sim_matrix[24], color = "Anchor"), linetype = "dashed", linewidth = 0.5) +
  xlab(expression(delta)) +
  ylab("Test MSE") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, hjust = 1)) +
  guides(
    color=guide_legend(override.aes=list(linetype=c("dashed","solid"),linewidth=c(0.02,1))))+
scale_color_manual(
  values = c("Anchor" = "red", "PH-Anchor" = "black"),
  labels = c("Anchor", "PH-Anchor"),
  c("Anchor", "PH-Anchor")) +
  theme(legend.position = c(0.92, 0.92),
        panel.background = element_blank(), 
        panel.grid = element_blank(), 
        legend.key = element_blank(),
        legend.background = element_rect(colour = 'black', fill = 'white', linetype = 'solid'),
        legend.title = element_blank(),
        legend.text = element_text(size = 6),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=1))
print(p)
#Saving plot
ggsave(filename = "anchor_transformation.pdf", path = path, plot = p, width = 8, height = 5)