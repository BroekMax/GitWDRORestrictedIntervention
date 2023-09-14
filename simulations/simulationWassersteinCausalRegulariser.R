####################################################################
# Working directory should be the root of the code folder
####################################################################
source("./src/wasserstein_causal_regulariser.R") #Wasserstein causal regulariser
path="./simulations/results/"

pkgs=c("MASS", "pcalg","ggplot2")
for(package in pkgs){
  if(!require(package,character.only = TRUE)) {install.packages(package)}
  library(package, character.only = TRUE)
}

#Function to generate environment data
env_generator=function(n,d,mu,Sigma,myDAG){
    mu_delta=runif(d,0,1)
    #Shift the mean and scale it
    mu_delta=5*c(mu_delta,0)
    #Error matrix for DAG
    eMat=mvrnorm(n,mu=mu_delta,Sigma)
    d.CnormMat=rmvDAG(n,myDAG,errMat = eMat)
    y=d.CnormMat[,11]
    X=d.CnormMat[,-11]
    return(list(y=y,X=X))
  }

#Function to generate test data that is mean-shifted 
shift_test_generator=function(n,d,mu,Sigma,myDAG,shift){
  #Shift the mean by amount shift
  mu_delta=runif(d,shift,shift+1)
  mu_delta=5*c(mu_delta,0)
  #Error matrix for DAG
  eMat=mvrnorm(n,mu=mu_delta,Sigma)
  d.CnormMat=rmvDAG(n,myDAG,errMat = eMat)
  y=d.CnormMat[,11]
  X=d.CnormMat[,-11]
  return(list(y=y,X=X))
}

n=1000
d=10
test_env=10
n_env=3
k=d+1
I=diag(k)
MSE_fixed=rep(0,test_env)
MSE_fixed2=rep(0,test_env)
MSE_pooled=rep(0,test_env)
MSE_wcr=rep(0,test_env)
MSE_cr=rep(0,test_env)
MSE_fixed3=rep(0,test_env)
MSE_wcr_oracle=rep(0,test_env)
MSE_cr_oracle=rep(0,test_env)
nsim=25
maxmatrix=matrix(0,nsim,8)
meanmatrix=matrix(0,nsim,8)

maxshift=11
meanmean=matrix(0,maxshift,8)
maxmean=matrix(0,maxshift,8)
gamma_oracle_tracker=rep(0,nsim)
gamma_oracle_average=rep(0,maxshift)

#Grid-search for causal regulariser
gammas_cr=c(seq(0.01,0.91,by=0.025),seq(0.92,0.99,by=0.01),seq(0.991,0.999,by=0.001),1)
set.seed(2)
for(shift in 1:maxshift){
  for(sim in 1:nsim){
    #Sample DAG in every iteration
    myDAG=randomDAG(n = 11, prob= 0.4, lB = -1, uB = 1)
    #Observational data
    G=matrix(runif((d+1)^2)*2-1, ncol=d+1) 
    Sigma=t(G) %*% G
    Sigma[d+1,-(d+1)]=0
    Sigma[-(d+1),d+1]=0
    #Encourage heterogeneity
    Sigma=10*Sigma/norm(Sigma,"2")
    mu=rep(0,d+1)
    #Error matrix for DAG
    eMat=mvrnorm(n,mu=mu,Sigma)
    d.CnormMat=rmvDAG(n,myDAG,errMat = eMat)
    y=d.CnormMat[,11]
    X=d.CnormMat[,-11]
    
    #Generate data for environments
    X_e =list()
    y_e =list() 
    for(i in 1:n_env){
      env=env_generator(n,10,mu,Sigma,myDAG)
      X_e[[i]] = env$X
      y_e[[i]] = env$y
    }
    #Construct matrix for pooled OLS, pooled Causal, and pooled Martingale
    pooled_X=do.call(rbind,X_e)
    pooled_X=rbind(X,pooled_X)
    pooled_y=do.call(c,y_e)
    pooled_y=c(y,pooled_y)
    
    Xtest_e =list()
    ytest_e =list() 
    for(i in 1:test_env){
      env=shift_test_generator(n,10,mu,Sigma,myDAG,2.5*shift-2.5)
      Xtest_e[[i]] = sweep(env$X,2,colMeans(pooled_X),"-")
      ytest_e[[i]] = env$y-mean(pooled_y)
    }
    
    for(i in 1:n_env){
      X_e[[i]] = sweep(X_e[[i]],2,colMeans(pooled_X),"-")
      y_e[[i]] = y_e[[i]]-mean(pooled_y)
    }

    X=sweep(X,2,colMeans(pooled_X),"-")
    y=y-mean(pooled_y)

    pooled_X=sweep(pooled_X,2,colMeans(pooled_X),"-")
    pooled_y=pooled_y-mean(pooled_y)
    lm(y~X)$coefficients[-1]
    lm(pooled_y~pooled_X)$coefficients[-1]
    coef_ols=lm(y~X)$coefficients[-1]
    coef_pooled=lm(pooled_y~pooled_X)$coefficients[-1]
    
    #Adjacency matrix to find matrix B for Martingale estimator
    Adjacency=matrix(0,k,k)
    B_true=matrix(0,k,k)
    B=matrix(0,k,k)
    o=1
    #Find start and end node to fill adjacency matrix
    for (edge in names(myDAG@edgeData@data)) {
      splitted=strsplit(edge, split = "\\|")[[1]]
      start_node=as.numeric(splitted[1])
      end_node=as.numeric(splitted[2])
      Adjacency[end_node,start_node]=1
      B_true[end_node,start_node]=unlist(myDAG@edgeData@data[[o]])
      o=o+1
    }
    #Regress each variable on its parents to fill B
    for(o in 1:d){
      parent=which(Adjacency[o,]!=0)
      if(sum(Adjacency[o,])>0){
        B[o,parent]=lm(pooled_X[,o]~pooled_X[,parent])$coefficients[-1]
      }
    }
    #Regress on parent of response to fill last row of B
    parent=which(Adjacency[k,]!=0)
    causal=lm(pooled_y~1)
    if(sum(Adjacency[k,])>0){
      causal=lm(pooled_y~pooled_X[,parent])
      B[k,parent]=causal$coefficients[-1]
    }
    #Construct matrix A needed for the Martingale estimator
    A=t(I-B)%*%(I-B)
    
    #CV on 6 environments, find best value on gamma on left out environments
    samp=sample(n_env,floor(0.7*n_env))
    
    X_etrain=X_e[samp]
    X_eval=X_e[-samp]
    y_etrain=y_e[samp]
    y_eval=y_e[-samp]
    
    pooledtrain_X=do.call(rbind,X_etrain)
    pooledtrain_y=do.call(c,y_etrain)
    
    pooledval_X=do.call(rbind,X_eval)
    pooledval_y=do.call(c,y_eval)

    #The value is monotonically increasing beyond a certain value, if value found do a finer search
    low=Inf
    #First search by 5, if the loss increases then the lowest value is around the one with the lowest value
    for(gamma in seq(0,100,by=5)){
      wcr=wasserstein_causal_regulariser(X,y,X_etrain,y_etrain,gamma=gamma,beta_init = coef_pooled,n_iter=50)
      nextvalue=mean((pooledval_X%*%wcr-pooledval_y)^2)
      if(nextvalue<=low){
        gamma_tracker=gamma
        low=nextvalue
      }
    }  
    
    low=Inf
    #Finer search
    for(gamma2 in seq(max(gamma_tracker-5,0),gamma_tracker+5,by=1)){
      wcr=wasserstein_causal_regulariser(X,y,X_etrain,y_etrain,gamma=gamma2,beta_init = coef_pooled,n_iter=50)
      nextvalue=mean((pooledval_X%*%wcr-pooledval_y)^2)
      if(nextvalue<=low){
        gamma_tracker=gamma2
        low=nextvalue
      }
    }
    
    low=Inf
    for(gamma3 in seq(max(gamma_tracker-1,0),gamma_tracker+1,by=0.1)){
      wcr=wasserstein_causal_regulariser(X,y,X_etrain,y_etrain,gamma=gamma3,beta_init = coef_pooled,n_iter=50)
      nextvalue=mean((pooledval_X%*%wcr-pooledval_y)^2)
      if(nextvalue<=low){
        gamma_tracker=gamma3
        low=nextvalue
      }
    }
    
    low=Inf
    #Finest search
    for(gamma4 in seq(max(gamma_tracker-0.1,0.001),gamma_tracker+0.1,by=0.01)){
      wcr=wasserstein_causal_regulariser(X,y,X_etrain,y_etrain,gamma=gamma4,beta_init = coef_pooled,n_iter=50)
      nextvalue=mean((pooledval_X%*%wcr-pooledval_y)^2)
      if(nextvalue<=low){
        gamma_tracker=gamma4
        low=nextvalue
      }
    }
    #Use gamma_tracker to estimate wcr
    beta_wcr=wasserstein_causal_regulariser(X,y,X_e,y_e,gamma=gamma_tracker,n_iter=50)
    
    low=Inf
    #Similar approach for causal regularisation but the value is never greater than 1
    for(gamma in seq(0,1,by=0.1)){
      cr=causal_regulariser(X,y,X_etrain,y_etrain,gamma=gamma)
      nextvalue=mean((pooledval_X%*%cr-pooledval_y)^2)
      if(nextvalue<=low){
        gamma_tracker=gamma
        low=nextvalue
      }
    }
    
    low=Inf
    for(gamma2 in seq(max(gamma_tracker-0.1,0),min(gamma_tracker+0.1,1),by=0.01)){
      cr=causal_regulariser(X,y,X_etrain,y_etrain,gamma=gamma2)
      nextvalue=mean((pooledval_X%*%cr-pooledval_y)^2)
      if(nextvalue<=low){
        gamma_tracker=gamma2
        low=nextvalue
      }
    }
    
    low=Inf
    for(gamma3 in seq(max(gamma_tracker-0.01,0),min(gamma_tracker+0.01,1),by=0.001)){
      cr=causal_regulariser(X,y,X_etrain,y_etrain,gamma=gamma3)
      nextvalue=mean((pooledval_X%*%cr-pooledval_y)^2)
      if(nextvalue<=low){
        gamma_tracker=gamma3
        low=nextvalue
      }
    }
    
    beta_cr=causal_regulariser(X,y,X_e,y_e,gamma=gamma_tracker)
    
    min_value_cr=Inf
    for(gamma in gammas_cr){
      beta_cr_oracle=causal_regulariser(X,y,X_e,y_e,gamma=gamma)
      for(i in 1:test_env){
        MSE_cr_oracle[i]=mean((Xtest_e[[i]]%*%beta_cr_oracle-ytest_e[[i]])^2)
        
      }
      if(mean(MSE_cr_oracle)<min_value_cr){
        min_gamma=gamma
        min_value_cr=mean(MSE_cr_oracle)
      }
    }
    beta_cr_oracle=causal_regulariser(X,y,X_e,y_e,gamma=min_gamma)
    
    #Do the same search but now for the oracle, look at the test MSE
    low=Inf
    for(gamma in seq(0,100,by=5)){
      wcr=wasserstein_causal_regulariser(X,y,X_e,y_e,gamma=gamma,beta_init = coef_pooled,n_iter=50)
      nextvalue=0
      for(i in 1:test_env){
        nextvalue=nextvalue+mean((Xtest_e[[i]]%*%wcr-ytest_e[[i]])^2)
      }
      if(nextvalue<=low){
        gamma_tracker=gamma
        low=nextvalue
      }
    }  
    low=Inf
    for(gamma2 in seq(max(gamma_tracker-5,0),gamma_tracker+5,by=1)){
      wcr=wasserstein_causal_regulariser(X,y,X_e,y_e,gamma=gamma2,beta_init = coef_pooled,n_iter=50)
      nextvalue=0
      for(i in 1:test_env){
        nextvalue=nextvalue+mean((Xtest_e[[i]]%*%wcr-ytest_e[[i]])^2)
      }
      if(nextvalue<=low){
        gamma_tracker=gamma2
        low=nextvalue
      }
    }  
    
    low=Inf
    for(gamma3 in seq(max(gamma_tracker-1,0),gamma_tracker+1,by=0.1)){
      wcr=wasserstein_causal_regulariser(X,y,X_e,y_e,gamma=gamma3,beta_init = coef_pooled,n_iter=50)  
      nextvalue=0
      for(i in 1:test_env){
        nextvalue=nextvalue+mean((Xtest_e[[i]]%*%wcr-ytest_e[[i]])^2)
      }
      if(nextvalue<=low){
        gamma_tracker=gamma3
        low=nextvalue
      }
    }
    #gamma_tracker=10
    low=Inf
    for(gamma4 in seq(max(gamma_tracker-0.1,0.001),gamma_tracker+0.1,by=0.01)){
      wcr=wasserstein_causal_regulariser(X,y,X_e,y_e,gamma=gamma4,beta_init = coef_pooled,n_iter=50)
      nextvalue=0
      for(i in 1:test_env){
        nextvalue=nextvalue+mean((Xtest_e[[i]]%*%wcr-ytest_e[[i]])^2)
      }
      #print(gamma4)
      #print(nextvalue)
      
      if(nextvalue<=low){
        gamma_tracker=gamma4
        low=nextvalue
      }
    }
    beta_wcr_oracle=wasserstein_causal_regulariser(X,y,X_e,y_e,gamma=gamma_tracker,beta_init = coef_pooled,n_iter=50)
    gamma_oracle_tracker[sim]=gamma_tracker
    
    #Tracking gamma=1 and gamma=10 to understand the behaviour
    fixed_wcr=wasserstein_causal_regulariser(X,y,X_e,y_e,beta_init = coef_pooled,gamma=1)
    fixed_wcr2=wasserstein_causal_regulariser(X,y,X_e,y_e,beta_init = coef_pooled,gamma=5)
    fixed_wcr3=wasserstein_causal_regulariser(X,y,X_e,y_e,beta_init = coef_pooled,gamma=10)
    #For each test environment find the MSE
    for(i in 1:test_env){
      
      MSE_cr_oracle[i]=mean((Xtest_e[[i]]%*%beta_cr_oracle-ytest_e[[i]])^2)
      MSE_wcr_oracle[i]=mean((Xtest_e[[i]]%*%beta_wcr_oracle-ytest_e[[i]])^2)
      
      MSE_cr[i]=mean((Xtest_e[[i]]%*%beta_cr-ytest_e[[i]])^2)
      MSE_wcr[i]=mean((Xtest_e[[i]]%*%beta_wcr-ytest_e[[i]])^2)
      MSE_fixed3[i]=mean((Xtest_e[[i]]%*%fixed_wcr3-ytest_e[[i]])^2)
      MSE_pooled[i]=mean((Xtest_e[[i]]%*%coef_pooled+lm(pooled_y~pooled_X)$coefficients[1]-ytest_e[[i]])^2)
      MSE_fixed[i]=mean((Xtest_e[[i]]%*%fixed_wcr-ytest_e[[i]])^2)
      MSE_fixed2[i]=mean((Xtest_e[[i]]%*%fixed_wcr2-ytest_e[[i]])^2)
    }
    #Mean MSE for each test environment
    meanmatrix[sim,1]=mean(MSE_fixed)
    meanmatrix[sim,2]=mean(MSE_pooled)
    meanmatrix[sim,3]=mean(MSE_fixed2)
    meanmatrix[sim,4]=mean(MSE_fixed3)
    meanmatrix[sim,5]=mean(MSE_wcr)
    meanmatrix[sim,6]=mean(MSE_cr)
    meanmatrix[sim,7]=mean(MSE_wcr_oracle)
    meanmatrix[sim,8]=mean(MSE_cr_oracle)
    
    #The worst-case test MSE
    maxmatrix[sim,1]=max(MSE_fixed)
    maxmatrix[sim,2]=max(MSE_pooled)
    maxmatrix[sim,3]=max(MSE_fixed2)
    maxmatrix[sim,4]=max(MSE_fixed3)
    maxmatrix[sim,5]=max(MSE_wcr)
    maxmatrix[sim,6]=max(MSE_cr)
    maxmatrix[sim,7]=max(MSE_wcr_oracle)
    maxmatrix[sim,8]=max(MSE_cr_oracle)
  }
  meanmean[shift,1]=mean(meanmatrix[,1])
  meanmean[shift,2]=mean(meanmatrix[,2])
  meanmean[shift,3]=mean(meanmatrix[,3])
  meanmean[shift,4]=mean(meanmatrix[,4])
  meanmean[shift,5]=mean(meanmatrix[,5])
  meanmean[shift,6]=mean(meanmatrix[,6])
  meanmean[shift,7]=mean(meanmatrix[,7])
  meanmean[shift,8]=mean(meanmatrix[,8])
  
  maxmean[shift,1]=mean(maxmatrix[,1])
  maxmean[shift,2]=mean(maxmatrix[,2])
  maxmean[shift,3]=mean(maxmatrix[,3])
  maxmean[shift,4]=mean(maxmatrix[,4])
  maxmean[shift,5]=mean(maxmatrix[,5])
  maxmean[shift,6]=mean(maxmatrix[,6])
  maxmean[shift,7]=mean(maxmatrix[,7])
  maxmean[shift,8]=mean(maxmatrix[,8])
  
  #Track best gamma
  gamma_oracle_average[shift]=mean(gamma_oracle_tracker)
}

meanmean
maxmean

gamma_oracle_average

#Plotting
smoothedmatrix=meanmean[,-c(1,4)]
shifts=seq(0,25,by=2.5)
for(b in 1:6){
  loess_fit=loess(smoothedmatrix[,b] ~ shifts,span=1.5)
  smoothedmatrix[,b]=predict(loess_fit)
}
result=as.data.frame(smoothedmatrix) 
result$shifter=shifts
p = ggplot(result, aes(x = shifter)) +
  geom_line(aes(y = V1, color = "WCR gamma")) +
  geom_line(aes(y = V2, color = "Pooled OLS")) +
  geom_line(aes(y = V3, color = "WCR CV")) +
  geom_line(aes(y = V4, color = "CR CV")) +
  geom_line(aes(y = V5, color = "Oracle WCR")) +
  geom_line(aes(y = V6, color = "Oracle CR")) +
  theme_minimal() +
  xlab("Perturbation strength") +
  ylab("Test MSE") +
  theme(axis.text.x = element_text(angle = 0, hjust = 1)) +
  scale_x_continuous(breaks = c(5,10,15,20,25))+
  scale_color_manual(
    values = c("Pooled OLS" = "green", 
               "WCR gamma" = "magenta", 
               "WCR CV" = "black", 
               "CR CV" = "#8B9400",
               "Oracle WCR" = "#8B0000",
               "Oracle CR" = "#4169E1"),
    breaks = c("WCR CV", "CR CV", "WCR gamma", "Pooled OLS", "Oracle WCR", "Oracle CR"),
    labels = c("WCR CV", "CR CV", "Pooled OLS", expression(paste("WCR ", gamma, "=5")), "Oracle WCR", "Oracle CR")
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
        panel.border = element_rect(colour = "black", fill=NA, linewidth=1))
print(p)
ggsave(filename = "plot_mean_wcr.pdf", path = path, plot = p, width = 8, height = 5)

#load("max.RData")
smoothedmatrix=maxmean[,-c(1,4)]
shifts=seq(0,25,by=2.5)
for(b in 1:6){
  loess_fit=loess(smoothedmatrix[,b] ~ shifts,span=1.5)
  smoothedmatrix[,b]=predict(loess_fit)
}
result=as.data.frame(smoothedmatrix) 
result$shifter=shifts
p = ggplot(result, aes(x = shifter)) +
  geom_line(aes(y = V1, color = "WCR gamma")) +
  geom_line(aes(y = V2, color = "Pooled OLS")) +
  geom_line(aes(y = V3, color = "WCR CV")) +
  geom_line(aes(y = V4, color = "CR CV")) +
  geom_line(aes(y = V5, color = "Oracle WCR")) +
  geom_line(aes(y = V6, color = "Oracle CR")) +
  theme_minimal() +
  xlab("Perturbation strength") +
  ylab("Worst-case test MSE") +
  theme(axis.text.x = element_text(angle = 0, hjust = 1)) +
  scale_x_continuous(breaks = c(5,10,15,20,25))+
  scale_color_manual(
    values = c("Pooled OLS" = "green", 
               "WCR gamma" = "magenta", 
               "WCR CV" = "black", 
               "CR CV" = "#8B9400",
               "Oracle WCR" = "#8B0000",
               "Oracle CR" = "#4169E1"),
    breaks = c("WCR CV", "CR CV", "WCR gamma", "Pooled OLS", "Oracle WCR", "Oracle CR"),
    labels = c("WCR CV", "CR CV", "Pooled OLS", expression(paste("WCR ", gamma, "=5")), "Oracle WCR", "Oracle CR")
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
        panel.border = element_rect(colour = "black", fill=NA, linewidth=1))

print(p)
ggsave(filename = "plot_max_wcr.pdf", path = path, plot = p, width = 8, height = 5)

#Plotting oracle gamma
smoothed=predict(loess(gamma_oracle_average~shifts,span=1))
smoothed
result=as.data.frame(smoothed)
result$shifter=shifts
p = ggplot(result, aes(x = shifter)) +
  geom_line(aes(y = smoothed, color = "WCR")) +
  theme_minimal() +
  xlab("Perturbation strength") +
  ylab(expression(gamma)) +
  theme(axis.text.x = element_text(angle = 0, hjust = 1)) +
  scale_x_continuous(breaks = c(5,10,15,20,25))+
  scale_color_manual(
    values = c("WCR" = "#8B9400"),
    breaks = c("WCR"),
    labels = c( expression(paste("WCR ", gamma, "-oracle")))
  ) +
  guides(
    color = guide_legend(override.aes = list(linewidth = rep(1, 1)))
  ) +
  theme(legend.position = c(0.175, 0.775),
        panel.background = element_blank(), 
        panel.grid = element_blank(), 
        legend.key = element_blank(),
        legend.background = element_rect(colour = 'black', fill = 'white', linetype = 'solid'),
        legend.title = element_blank(),
        legend.text = element_text(size = 6),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=1))
print(p)
ggsave(filename = "gamma_oracle.pdf", path = path, plot = p, width = 8, height = 5)
