####################################################################
# Working directory should be the root of the code folder
####################################################################
source("./src/ph_anchor.R") #Anchor estimators

####################################################################
# This is a simulation from Section 6.1 with PH-Anchor-RWP
####################################################################

simulation=function(nsim,verbose=FALSE){
  set.seed(417)
  #Number of test points
  ntest=100
  #Sample size
  n=100
  #Dimension
  d=10
  #Identity matrix for anchor transformation
  Id=diag(n)
  gamma=5
  #Vectors to store the MSE
  MSE_test_h_vec=rep(0,nsim)
  MSE_test_anch_vec=rep(0,nsim)
  MSE_test_anch_ph_vec=rep(0,nsim)
  MSE_test_anch_pha_rwp_vec=rep(0,nsim)
  MSE_test_l_vec=rep(0,nsim)
  
  #Track the regularisation parameters
  lambdas_min=rep(0,nsim)
  epsilons_ancha=rep(0,nsim)
  
  #Grid for cross-validation
  lambdas=seq(0.001,0.2,0.005)
  
  #Delta for huber
  delta_h = 0.01
  #Delta for PH
  delta = 0.01
  #Epsilon for PH-Anchor
  epsilon=0.1
  #Setting 7 is included but not mentioned in the thesis
  nsettings=10
  
  MSE_matrix=matrix(0,nsettings,5)
  Id=diag(1,n)
  for(setting in 1:nsettings){
    #High-dimensional setting has different parameters
    if(setting > 6) {
      delta=0.0125
      epsilon=0.05
      d=110
      lambdas=seq(0.01,1,0.025)
    }
    for(sim in 1:nsim){
      #Simulation generation, see paper for details about settings
      A1=rnorm(n,0,1)
      A2=rnorm(n,0,1)
      Atest1=rnorm(ntest,5,1)
      Atest2=rnorm(ntest,5,1)
      H=rnorm(n)+A2*(setting==8)
      Htest=rnorm(ntest)+Atest2*(setting==8)
      if(setting==2){
        Atest1=rnorm(ntest,5,4)
        Atest2=rnorm(ntest,5,4)
      }
      M1=rnorm(d)
      M1=(M1-mean(M1))/sd(M1)
      M2=rnorm(d)
      M2=(M2-mean(M2))/sd(M2)
      X=matrix(rnorm(n*d),nrow=n)
      Xtest=matrix(rnorm(ntest*d),nrow=ntest)
      for(i in 1:d){
        X[,i]=X[,i]+M1[i]*A1+M2[i]*A2+H*(setting==3|setting==6|setting==8)
        Xtest[,i]=Xtest[,i]+M1[i]*Atest1+M2[i]*Atest2+Htest*(setting==3|setting==6|setting==8)
        Xtest[,i]=Xtest[,i]-mean(X[,i])
        X[,i]=X[,i]-mean(X[,i])
      }
      y=X[,1]+2*X[,2]+rnorm(n)+H*(setting==3|setting==6|setting==8)
      ytest=Xtest[,1]+2*Xtest[,2]+rnorm(ntest)+Htest*(setting==3|setting==6|setting==8)
      if(setting==5|setting==6|setting==9){
        y=y+A1
        ytest=ytest+Atest1
      }
      #Outliers from setting 4 onwards
      if(setting>=4){
        y[7]=y[7]+10
      }
      if(setting==10){
        y[9]=y[9]-5
        y[11]=y[11]+10
      }
      
      
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
      
      #Computing the estimators and the MSE on the test data
      #Anchor
      anchor = anchor_regression(yanch,Xanch,lambdas = lambdas)
      a0_anch = anchor$a0
      coef_anch = anchor$beta
      ypred_anch = (Xtest %*% coef_anch + a0_anch)
      MSE_test_anch = mean((ypred_anch - ytest)^2)
      
      #Huber
      huber=huber_anchor(yanch,Xanch,delta_h)
      coef_h=huber$coef_h
      a0_h=huber$a0_h
      ypred_h=(Xtest%*%coef_h+a0_h)
      MSE_test_h=mean((ypred_h-ytest)^2)
      
      #PH Anchor with RWP 
      ph_anchor_rwp=pseudo_huber_anchor_RWP(yanch,Xanch,delta)
      coef_anch_pha_rwp=ph_anchor_rwp$coef
      a0_anch_pha_rwp=ph_anchor_rwp$a0
      ypred_anch_pha_rwp=((Xtest%*%coef_anch_pha_rwp)+a0_anch_pha_rwp)
      MSE_test_anch_pha_rwp=mean((ypred_anch_pha_rwp-ytest)^2)
      
      #PH Anchor without RWP
      ph_anchor=pseudo_huber_anchor(yanch,Xanch,delta=1.345,epsilon)
      coef_anch_pha=ph_anchor$coef
      a0_anch_pha=ph_anchor$a0
      ypred_anch_pha=((Xtest%*%coef_anch_pha)+a0_anch_pha)
      MSE_test_anch_pha=mean((ypred_anch_pha-ytest)^2)
      
      #Sqrt Lasso with RWP choice of Blanchet et al.
      sqrt_lasso_rwp=sqrt_lasso_RWP(yanch,Xanch)
      coef_l=sqrt_lasso_rwp$beta
      a0_l=sqrt_lasso_rwp$a0
      ypred_l=((Xtest%*%coef_l)+a0_l)
      MSE_test_l=mean((ypred_l-ytest)^2)
      
      #Recording the MSE
      MSE_test_anch_vec[sim]=MSE_test_anch
      MSE_test_h_vec[sim]=MSE_test_h
      MSE_test_anch_ph_vec[sim]=MSE_test_anch_pha
      MSE_test_anch_pha_rwp_vec[sim]=MSE_test_anch_pha_rwp
      MSE_test_l_vec[sim]=MSE_test_l

      if(verbose==TRUE){
        if(sim%%10==0){
          cat("Setting:", setting, ", iteration:", sim, "\n") #Iteration number progress
        }
      }
    }
    #Averaging over each setting
    MSE_matrix[setting,1] = mean(MSE_test_anch_pha_rwp_vec)
    MSE_matrix[setting,2] = mean(MSE_test_anch_ph_vec)
    MSE_matrix[setting,3] = mean(MSE_test_h_vec)
    MSE_matrix[setting,4] = mean(MSE_test_l_vec)
    MSE_matrix[setting,5] = mean(MSE_test_anch_vec)
  }
  return(MSE_matrix)
}

result=simulation(100)
print(result)