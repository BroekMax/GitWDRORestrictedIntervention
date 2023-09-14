####################################################################
# Code that has the functions for the Iterative scheme
####################################################################

#Helper function that is an easier projection
pi_W=function(beta, lambda, K1, K2, R_beta){
  if (is_W(beta, lambda)) {
    return(list(beta = beta, lambda = lambda))
  }
  #Check if in the projection space
  else{if (lambda < 0.5 | lambda>20) {
    #Project to boundary
    lambda[lambda>20]=20
    lambda[lambda<0.5]=0.5
  }
    #Project to boundary
    beta[beta > 5]=5
    beta[beta < -5]=-5
    return(list(beta = beta, lambda = lambda))}
}

#Helper function that projects the parameters from gradient descent to the projection space
pi_W=function(beta, lambda, K1, K2, R_beta) {
  #Compute the norm of the parameter
  beta_norm=sqrt(sum(beta[1:d]^2))
  #Case1: check if in feasible space
  if (is_W(beta, lambda)) {
    return(list(beta = beta, lambda = lambda))
  }
  #Case 2:
  if (beta_norm <= K2 * R_beta / K1 && lambda > K2 * R_beta) {
    return(list(beta = beta, lambda = K2 * R_beta))
  }
  
  #Case 3:
  if (lambda < -beta_norm / K1) {
    return(list(beta = c(rep(0,k-1),-1), lambda = 0))
  }
  
  #Case 4:
  if (-beta_norm / K1 <= lambda && lambda < min(K1 * beta_norm, K2 * R_beta * (1 + 1/K1^2) - beta_norm / K1)) {
    return(list(beta = as.vector((beta_norm + lambda) / (1 + K2)) * beta / beta_norm , 
                lambda = (K1 * beta_norm + K2 * lambda) / (1 + K2)))
  }
  
  #Case 5 (otherwise):
  return(list(beta = K2 * R_beta / K1 * (beta / beta_norm), lambda = K2 * R_beta))
}

#Helper function checks if in the projection space for the projected gradient descent step
is_W=function(beta, lambda) {
  if(!any(abs(beta) > 5) && lambda > 0.5){
    return(TRUE)
  }
  return(FALSE)
}

#Helper function to compute the gradient of the parameter for the gradient scheme
gradient_beta=function(x, beta, lambda, epsilon, A,batch) {
  gradient=numeric(length(beta))
  #Make distinction between gradient descent and stochastic gradient descent (batch=1)
  for(i in 1:batch){
    if(batch>1){
      xi=x[i,]
    }
    else{
      xi=x
    }
    v=lambda - sqrt(epsilon) * t(beta) %*% solve(A) %*% beta
    numerator=2 * as.vector(lambda) * as.vector(as.numeric(t(beta) %*% xi)) * as.vector(xi) * as.numeric(v)
    denominator=v^2
    grad=(numerator - as.numeric(lambda * (t(beta) %*% xi)^2) * (-2 * sqrt(epsilon) * solve(A) %*% beta)) / as.numeric(denominator)
    #Add the gradient of all samples in the batch
    gradient=gradient+grad
  }
  #Average the gradient over the batch
  return(gradient/batch)
}

#Helper function to compute the gradient of the multiplier in the dual
gradient_lambda=function(x, beta, lambda, epsilon, A,batch) {
  gradient=0
  #Find average gradient over batch
  for(i in 1:batch){
    if(batch>1){
      xi=x[i,]
    }
    else{
      xi=x
    }
    v=lambda-sqrt(epsilon)*t(beta)%*%solve(A)%*%beta
    grad=sqrt(epsilon) + ((t(beta) %*% xi)^2/v)-(lambda*(t(beta)%*%xi)^2) / v^2
    #Add the gradient of all samples in the batch
    gradient=grad+gradient
  }
  #Average the gradient over the batch
  return(gradient/batch)
}

####################################################################
# Function computing the iterative scheme estimator
# Inputs:
# X: Design matrix
# y: Response vector
# A: Matrix t(I-B)%*%(I-B)
# n_iter: Maximum number of iterations if not converged
# burn_in: Number of iterations before main scheme, can affect trajectory
# beta_init: Initial value for estimator
# lambda_init: Initial value for Lagrange multiplier
# alpha: Parameter affecting the learning rate, needs to be positive
# tau: Parameter affecting the learning rate, should be in (0.5,1)
# epsilon: Radius ambiguity set, functions as regularisation parameter
# Output:
# last_beta: The last value of the GD
# lastlambda: Last value of GD for multiplier
# beta: Estimator, last value of the trajectory
# lambda: Estimate for Lagrange multiplier, last value of trajectory
####################################################################
#Iterative scheme estimator
iterative_scheme=function(X,y,A,n_iter,burn_in=10,beta_init=NULL,lambda_init=10,alpha=0.55,tau=0.9,epsilon=0.05){
  #Values used in the projection
  K1=1
  R_beta=8 
  K2=1
  d=ncol(X)
  n=nrow(X)
  k=d+1
  xy=cbind(X,y)
  converged=FALSE
  batch=n
  
  #Initialising lambda and beta
  if (is.null(beta_init)) {
    beta <- lm(xy[,k]~xy[,-k])$coefficients[-1]
  } else {
    beta <- beta_init
  }
  if (is.null(lambda_init)) {
    lambda <- 5
  } else {
    lambda <- lambda_init
  }
  beta=c(beta,-1)
  grad_beta=numeric(k)
  grad_lambda=0
  traj_average_beta=matrix(0,(n_iter-burn_in),k)
  traj_average_lambda=rep(0,(n_iter-burn_in))
  #Burn_in if initial value might affect the trajectory too much
  if(burn_in>0){
    for(i in 1:burn_in){
      #Decreasing learning rate
      alpha_i=alpha*i^(-tau)
      #Sample batch from data
      samp=sample(n,batch)
      xi=xy[samp,]
      #Calculate gradients
      grad_beta=gradient_beta(xi,beta,lambda,epsilon,A,batch)
      grad_beta[k]=0
      grad_lambda=gradient_lambda(xi,beta,lambda,epsilon,A,batch)
      #Gradient descent
      beta=beta-alpha_i*grad_beta
      lambda=lambda-alpha_i*grad_lambda
      #Projection
      pi=pi_W(beta,lambda,K1,K2,R_beta)
      beta=pi$beta
      #Project value of response to -1
      beta[k]=-1
      lambda=pi$lambda
    }
  }
  
  j=0
  #Iterate till convergence or max number of iterations
  while(converged==FALSE && j<n_iter){
    #iteration counter, j counts only the iterations without burn_in
    j=j+1
    #Previous values needed for convergence criterion
    beta_prev=beta
    lambda_prev=lambda
    #Decreasing learning rate
    alpha_j=alpha*j^(-tau)
    samp=sample(n,batch)
    xi=xy[samp,]
    #Gradients
    grad_beta=gradient_beta(xi,beta,lambda,epsilon,A,batch)
    grad_beta[k]=0
    grad_lambda=gradient_lambda(xi,beta,lambda,epsilon,A,batch)
    #Gradient descent
    beta=beta-alpha_j*grad_beta
    lambda=lambda-alpha_j*grad_lambda
    #Projection
    pi=pi_W(beta,lambda,K1,K2,R_beta)
    beta=pi$beta
    beta[k]=-1
    lambda=pi$lambda
    traj_average_beta[j,]=(j-1)/j*as.vector(traj_average_beta[max(j-1,1),])+1/j*beta
    traj_average_lambda[j]=(j-1)/j*traj_average_lambda[max(j-1,1)]+1/j*lambda
    #Change in iterations to determine convergence
    beta_change = sum(abs(beta - beta_prev))
    lambda_change = abs(lambda - lambda_prev)
    
    #Convergence criterion
    converged = (beta_change < 0.00001) & (lambda_change < 0.01) 
  }
  #Take last value of trajectory
  return(list(lastbeta = beta, lastlambda = lambda,beta=traj_average_beta[j,],lambda=traj_average_lambda[j]))
}