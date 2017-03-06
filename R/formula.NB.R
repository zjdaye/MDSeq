# formula for NB distribution

# Negative binomial distribution
loglik.MD <- function(beta, gamma, y, X, U, logoffsets, max_val=1e300, tol=1e-12) {
  yeq0 = (y==0); ygt0 = !yeq0
  mu  = pmin( exp(X%*%beta  + ygt0*logoffsets), max_val); 
  phi = pmax( exp(U%*%gamma + ygt0*logoffsets), 1+tol); 
  inv_phi1 = 1/(phi-1); theta = mu*inv_phi1; logphi = log(phi)
  return( sum(lgamma(y+theta) - lgamma(theta) - lgamma(y+1) - y*(logphi+log(inv_phi1)) - theta*logphi) )
}

# Closed form gradient for loglik of nb glm
loglik.MD.grad <- function(beta, gamma, y, X, U, logoffsets, max_val=1e300, tol=1e-12) {
  p = dim(X)[2]; q = dim(U)[2]
  yeq0 = (y==0); ygt0 = !yeq0
  mu  = pmin( exp(X%*%beta  + ygt0*logoffsets), max_val); 
  phi = pmax( exp(U%*%gamma + ygt0*logoffsets), 1+tol); 
  
  # Temporary variables for efficient computation
  phi1 = (phi-1); inv_phi1 = 1/phi1; theta = mu*inv_phi1
  A = digamma(y+theta) - digamma(theta) - log(phi)
  
  Li.beta = theta * A; grad.beta = t(X) %*% Li.beta
  Li.gamma = inv_phi1 * ( y - phi*Li.beta - theta*phi1 ); grad.gamma = t(U) %*% Li.gamma
  
  return(list(grad.beta=grad.beta, grad.gamma=grad.gamma))
}

# Closed form hessian for loglik of nb glm
loglik.MD.hessian <- function(beta, gamma, y, X, U, logoffsets, max_val=1e300, tol=1e-12) 
{
  p = dim(X)[2]; q = dim(U)[2]; n = length(y)
  yeq0 = (y==0); ygt0 = !yeq0
  mu  = pmin( exp(X%*%beta  + ygt0*logoffsets), max_val); 
  phi = pmax( exp(U%*%gamma + ygt0*logoffsets), 1+tol); 
  
  # Temporary variables for efficient computation
  phi1 = phi-1; inv.phi1 = 1/phi1; theta = mu*inv.phi1
  phi_inv.phi1 = phi*inv.phi1; phi_inv.phi1_2 = phi*(inv.phi1^2)
  A = digamma(y+theta) - digamma(theta) - log(phi)
  B = trigamma(y+theta) - trigamma(theta)
  
  beta_multiplier = theta*(A+theta*B)
  gamma_multiplier = phi_inv.phi1_2*( theta*(2*phi1 + A) + phi*beta_multiplier - y )
  joint_multiplier = - (theta + phi_inv.phi1*beta_multiplier)
  
  hessian_beta = t(X)%*%(X*array(beta_multiplier,c(n,p)))
  hessian_gamma = t(U)%*%(U*array(gamma_multiplier,c(n,q)))
  hessian_joint = t(X)%*%(U*array(joint_multiplier,c(n,q)))
  
  hessian = rbind(cbind(hessian_beta, hessian_joint), cbind(t(hessian_joint), hessian_gamma))
  rownames(hessian) = colnames(hessian) = NULL
  
  return(hessian)
}


##################################################################
# helper for constrained 

loglik.MD_helper <- function(pars, y, X, U, logoffsets, max_val=1e300, tol=1e-12) {
  p = dim(X)[2]; q = dim(U)[2]
  return ( -loglik.MD(beta=pars[1:p], gamma=pars[(p+1):(p+q)], y=y, X=X, U=U, logoffsets=logoffsets, max_val=max_val, tol=tol) )
}

loglik.MD.grad_helper <- function(pars, y, X, U, logoffsets, max_val=1e300, tol=1e-12) {
  p = dim(X)[2]; q = dim(U)[2]
  obj = loglik.MD.grad(beta=pars[1:p], gamma=pars[(p+1):(p+q)], y=y, X=X, U=U, logoffsets=logoffsets, max_val=max_val, tol=tol)
  return( -c(obj$grad.beta, obj$grad.gamma) )
}