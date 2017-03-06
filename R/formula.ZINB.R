# formula for ZINB

# Zero-inflated Negative binomial distribution
loglik.ZIMD <- function(beta, gamma, s, y, X, U, logoffsets, max_val=1e300, tol=1e-12) 
{
  yeq0 = (y==0); ygt0 = !yeq0
  mu  = pmin( exp(X%*%beta  + ygt0*logoffsets), max_val); 
  phi = pmax( exp(U%*%gamma + ygt0*logoffsets), 1+tol); 
  
  # temporal variable that helps caculation
  theta <- mu/(phi-1)
  
  inv.phi = 1/phi; inv.phi.theta = inv.phi^theta
  
  # y == 0
  l1 <- sum(log( s + (1 - s)*(inv.phi.theta[yeq0])))
  # y > 0
  l2 <- sum(ygt0)*log(1 - s) + sum( lgamma(y[ygt0] + theta[ygt0]) - lgamma(theta[ygt0]) - lgamma(y[ygt0] + 1) + log(inv.phi.theta[ygt0]) + y[ygt0]*log(1 - inv.phi[ygt0]))
  
  return(l1+l2)  
}

loglik.ZIMD_helper <- function(x, y, X, U, logoffsets, ...) {
  p = dim(X)[2]; q = dim(U)[2]
  return( loglik.ZIMD(beta=x[1:p], gamma=x[(p+1):(p+q)], s=x[p+q+1], y=y, X=X, U=U, logoffsets=logoffsets, ...) )
}

# Closed form hessian for loglik of nb glm
loglik.ZIMD.hessian <- function(beta, gamma, s, y, X, U, logoffsets, max_val=1e300, tol=1e-12) 
{
  p <- dim(X)[2]; q <- dim(U)[2]; n <- length(y)
  yeq0 = (y==0); ygt0 = !yeq0
  mu  = pmin( exp(X%*%beta  + ygt0*logoffsets), max_val); 
  phi = pmax( exp(U%*%gamma + ygt0*logoffsets), 1+tol); 
  
  # Temporary variables for efficient computation
  yeq0 = as.numeric(y==0); yneq0 = 1-yeq0
  phi1 = phi-1; inv_phi1=1/phi1; inv_phi1_2 = inv_phi1^2; theta <- mu*inv_phi1
  one_s = 1-s; logphi = log(phi); phi_logphi = phi*logphi; 
  phi_to_neg_theta = phi^(-theta); one_phi_to_neg_theta = 1 - phi_to_neg_theta
  ytheta = y+theta; A = digamma(ytheta)-digamma(theta)-logphi; B = trigamma(ytheta)-trigamma(theta)
  invC = 1/(s+one_s*phi_to_neg_theta); yeq0_invC = yeq0*invC; 
  D = yeq0_invC*one_s*phi_to_neg_theta; D1 = D-1; 
  thetaB = theta*B; theta_D1 = theta*D1; 
  
  beta_multiplier = theta*( -yeq0*D*logphi * ( 1 + theta_D1*logphi ) + yneq0*(A+thetaB) )
  gamma_multiplier1 = yeq0 * ( -D*theta*inv_phi1_2 ) * ( phi1*( theta_D1*phi1-2*phi ) + phi_logphi*( 1+phi+theta_D1*( phi_logphi-2*phi1 ) ) )
  gamma_multiplier2 = yneq0 * ( phi*inv_phi1_2 )*( theta*( 2*phi1+(phi+1)*A+phi*thetaB ) - y )
  joint_multiplier = (theta*inv_phi1) * ( yeq0*D*( phi_logphi-phi1 )*( 1+theta_D1*logphi ) + yneq0*( 1-phi*(1+A+thetaB) ) )
  gamma_s_multiplier = yeq0_invC*theta*phi_to_neg_theta*( 1-inv_phi1*phi_logphi )*( 1+one_s*invC*one_phi_to_neg_theta )
  beta_s_multiplier = yeq0_invC*theta*logphi*( D-D1*phi_to_neg_theta )
  hessian.s = -sum( yeq0*(one_phi_to_neg_theta*invC)^2 + yneq0/one_s^2 )
  
  tX = t(X); tU = t(U)
  hessian.beta = tX%*%(X*array(beta_multiplier,c(n,p)))
  hessian.gamma = tU%*%(U*array(gamma_multiplier1+gamma_multiplier2,c(n,q)))
  hessian.joint = tX%*%(U*array(joint_multiplier,c(n,q)))
  hessian.beta_s = tX%*%beta_s_multiplier
  hessian.gamma_s = tU%*%gamma_s_multiplier
  
  hessian = rbind( cbind(hessian.beta, hessian.joint, matrix(hessian.beta_s)),
                   cbind(t(hessian.joint), hessian.gamma, matrix(hessian.gamma_s)),
                   c(hessian.beta_s, hessian.gamma_s, hessian.s) )
  
  colnames(hessian) <- rownames(hessian) <- NULL
  
  return(hessian)
}


################################################################
# Complete likelihood
loglik.ZIMD.complete <- function(beta, gamma, y, X, U, logoffsets, psi, max_val=1e300, tol=1e-12) 
{
  yeq0 = (y==0); ygt0 = !yeq0
  mu  = pmin( exp(X%*%beta  + ygt0*logoffsets), max_val); 
  phi = pmax( exp(U%*%gamma + ygt0*logoffsets), 1+tol); 
  
  # temporal variable that helps caculation
  theta <- mu/(phi-1)

  inv.phi = 1/phi; inv.phi.theta = inv.phi^theta
  
  # y == 0
  l1 <- sum( (1-psi[yeq0]) * log(inv.phi.theta[yeq0]) )
  # y > 0
  l2 <- sum( lgamma(y[ygt0] + theta[ygt0]) - lgamma(theta[ygt0]) - lgamma(y[ygt0] + 1) + log(inv.phi.theta[ygt0]) + y[ygt0]*log(1 - inv.phi[ygt0]))
  return(l1+l2)  
}

# close form gradient
loglik.ZIMD.complete.grad <- function(beta, gamma, y, X, U, logoffsets, psi, max_val=1e300, tol=1e-12) 
{
  X <- as.matrix(X)
  U <- as.matrix(U)
  
  yeq0 = (y==0); ygt0 = !yeq0
  mu  = pmin( exp(X%*%beta  + ygt0*logoffsets), max_val); 
  phi = pmax( exp(U%*%gamma + ygt0*logoffsets), 1+tol); 
  
  # Temporary variables for efficient computation
  theta <- mu/(phi-1)

  logphi=log(phi); inv.phi1 = 1/(phi-1); phi.inv.phi1 = phi*inv.phi1
  inv.phi = 1/phi; inv.phi.theta = inv.phi^theta
  
  B <- digamma(y+theta) - digamma(theta) - logphi
  # multipliers
  Li.beta <- theta* ( ygt0*B - yeq0 * (1-psi) * logphi )
  Li.gamma = -theta* ( yeq0 * ( (1-psi) * (-phi.inv.phi1*logphi + 1) ) + ygt0 * (phi.inv.phi1*B + 1 - y*inv.phi1/theta) )

  # gradient
  grad.beta <- t(X) %*% Li.beta
  grad.gamma <- t(U) %*% Li.gamma
  return(list(grad.beta=grad.beta, grad.gamma=grad.gamma))
}


loglik.ZIMD.complete_helper <- function(par, y, X, U, logoffsets, psi, ...) {
  p <- dim(X)[2]; q <- dim(U)[2]
  return( -loglik.ZIMD.complete(beta=par[1:p], gamma=par[(p+1):(p+q)], y=y, X=X, U=U, logoffsets=logoffsets, psi=psi, ...) )
}

loglik.ZIMD.complete.grad_helper <- function(par, y, X, U, logoffsets, psi, ...) {
  p <- dim(X)[2]; q <- dim(U)[2]
  obj <- loglik.ZIMD.complete.grad(beta=par[1:p], gamma=par[(p+1):(p+q)], y=y, X=X, U=U, logoffsets=logoffsets, psi=psi, ...)
  return( -c(obj$grad.beta, obj$grad.gamma) )
}
