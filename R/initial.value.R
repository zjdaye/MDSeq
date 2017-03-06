# Computes initial values based on robust moment estimates
# If contrast == NULL, then no contrasts
get.initial.pars <- function(y, contrast=list(mean=NULL,dispersion=NULL), X, U, logoffsets, btol, qtol, trim)
{
  require(quadprog)
  n=length(y); if (!is.null(X)) {p = ncol(X)} else{p = 0}; if (!is.null(U)) {q = ncol(U)} else{q = 0}
  
  #########################################
  if (!is.null(contrast$mean)) {
    unique_contrast=unique(contrast$mean); unique_len=dim(unique_contrast)[1]; 
    mu = rep(NA, n)
    for (j in 1:unique_len) {
      contrastTmp = unique_contrast[j,];                 
      idx = which( apply(contrast$mean, 1, function(x) all(x == contrastTmp)) ); 
      mu[idx] = mean(y[idx], trim=trim)
    }
    mu <- pmax(mu, btol)
    beta0 = lm((log(mu) - logoffsets) ~ contrast$mean - 1)$coeff
    names(beta0) = colnames(contrast$mean)
  } else { mu.total = max(mean(y/exp(logoffsets), trim=trim),btol); beta0 = log(mu.total); names(beta0)="(Intercept)" } # in the null, there is no contrast
  
  if (!is.null(contrast$dispersion)) {
    unique_contrast1=unique(contrast$dispersion); unique_len=dim(unique_contrast1)[1]; 
    mu = rep(NA, n); D = rep(NA, n)
    for (j in 1:unique_len) {
      contrastTmp = unique_contrast1[j,]; 
      idx = which( apply(contrast$dispersion, 1, function(x) all(x == contrastTmp)) ); 
      yTmp = y[idx]; D[idx] = mad(yTmp)^2/mean(yTmp, trim=trim)
    }
    D <- pmax(D, 1+btol)
    gamma0 = lm((log(D)-logoffsets) ~ contrast$dispersion - 1)$coeff
    #if (anyNA(gamma0)) {gamma0[-1] <- NA}
    names(gamma0) = colnames(contrast$dispersion)
  } else { D.total = max(mad(y/exp(logoffsets))^2/mean(y/exp(logoffsets), trim=trim),1+btol); gamma0 = log(D.total); names(gamma0)="(Intercept)" } # in the null, there is no contrast
  
  # default values when values missing  
  beta0[is.na(beta0) | is.infinite(beta0) | is.nan(beta0)] = btol
  gamma0[is.na(gamma0) | is.infinite(gamma0) | is.nan(gamma0)] = btol
  
  # Pad coefficients on covariates with 0's, e.g. assuming covariates are centered
  beta0 = c(beta0, rep(0,times=p-length(beta0)))
  gamma0 = c(gamma0, rep(0,times=q-length(gamma0)))
  
  #########################################
  ## Reset initial values to feasible region
  reset.obj = reset.initial.pars(list(beta=beta0, gamma=gamma0), y, X, U, logoffsets, btol, qtol)
  
  return(reset.obj)
}


##################################################################
# Computes initial values based on robust moment estimates
# If contrast == NULL, then no contrasts
reset.initial.pars <- function(initial.pars, y, X, U, logoffsets, btol, qtol) 
{
  n=length(y); p = ncol(X); q = ncol(U)
  #########################################
  ## Reset initial values to feasible region
  x0 = c(initial.pars$beta, initial.pars$gamma); min_mu = btol; min_phi = 1 + btol; 
  max_y = max(y/exp(logoffsets)); max_mu = max_y + 1; max_phi = max_y^2/(4*min_mu)+1
  
  # Initial linear inequality variables
  # Set A = (-X  0
  #           X  0
  #           0  U
  #           0 -U)
  #         X%*%beta > min_mu
  #         U%*%gamma > min_phi
  #        -U%*%gamma > -max_phi
  Zero_p = array(0, c(n,p)); Zero_q = array(0, c(n,q))
  Amat = as.matrix(rbind( cbind(X, Zero_q), 
                          cbind(-X, Zero_q),
                          cbind(Zero_p, U),
                          cbind(Zero_p, -U)))
  bvec = c(rep(-log(max_mu), n)- logoffsets, rep(log(min_mu),n)- logoffsets,
           rep(log(min_phi),n)- logoffsets,rep(-log(max_phi),n)- logoffsets)
  # Set optimization initial value x to be parameters 
  # closest to give initial value x0 under constraints by 
  #   min || x - x0 || subject to Ax >= 0; equivalently, 
  #   => min x^Tx - 2x0^T x subject to Ax >= 0
  bvecRobust = c(rep(-log(max_mu) - qtol, n)- logoffsets, rep(log(min_mu)+ qtol ,n)- logoffsets,
                 rep(log(min_phi) + qtol,n)- logoffsets,  rep(-log(max_phi) - qtol,n)- logoffsets) # Add qtol to constraint, since QP can violate constraint
  if ( any( (Amat %*% as.vector(x0)) <= bvec) ) {# If any inequality constraint is violate
    obj.qp = solve.QP(Dmat=diag(p+q), dvec=x0, Amat=t(Amat), bvec=bvecRobust)
    theta = obj.qp$solution
  } else { theta=x0 }
  
  return(list(theta = theta,
              Amat = Amat,
              bvec = bvec))
}
