#' Checking outliers
#' 
#' Identifies outliers influential upon the effects of a covariate 
#' or a set of covariates for an individual gene.  
#' 
#' @param y vector containing the counts of a gene.
#' @param X vector or matrix containing covariates in the mean GLM.
#' @param U vector or matrix containing covariates in the dispersion GLM, 
#'          which can be different from X.
#' @param contrast list containing contrast matices for mean and dispersion.
#' @param offsets numeric vector of same size as y giving offsets for use in log-linear models.
#' @param check which covariate(covariates) to check. if check="contrast", outliers 
#'              influential on the effects of contrasts will be checked; if check="other", outliers 
#'              influential on the effects of all other covariates (except contrasts) will be checked; 
#'              if check="all", outliers influential on the effects of all covariates will be checked. 
#'              In addition, outliers influential towards specific covariate(/covariates) can be checked. 
#'              User needs to provide the name(s) of targeted covariate(s).  
#' @param cut.off cutoff point defining the significance level of outlier selections. 
#' @param ZI logical, if TRUE (default), the effect of zero-inflation will be considered. 
#' @param conserv logical, if TRUE (default), outlier detection will skip any gene with estimation 
#'                         error for either coefficients or standard errors of coefficients and label it as an error; 
#'                         if FALSE, outlier detection will only skip any gene with estimation error 
#'                         on coefficients. 
#' @param verbose logical, if TRUE, progress information will be reported.
#' @param output logical, if TRUE, one-step estimator of difference in deviance will be reported.
#' @param ... other arguments could be passed in.
#' 
#' @return list containing of the result of outlier checking.
#' \item{outliers}{indicator of outlier}
#' \item{I2}{one-step estimator of difference in deviance}
#' 
#' @export
#' 

check.outlier <- function(y, X = NULL, U = NULL, contrast = list(mean=NULL,dispersion=NULL), offsets = NULL,
                          check = 'contrast', cut.off=outlier.cutoff(0.05), 
                          ZI = TRUE, conserv = TRUE, verbose = TRUE, output = FALSE, ...)
{
  if (!any(check %in% c(colnames(X), colnames(U), 'contrast', 'all', 'other'))) 
      stop('Could not find the variable! Please choose \'contrast\' or other valid covariate.')
  
  check.obj = check_inputs(y, contrast, X, U, initial.pars=list(beta=NULL,gamma=NULL), offsets, ...)
  
  y = check.obj$y; contrast=check.obj$contrast; X=check.obj$X; U=check.obj$U; offsets=check.obj$offsets
  
  n <- length(y)
  
  beta.H0 <- gamma.H0 <- beta.H1 <- gamma.H1 <- X.H0 <- U.H0 <- NA
  
  fit.H1 <- try(glm.ZIMD(y, X=X, U=U, contrast = contrast, offsets=offsets, ZI = ZI, ...), silent = T)
  
  if ('all' %in% check) {
    
    if(verbose) message('Checking outliers related to contrast and covariates...')
    contrast.H0 <- list(mean=NULL,dispersion=NULL)
    X.H0 <- U.H0 <- NULL
    
  } else {
    
    if ( ('contrast' %in% check) && (!(is.null(contrast$mean) || is.null(contrast$dispersion))) ) {
      
      if(verbose) message('Checking outliers related to contrast...')
      
      contrast.H0 <- list(mean=NULL,dispersion=NULL)
      X.H0 <- X
      U.H0 <- U
      
    } else {
      
      contrast.H0 <- contrast
      X.H0 <- X
      U.H0 <- U
      
      if ('other' %in% check) {
        
        if(verbose) message('Checking outliers related to all covariates...')
        X.H0 <- U.H0 <- NULL
        
      } else {
        
        if ( any(check %in% colnames(X)) || any(check %in% colnames(U))) {   
          
          if (any(check %in% colnames(X))) {
            if(verbose) message('Checking outliers related to covariate (covariate for mean)...')
            X.H0 <- as.matrix(X[,colnames(X)[!(colnames(X) %in% check)]])
            colnames(X.H0) <- colnames(X)[!(colnames(X) %in% check)]
          }
          
          if (any(check %in% colnames(U))) {
            if(verbose) message('Checking outliers related to covariate (covariate for dispersion)...')
            U.H0 <- as.matrix(U[,colnames(U)[!(colnames(U) %in% check)]])
            colnames(U.H0) <- colnames(U)[!(colnames(U) %in% check)]
          }
          
        }
        
      }
      
    }
    
  }
  
  fit.H0 <- try(glm.ZIMD(y, X=X.H0, U=U.H0, contrast = contrast.H0, offsets=offsets, ZI = ZI,...), silent = T)
  
  if (!inherits(fit.H0,'try-error') && !inherits(fit.H1,'try-error')) {
    
    beta.H0 <- fit.H0$coefficients$beta
    gamma.H0 <- fit.H0$coefficients$gamma
    beta.H1 <- fit.H1$coefficients$beta
    gamma.H1 <- fit.H1$coefficients$gamma
    
    se.beta.H0 <- fit.H0$se$se.beta
    se.gamma.H0 <- fit.H0$se$se.gamma
    se.beta.H1 <- fit.H1$se$se.beta
    se.gamma.H1 <- fit.H1$se$se.gamma
    
  } 
  
  if (anyNA(c(beta.H0,beta.H1,gamma.H0,gamma.H1))) {
    
    stop('Cannot check outliers!')
    
  } else {
    
    if (conserv && anyNA(c(se.beta.H0,se.beta.H1,se.gamma.H0,se.gamma.H1))) {
      
      stop('Cannot check outliers!')

    } else {
      
      rg2.H1 <- rg2est(beta.H1, gamma.H1, y, X, U, contrast, offsets, ZI)
      
      rg2.H0 <- rg2est(beta.H0, gamma.H0, y, X.H0, U.H0, contrast.H0, offsets, ZI)
      
      I2 <- rg2.H1 - rg2.H0
      
      outliers <- which( abs(I2) > cut.off )
      
      if (output) {
        return(list(outliers=outliers, I2=I2))
      } else {
        return(as.integer(outliers))
      }
      
    }
  }
}





############################################################################
rg2est <- function(beta, gamma, y, X, U, contrast, offsets, ZI)
{
  n <- length(y)
  
  if (!is.null(contrast$mean)) {
    X <- as.matrix(cbind(contrast$mean, X))
  } else {
    X <- as.matrix(cbind(rep(1, n), X))
    colnames(X)[1] <- '(Intercept)'
  }
  if (!is.null(contrast$dispersion)) {
    U <- as.matrix(cbind(contrast$dispersion, U))
  } else {
    U <- as.matrix(cbind(rep(1, n), U))
    colnames(U)[1] <- '(Intercept)'
  }
  
  rg2.out <- rep(NA, n)
  y1 <- y
  
  px <- dim(X)[2]
  pu <- dim(U)[2]
  
  if (ZI) {
    X <- X[y>0,]
    U <- U[y>0,]
    offsets <- offsets[y>0]
    y <- y[y>0]
  }
  
  n <- length(y)
  logoffsets <- log(offsets)
  
  # mu and D depends on v, d, beta, gamma, contr and X
  # nx1 vectors
  if (px==1) {
    mu <- exp(X*beta + logoffsets) 
  } else {
    mu <- exp(X%*%beta + logoffsets)            # mean
  }
  if (pu==1) {
    D <- exp(U*gamma + logoffsets)
  } else {
    D  <- exp(U%*%gamma + logoffsets)           # dispersion
  }
  
  V  <- mu*D                     # variance
  
  # weight - iterated weighted least squares
  # nx1 vector
  w <- mu/D
  
  # Weight matrix
  #W <- diag(w[1:n])
  W <- diag(w)
  
  # Hat matrix
  sqrtDiagW <- sqrt(w)
  sqrtWX <- X*array(sqrtDiagW, c(n,px))
  H <- sqrtWX %*% solve( t(sqrtWX) %*% sqrtWX ) %*% t(sqrtWX)
  h <- diag(H)
  oneh = 1-h
  
  # deviance
  dev <- res <- r.y <- r  <- p <- loglik.y <- loglik.mu <- rep(NA,n)
  
  D1 = D-1
  r.y <- y/D1                          
  r   <- pmax( mu/D1, 1e-5 )                         
  p   <- pmax( D1/D, 1e-5 )
  
  loglik.y  <- logpnb(r.y, p, y)
  loglik.y[is.na(loglik.y)] <- 0      #when r=0, loglik=NaN, need to change it to 0
  loglik.mu <- logpnb(r, p, y)
  loglik.mu[is.na(loglik.mu)] <- 0
  
  # dev takes the same sign as y - mu_hat
  res = y-mu
  dev <- sign( res ) * sqrt( 2*D*( abs( loglik.y - loglik.mu )))
  
  # Pearson residuals, standardised to have unit asymptotic variance
  rp  <- res/sqrt(V*oneh)
  
  # standardised deviance residuals
  rd  <- dev/sqrt(D*oneh)
  
  # likelihood ratio statistic reduction when ith case is delete
  rg2 <- oneh * rd^2 + h * rp^2
  
  if (ZI) { 
    rg2.out[y1>0] <- rg2
  } else {
      rg2.out <- rg2
    }
  
  return(rg2.out)
}


# log probability mass function - NB
logpnb <- function(r, p, k){
  suppressWarnings( loglik <- lgamma(r+k)-lgamma(r)-lgamma(k+1)+k*log(p)+r*log(1-p) )
  return(loglik)
}



