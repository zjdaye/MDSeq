#' Zero-Inflated Mean-Dispersion Generalized Linear Models 
#'
#' Fit zero-inflated mean-dispersion generalized linear models (GLMs) on 
#' both the mean and dispersion for an individual gene. Maximum 
#' likelihood estimations will be calculated via EM algorithm 
#' incorporated with linearly constrained optimization 
#' using an adaptive barrier algorithm.
#' 
#' @param y vector contains the counts of one gene
#' @param X vector or matrix contains covariates in the mean GLM
#' @param U vector or matrix contains covariates in the dispersion GLM, 
#'          which can be different from X.
#' @param contrast list containing design matices for mean and dispersion.
#' @param initial.pars list giving initial values for parameters of mean and/or dispersion.
#' @param offsets numeric vector of same size as y giving offsets for use in log-linear models.
#' @param constrOptim.maxIter maximum iterations of the barrier algorithm in constrained optimization.
#' @param constrOptim.tol non-negative number; the relative convergence 
#'                        tolerance of the barrier algorithm in constrained optimization.
#' @param control other control parameters for constrained optimization.
#' @param ZI logical, if True (defualt), the effect of zero-inflation will be considered. 
#'            If ZI=False, the effect of zero-inflation will be ignored in computation. 
#' @param EM.maxIter maximum iteration for EM algorithm.
#' @param EM.tol convergence threshold for EM algorithm.
#' 
#' 
#' @return an object of class 'glmZIMD' containing the following information:
#' \item{call}{formula}
#' \item{contrasts}{contrast}
#' \item{coefficients}{a list containing estimated coefficients from the GLM fits, in the natural log scale}
#' \item{se}{a list containing standard errors of coefficients}
#' \item{s}{zero-inflation parameter and standard error}
#' \item{loglik}{log-likelihood}
#' \item{hessian}{hessian matrix}
#' \item{convergence}{an integer code returned by the optimizer. 0 indicates successful completion}
#' \item{outer.iterations}{number of iterations of the barrier algorithm}
#' \item{EM.convergence}{EM indicator}
#' \item{EM.iter}{numbers of EM iterations}
#' 
#' @examples
#' # simulate 1000 genes, each has 200 observations.
#' n <- 200; nsim <- 1000; seed <- 1234
#' 
#' # set covariate coefficient for mean and dispersion
#' beta <- c(5, log(2), 0.5, 0.5); gamma <- c(4, 0, 0.5, 0.5)
#' 
#' # set proportion of excessive-zero
#' s <- 0.2
#' 
#' # set design matrix for mean and dispersion 
#' factors <- c('Control', 'Case')
#' group <- gl(n=2, k=n/2, labels=factors)
#' design <- get.model.matrix(group, contrast.type="contr.treatment")
#' 
#' # simulate data nsim (genes) X n (observations)
#' ysim <- sim.ZIMD(n=n, nsim=nsim, beta=beta, gamma=gamma, contrast = design, s=s, seed=seed)
#' 
#' # one example and covariates for mean and dispersion
#' test <- ysim$y[1,]
#' X <- ysim$X[,3:4]; colnames(X) = c("X3", "X4")
#' U <- ysim$U[,3:4]; colnames(U) = c("U3", "U4")
#' 
#' # estimation for a single gene
#' # consider zero-inflation
#' fit1 <- glm.ZIMD(test, X=X, U=U, ZI=T)
#' # ignore zero-inflation
#' fit2 <- glm.ZIMD(test, X=X, U=U, ZI=F)
#' 
#' fit1
#' fit2
#' 
#' @export
#' 



glm.ZIMD <- function(y, X = NULL, U = NULL, contrast=list(mean=NULL,dispersion=NULL), 
                     initial.pars=list(beta=NULL,gamma=NULL), offsets=NULL, constrOptim.maxIter=100, 
                     constrOptim.tol=1e-05, ZI=TRUE, EM.maxIter=200, EM.tol=1e-04, control=list())
{
  require(gtools)
  # check inputs
  check.obj = check_inputs(y, contrast, X, U, initial.pars, offsets, control)
  y = check.obj$y; contrast=check.obj$contrast; X=check.obj$X; U=check.obj$U; control=check.obj$control; offsets=check.obj$offsets
  logoffsets = log(offsets)
  
  n <- length(y)
  dx <- du <- 0 
  
  # Standardize covariates for numerical stability
  scale.X <- scale.U <- NULL
  if (!is.null(X)){
    X <- as.matrix(X); dx <- ncol(X)
    scale.X = scale(X); X.colmean=attr(scale.X, "scaled:center"); X.colsd = attr(scale.X, "scaled:scale") 
    if (!is.null(initial.pars$beta)) {
      idx = (length(initial.pars$beta)-length(X.colsd)+1):length(initial.pars$beta)
      initial.pars$beta[1] = initial.pars$beta[1] + sum( X.colmean*initial.pars$beta[idx] )
      initial.pars$beta[idx] = initial.pars$beta[idx]*X.colsd
    }
  }
  if (!is.null(contrast$mean)) {
    X <- as.matrix(cbind(contrast$mean, scale.X))
  } else {
    X <- as.matrix(cbind(rep(1, n), scale.X))
    colnames(X)[1] <- '(Intercept)'
  }
  
  if(!is.null(U)){
    U <- as.matrix(U); du <- ncol(U)
    scale.U = scale(U); U.colmean=attr(scale.U, "scaled:center"); U.colsd = attr(scale.U, "scaled:scale") 
    if(!is.null(initial.pars$gamma)) {
      idx = (length(initial.pars$gamma)-length(U.colsd)+1):length(initial.pars$gamma)
      initial.pars$gamma[1] = initial.pars$gamma[1] + sum( U.colmean*initial.pars$gamma[idx] )
      initial.pars$gamma[idx] = initial.pars$gamma[idx]*U.colsd
    }
  }
  if (!is.null(contrast$dispersion)) {
    U <- as.matrix(cbind(contrast$dispersion, scale.U))
  } else {
    U <- as.matrix(cbind(rep(1, n), scale.U))
    colnames(U)[1] <- '(Intercept)'
  }
  
  p=ncol(X); q=ncol(U) # Dimensions
  
  beta <- se.beta <- rep(NA,p); se.gamma <- gamma <- rep(NA,q); par <- c(beta, gamma)
  hes <- se <- loglik <- convergence <- outer.iter <- iter <- EM.convergence <- se.s <- NA
  ygt0 = y>0 # Boolean for indices of nonzeros
  yeq0 = !ygt0 # Boolean for indices of 0's
  
  # Does not use zero-inflation
  if ( all(ygt0) || (ZI == FALSE) ) {
    s = 0  # no zero-inflation
    # initial values for constrained optimization using moment estimates
    if(is.null(initial.pars$beta) & is.null(initial.pars$gamma)) { 
      ini = get.initial.pars(y, contrast, X, U, logoffsets, control$btol, control$qtol, control$trim) 
    } else { ini = reset.initial.pars(initial.pars, y, X, U, logoffsets, control$btol, control$qtol) }
    
    # Constrained optimization
    obj <- try( constrOptim(theta=ini$theta, f=loglik.MD_helper, grad=loglik.MD.grad_helper,
                            ui=ini$Amat, ci=ini$bvec, method="BFGS", y=y, X=X, U=U, logoffsets=logoffsets,
                            outer.iterations = constrOptim.maxIter, outer.eps = constrOptim.tol,
                            control=list(reltol=control$reltol, maxit=control$maxit),
                            max_val=control$max_val, tol=control$tol), silent = TRUE)
    if(inherits(obj,'try-error')) {
      obj <- try( constrOptim(theta=ini$theta, f=loglik.MD_helper, grad=loglik.MD.grad_helper,
                              ui=ini$Amat, ci=ini$bvec, method="Nelder-Mead", y=y, X=X, U=U, logoffsets=logoffsets,
                              outer.iterations = constrOptim.maxIter, outer.eps = constrOptim.tol,
                              control=list(reltol=control$reltol, maxit=control$maxit),
                              max_val=control$max_val, tol=control$tol), silent = TRUE)      
    }
    
  } else {
    # has zero-inflation
    s <- sum(yeq0)/n  # initial guess of the proportion of excess zero
    
    # initial values for constrained optimization 
    if(is.null(initial.pars$beta) & is.null(initial.pars$gamma)) { 
      contrast_r0 = list(mean=contrast$mean[ygt0,], dispersion=contrast$dispersion[ygt0,])
      ini = get.initial.pars(y[ygt0], contrast_r0, as.matrix(X[ygt0,]), as.matrix(U[ygt0,]), logoffsets[ygt0], control$btol, control$qtol, control$trim) 
    } else { ini = reset.initial.pars(initial.pars, round(y[ygt0]/offsets[ygt0]), as.matrix(X[ygt0,]), as.matirx(U[ygt0,]), logoffsets[ygt0], control$btol, control$qtol) }
    
    beta_0 <- ini$theta[1:p]
    gamma_0 <- ini$theta[(p+1):(p+q)]
    
    status <- 0
    loglik <- -Inf
    
    for(iter in 1:EM.maxIter){
      loglik_pre <- loglik
      
      # E-step
      mu  = pmin( exp(X%*%beta_0  + ygt0*logoffsets), control$max_val); 
      phi = pmax( exp(U%*%gamma_0 + ygt0*logoffsets), 1+control$tol); 
      theta <- mu/(phi-1)
      
      # probability at y==0
      p_0 <- phi^(-theta)
      
      # probability in zero-state
      psi <- rep(0, n)
      psi[yeq0] <- 1/(1 + (1/s - 1) * p_0[yeq0])
      
      # M-step: updating s
      s <- sum(psi)/n
      
      # M-step: updating beta and gamma
      obj <- try( constrOptim(theta=ini$theta, f=loglik.ZIMD.complete_helper, grad=loglik.ZIMD.complete.grad_helper,
                              ui=ini$Amat, ci=ini$bvec, method="BFGS", y=y, X=X, U=U, logoffsets=logoffsets, psi=psi,
                              outer.iterations = constrOptim.maxIter, outer.eps = constrOptim.tol,
                              control=list(reltol=control$reltol, maxit=control$maxit),
                              max_val=control$max_val, tol=control$tol), silent = TRUE)
      
      if(inherits(obj,'try-error')) {
        obj <- try( constrOptim(theta=ini$theta, f=loglik.ZIMD.complete_helper, grad=loglik.ZIMD.complete.grad_helper,
                                ui=ini$Amat, ci=ini$bvec, method="Nelder-Mead", y=y, X=X, U=U, psi=psi,
                                outer.iterations = constrOptim.maxIter, outer.eps = constrOptim.tol,
                                control=list(reltol=control$reltol, maxit=control$maxit),
                                max_val=control$max_val, tol=control$tol), silent = TRUE)
      } 
      
      # if no result from optimization, program will go to next iteration
      if(!inherits(obj,'try-error')){
        pars <- obj$par
        loglik <- loglik.ZIMD_helper(c(pars,s),y=y,X=X,U=U,logoffsets=logoffsets)
        diff <- (abs(loglik-loglik_pre) < EM.tol)  # convergence criteria
      }else{
        diff <- FALSE 
      }
      
      # if convergence criteria met, iteration stops.
      if(diff) {status = 1; break;}
    }
    
    EM.convergence <- 1-status
  }
  
  # results
  if (!inherits(obj,'try-error')) {
    par = obj$par; beta = par[1:p]; gamma = par[(p+1):(p+q)]
    convergence <- obj$convergence; outer.iter <- obj$outer.iterations
    
    if (s==0) { loglik <- -1 * obj$value; hes <- loglik.MD.hessian(beta=beta, gamma=gamma, y=y, X=X, U=U, logoffsets=logoffsets)
    } else { hes <- loglik.ZIMD.hessian(beta=beta, gamma=gamma, s=s, y=y, X=X, U=U, logoffsets=logoffsets) }
    
    se <- try(suppressWarnings(sqrt(-diag(solve(hes)))), silent = TRUE)
    if(inherits(se,'try-error')) se <- rep(NA,(p+q))
  }
  
  se.beta=zvalue.beta=pvalue.beta=rep(NA, p); se.gamma=zvalue.gamma=pvalue.gamma=rep(NA, q); 
  #if(!anyNA(se)) {se.beta <- se[1:p]; se.gamma <- se[(p+1):(p+q)]}
  # no NA's are allowed in either se.beta and se.gamma
  #if(!anyNA(se[1:p])) se.beta <- se[1:p]
  #if(!anyNA(se[(p+1):(p+q)])) se.gamma <- se[(p+1):(p+q)]
  se.beta = se[1:p]; se.gamma = se[(p+1):(p+q)]  # Can keep NA on individual coefficients
  
  p0 = p - dx; q0 = q - du;
  names(beta) <- colnames(X)
  if (!is.null(contrast$mean)) { names(beta)[2:p0] <- paste('Contrast$mean', names(beta)[2:p0], sep='.')  }
  names(gamma) <- colnames(U)
  if (!is.null(contrast$dispersion)) { names(gamma)[2:q0] <- paste('Contrast$dispersion', names(gamma)[2:q0], sep='.')  }
  names(se.beta) <- names(zvalue.beta) <- names(pvalue.beta) <- names(beta)
  names(se.gamma) <- names(zvalue.gamma) <- names(pvalue.gamma) <- names(gamma)
  #if (!is.null(contrast$mean)) { names(beta)=c('Intercept', 'Contrast$mean', colnames(X)[(p0+1):p])
  #} else{ names(beta)=c('Intercept', colnames(X)[(p0+1):p]) }
  #if (!is.null(contrast$dispersion)) { names(gamma) = c('Intercept', 'Contrast$dispersion', colnames(U)[(q0+1):q]) 
  #} else{ names(gamma)=c('Intercept', colnames(U)[(q0+1):q])}
  #names(se.beta)=names(zvalue.beta)=names(pvalue.beta)=names(beta)
  #names(se.gamma)=names(zvalue.gamma)=names(pvalue.gamma)=names(gamma)
  
  # Re-scale standardized estimates to original scale
  vc <- try(solve(-hes),silent = T) # variance-covariance matrix
  if(!inherits(vc, 'try-error')) {
    if((p > p0) && !anyNA(se.beta) ){ 
      beta[(p0+1):p] = beta[(p0+1):p]/X.colsd; se.beta[(p0+1):p] = se.beta[(p0+1):p]/X.colsd # estimates of covariates
      # estimates of intercept
      beta[1] <- beta[1] - beta[(p0+1):p] %*% X.colmean
      vcX <- vc[1:p,1:p]; if (!is.null(contrast$mean)) { vcX <- vcX[c(1,((p0+1):p)),c(1,((p0+1):p))] }
      aX <- c(1,-X.colmean/X.colsd)
      se.beta[1] <- sqrt(t(aX) %*% vcX %*% aX)
    }
    if( (q > q0) && !anyNA(se.gamma) ){ # gamma
      gamma[(q0+1):q] = gamma[(q0+1):q]/U.colsd; se.gamma[(q0+1):q] = se.gamma[(q0+1):q]/U.colsd # estimates of covariates
      # estimates of intercept
      gamma[1] <- gamma[1] - gamma[(q0+1):q] %*% U.colmean
      vcU <- vc[(p+1):(p+q),(p+1):(p+q)]; if (!is.null(contrast$dispersion)) { vcU <- vcU[c(1,((q0+1):q)),c(1,((q0+1):q))]}
      aU <- c(1,-U.colmean/U.colsd)
      se.gamma[1] <- sqrt(t(aU) %*% vcU %*% aU)
    }
    zvalue.beta = beta/se.beta; zvalue.gamma = gamma/se.gamma
    pvalue.beta = 2*pnorm(-abs(zvalue.beta)); pvalue.gamma = 2*pnorm(-abs(zvalue.gamma));
  }
  
  allcM <- allcD <- NULL
  if ((p0 > 1) && !anyNA(se.beta)) {
    VCM <- vc[1:p0, 1:p0]
    betaM <- beta[1:p0]
    cmatM <- attributes(contrast$mean)$contrast$group
    levM <- 1:nrow(cmatM); names(levM) <- rownames(cmatM)
    allpM <- permutations(n=nrow(cmatM),r=2,v=names(levM), repeats.allowed = T)
    allcM <- as.vector(t(apply(allpM, 1, VCf, levM, cmatM, betaM, VCM)))
    nameM <- apply(allpM, 1, function(x) paste(x, collapse = 'vs'))
    names(allcM) <- c(nameM, paste('se', nameM, sep = '.'))
  }
  if ((q0 > 1) && !anyNA(se.gamma)){
    VCD <- vc[(p+1):(p+q0), (p+1):(p+q0)]
    gammaD <- gamma[1:q0]
    cmatD <- attributes(contrast$dispersion)$contrast$group
    levD <- 1:nrow(cmatD); names(levD) <- rownames(cmatD)
    allpD <- permutations(n=nrow(cmatD),r=2,v=names(levD), repeats.allowed = T)
    allcD <- as.vector(t(apply(allpD, 1, VCf, levD, cmatD, gammaD, VCD)))
    nameD <- apply(allpD, 1, function(x) paste(x, collapse = 'vs'))
    names(allcD) <- c(nameD, paste('se', nameD, sep = '.'))
  }
  
  coefficients.contrast <- list(beta=allcM, gamma=allcD)
  
  coefficients <- list(beta=beta, gamma=gamma)
  stderr <- list(se.beta=se.beta, se.gamma=se.gamma)
  zvalue = list(zvalue.beta=zvalue.beta, zvalue.gamma=zvalue.gamma)
  pvalue = list(pvalue.beta=pvalue.beta, pvalue.gamma=pvalue.gamma)
  
  slist <- list(s=s, se.s=se[p+q+1])
  
  rslt <- list(coefficients=coefficients,
               se=stderr,
               zvalue=zvalue,
               pvalue=pvalue, 
               loglik=loglik,
               hessian=hes,
               convergence=convergence,
               outer.iterations=outer.iter,
               coefficients.contrast=coefficients.contrast)
  
  rslt$EM.convergence <- EM.convergence
  rslt$EM.iter <- iter
  rslt$s <- slist
  
  rslt$call <- match.call()
  rslt$contrasts <- contrast
  
  class(rslt) <- "glmZIMD"
  
  return(rslt)
}





VCf <- function(x, lev, cmat, beta, variance) {
  beta <- beta[-1]
  variance <- variance[-1,-1]
  x <- lev[x]
  theta <- beta %*% (cmat[x[1],] - cmat[x[2],]) 
  var = t(cmat[x[1],] - cmat[x[2],]) %*% variance %*% (cmat[x[1],] - cmat[x[2],])
  return(c(theta,sqrt(var)))
}