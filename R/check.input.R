# check inputs
check_inputs <- function(y, contrast, X, U, initial.pars, offsets, control = list())
{
  if (all(y==0)) stop('All input values are zero!')
  # Set control parameters
  # defaults, will be reset to new values
  con = list(btol=1e-5, qtol=1e-4, trim=0.05, # passed to get.initial.pars
             maxit=100, reltol=1e-8, # passed to optim
             max_val=1e300, tol=1e-12 # passed to log-likelihood, gradient, and hessian
  )
  idx = match(names(control), names(con))
  if (length(idx)>0) { for (j in 1:length(idx)) { idxTmp = idx[j]; if (is.na(idxTmp)) next; con[idxTmp] = control[j] } }
  
  if (!is.null(contrast$mean)) cM <- attributes(contrast$mean)$contrast$group
  if (!is.null(contrast$dispersion)) cD <- attributes(contrast$dispersion)$contrast$group
  
  y = as.integer(floor(y+.5)); n = length(y);
  if (!is.null(X)) { X <- as.matrix(X) }; if(!is.null(U)) { U <- as.matrix(U) };
  # Check if matching sample sizes
  if (!is.null(X)) { if (n!=nrow(X)) { stop("Sample size does not match between response and covariates!") } }
  if (!is.null(U)) { if (n!=nrow(U)) { stop("Sample size does not match between response and covariates!") } }
  if (!is.null(contrast$mean)) { if (n!=nrow(as.matrix(contrast$mean))) { stop("Sample size does not match between response and contrasts!") }  }
  if (!is.null(contrast$dispersion)) { if (n!=nrow(as.matrix(contrast$dispersion))) { stop("Sample size does not match between response and contrasts!") }  }
  if (!is.null(offsets)) { if (n!=length(offsets)) {stop("Sample size does not match between response and offsets!")}}
  if (!is.null(offsets)) { if (any(offsets <= 0)) {stop("Offsets must be greater than 0!")}}
  
  # set default offsets
  if (is.null(offsets)) { offsets <- rep(1, length=n) }
  
  # Check for missing values and remove them
  idx.y = which(is.na(y) | y<0 | is.infinite(y) | is.nan(y)); 
  idx.X = idx.U = NULL
  if (!is.null(X)) {idx.X = which(apply( is.na(X) | is.infinite(X) | is.nan(X), 1, any))}; 
  if (!is.null(U)) {idx.U = which(apply( is.na(U) | is.infinite(U) | is.nan(U), 1, any))}
  idx = c(idx.y, idx.X, idx.U)
  if (length(idx)>0) { 
    warning(c("In checking inputs, ",
              length(idx), " observations with missing values have been removed."), call.=F); 
    y=y[-idx]; X=X[-idx,]; U=U[-idx,]; n=n-length(idx); 
    if (!is.null(contrast$mean)) { contrast$mean=as.matrix(contrast$mean)[-idx,]}
    if (!is.null(contrast$dispersion)) { contrast$dispersion=as.matrix(contrast$dispersion)[-idx,]}
    offsets = offsets[-idx]
  }
  
  if (!is.null(initial.pars$beta)) { 
    if (length(initial.pars$beta)!=(ncol(contrast$mean)+ncol(X))) { stop("Dimension of initial.pars$beta incorrect!") }           # both beta and gamma need to be given? 
    if (is.null(initial.pars$gamma)) { stop("Missing initial.pars$gamma") }                                                               # if only contrast$mean was specified
  }
  if (!is.null(initial.pars$gamma)) { 
    if (length(initial.pars$gamma)!=(ncol(contrast$dispersion)+ncol(U))) { stop("Dimension of initial.pars$beta incorrect!") } 
    if (is.null(initial.pars$beta)) { stop("Missing initial.pars$beta") }                                                                 # if only contrast$dispersion was specified
  }
  
  if (!is.null(contrast$mean)) attributes(contrast$mean)$contrast$group <- cM
  if (!is.null(contrast$dispersion)) attributes(contrast$dispersion)$contrast$group <- cD
  
  return(list(y=y, contrast=contrast, X=X, U=U, offsets=offsets, control=con))
}
