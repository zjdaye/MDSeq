#' Goodness-of-fit test for zero-inflation
#' 
#' This function will perform likelihood ratio test for the presence of technical excess zeros 
#' at each gene.
#' 
#' @param counts matrix of count data
#' @param X vector or matrix contains covariates in the mean GLM.
#' @param U vector or matrix contains covariates in the dispersion GLM.
#' @param contrast list containing contrast for mean and dispersion.
#' @param alpha threshold of significance level.
#' @param initial.pars list giving initial values for parameters of 
#'                    mean and/or dispersion for each gene.
#' @param offsets numeric vector of same size as rows of counts giving offsets for use in log-linear models.
#' @param verbose logical, if True, progress will be reported.
#' @param mc.cores integer equals to or greater than 1 (default), indicating how many cores 
#'                  will be used in computation. If mc.cores > 0, then parallel processes 
#'                  will be used, based on the 'parallel' package.
#' @param ... other arguments could be passed in.
#' 
#' @return a data frame containing testing result of zero-inflation.
#' \item{ZI}{zero-inflation indicator.
#'           0, not zero-inflated.
#'           1, zero-inflated.
#'           2, error when testing}
#' \item{pval}{p-value of likelihood ratio test}
#' 
#' @export
#' 



gof.ZI <- function(counts, X = NULL, U = NULL, contrast = list(mean=NULL,dispersion=NULL), alpha=0.05,
                   initial.pars = list(beta = NULL, gamma = NULL), offsets = NULL, verbose = TRUE, mc.cores = 1, ...)
{
  startTime <- proc.time()
  
  counts <- as.data.frame(counts)
  
  if (!is.null(X)) {
    if (ncol(counts) != nrow(X)) stop('Sample size does not match between response and covariates!')
    X <- as.matrix(X)
  }
  if (!is.null(U)) {
    if (ncol(counts) != nrow(U)) stop('Sample size does not match between response and covariates!')
    U <- as.matrix(U)
  }

  nsample <- ncol(counts)
  
  ini.beta.idx <- ini.gamma.idx <- 0
  
  if (!is.null(initial.pars$beta) || !is.null(initial.pars$gamma)) {
    
    if (!is.null(initial.pars$beta)) {
      
      ini.beta <- as.matrix(initial.pars$beta)
      
      if ((nrow(ini.beta) != nrow(counts)) || (ncol(U) != ncol(ini.beta))) {
        stop('Dimension of initial.pars$beta incorrect!')
      } else {
        counts <- cbind(counts, ini.beta) 
        ini.beta.idx <- 1
      }
    }
    
    if (!is.null(initial.pars$gamma)) {
      
      ini.gamma <- as.matrix(initial.pars$gamma)
      
      if ((nrow(ini.gamma) != nrow(counts)) || (ncol(U) != ncol(ini.gamma))) {
        stop('Dimension of initial.pars$gamma incorrect!')
      } else {
        counts <- cbind(counts, ini.gamma) 
        ini.gamma.idx <- 1
      }
    }
  }
  
  if (mc.cores > 1) {
    
    require(parallel)
    
    no_cores <- detectCores() - 1
    if(mc.cores > no_cores) { 
      mc.cores <- no_cores
      if (verbose) message("Only ", mc.cores, " threads can be used!")
    } else {
      if (verbose) message(mc.cores, " threads are using!")  
    }
    
    countsList <- as.list(as.data.frame(t(counts)))
    
    if (.Platform$OS.type == "unix") {
      
      Dat <- mclapply(countsList, ftestzi, alpha, ini.beta.idx, ini.gamma.idx, nsample, X, U, contrast, offsets,
                      mc.cores = mc.cores)
      
    } else {
      type <- if (exists("mcfork", mode="function")) "FORK" else "PSOCK"
      cl <- makeCluster(mc.cores, type=type)
      setDefaultCluster(cl)
      clusterEvalQ(NULL, library(MDSeq))
      clusterEvalQ(NULL, library(quadprog))
      Dat <- parLapply(NULL, countsList, ftestzi, alpha, ini.beta.idx, ini.gamma.idx, nsample, X, U, contrast, offsets)
      stopCluster(cl)
      
    }
    
    Dat <- data.frame(matrix(unlist(Dat), nrow=nrow(counts), byrow=T))
    
  } else {
    
    if (verbose) message("Only using 1 thread!") 
    
    Dat <- as.data.frame(t(apply(counts, 1, ftestzi, alpha, ini.beta.idx, ini.gamma.idx, nsample, X, U, contrast, offsets)))
    
  }
  
  rownames(Dat) <- rownames(counts)
  colnames(Dat) <- c('ZI', 'LR', 'pval')
  
  totalTime <- round((proc.time() - startTime)[3], 2)
  if (verbose) message("Total time elapsed:",totalTime, " seconds")
  
  return(Dat)
  
}




ftestzi <- function(x, alpha, ini.beta.idx, ini.gamma.idx, nsample, X, U, contrast, offsets, ...)
{
  x <- as.integer(x)
  count <- x[1:nsample]
  ans <- pval <- dev <- NA
  
  if (all(count[!is.na(count)] > 0)) {
    
    ans <- 0
    dev <- 0
    pval <- 1
    
  } else {
    ibeta <- igamma <- NULL
    
    if (ini.beta.idx || ini.gamma.idx) {
      if (ini.beta.idx) {
        ibeta <- x[(nsample+1):(nsample+ncol(X)+1)]
      }
      if (ini.gamma.idx) {
        igamma <- x[(length(x)-ncol(U)+1):length(x)]
      }
    }
    
    suppressWarnings(fit <- try(glm.ZIMD(count, X, U, contrast, initial.pars = list(beta=ibeta, gamma=igamma), offsets, ZI=TRUE), ...))
    
    suppressWarnings(fit0 <- try(glm.ZIMD(count, X, U, contrast, initial.pars = list(beta=ibeta, gamma=igamma), offsets, ZI=FALSE), ...))
    
    if (inherits(fit,'try-error') || inherits(fit0,'try-error')) {
      
      ans <- 2
      
    } else {
      
      loglik <- fit$loglik
      loglik0 <- fit0$loglik
      
      if (is.na(loglik) || is.na(loglik0)) {
        
        ans <- 2
        
      } else {
        
        dev <- 2*abs(loglik - loglik0)
        pval <- pchisq(dev, df = 1, lower.tail = F)
        if (pval < alpha) {
          ans <- 1
        } else {
          ans <- 0
        }
      }
    }
  }
  return(c(ZI=as.integer(ans), LR=dev, pval=pval))
}
