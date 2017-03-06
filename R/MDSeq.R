#' Differential mean and variability analysis with mean-dispersion GLMs.
#' 
#' This is a parallel version of glm.ZIMD.
#' 
#' @param count matrix of counts, each row containing expression values for one gene.
#' @param X vector or matrix contains covariates in the mean GLM.
#' @param U vector or matrix contains covariates in the dispersion GLM.
#' @param contrast list containing contrast for mean and dispersion.
#' @param initial.pars list giving initial values for parameters of mean and/or dispersion for each gene.
#' @param offsets numeric vector of same size as rows of counts giving offsets for use in log-linear models.
#' @param test.ZI whether test zero-inflation due to technical excess zeros.
#' @param ZI.alpha specific cutoff point indicating significance of zero-inflation.
#' @param ZI vector of logical values indicating zero-inflation for each gene. 
#'            If True (defualt), the effect of zero-inflation will be considered. 
#'            If ZI=False, the effect of zero-inflation will be ignored in computation. 
#' @param verbose logical, if True, progress will be reported.
#' @param mc.cores integer equals to or greater than 1 (default), indicating how many cores 
#'                  will be used in computation. If mc.cores > 0, then parallel processes 
#'                  will be used, based on the 'parallel' package.
#' @param ... other arguments could be passed in.
#' 
#' @return a class 'ZIMD' object containing testing results and other information:
#' \item{Dat}{matrix containing testing results. See function 'glm.ZIMD'.}
#' \item{X}{covariate in the mean GLM.}
#' \item{U}{covariate in the dispersion GLM.}
#' \item{counts}{counts matrix of gene expression.}
#' \item{contrast}{list containing contrast matrix used in the model.}
#' \item{test.ZI}{if 1, zero-inflation will be tested using likelihood ratio test.}
#' \item{mc.cores}{how many cores were used in computation.}
#' \item{totalTime}{total time elapsed in computation.}
#' 
#' @seealso testing for single gene is done by glm.ZIMD.
#' 
#' @examples 
#' library(MDSeq)
#' data(sampleData)
#' 
#' #expression
#' dat <- sample.exprs
#' dim(dat)
#' 
#' # covariates
#' X <- sample.pheno[,c("X1","X2")]
#' # group information
#' group <- sample.pheno$group
#' 
#' # lowly expressed gene filtered by mean cpm value across all samples
#' dat.filtered <- filter.counts(dat, mean.cpm.cutoff = 0.1) 
#' dim(dat.filtered)
#' 
#' # design matrix and constrast setting
#' group <- factor(sample.pheno$group, labels = c("Control", "Case"))
#' # make design matrix with proper contrast setting
#' groups <- get.model.matrix(group)
#' 
#' # normalization factor was calculated by using TMM method
#' # and used as an offset in mean-disperison GLM
#' # require(edgeR)
#' cnf <- calcNormFactors(dat.filtered, method="TMM") 
#' libsize <- colSums(dat.filtered)              #normalization factor
#' rellibsize <- libsize/exp(mean(log(libsize))) #relative library size
#' nf <- cnf * rellibsize                        #final normalization factor including library size
#'
#' # check outliers using parallel process with 4 threads
#' # the first 100 genes are used as an example
#' dat.checked <- remove.outlier(dat.filtered[1:100,], X=X, U=X, contrast = groups,
#'                               offsets = nf, mc.cores = 4)
#' 
#' # status of outlier checking
#' table(dat.checked$outlier$status)
#' 
#' # frequency distribtuion of outliers
#' table(dat.checked$outlier$num.outliers)
#' 
#' # remove genes with status flag other than 0
#' counts <- dat.checked$count[dat.checked$outliers$status==0,]
#' dim(counts)
#' 
#' # differential analysis
#' # using parallel process with 4 threads
#' fit <- MDSeq(counts, X=X, U=X, contrast = groups, offsets = nf, mc.cores = 4)
#' 
#' # testing with a given log2-fold-change threshold tau
#' # Ha: |log2FC| > tau
#' # tau = 1 
#' result <- extract.ZIMD(fit, compare = list(A="Case", B="Control"), log2FC.threshold = 1)
#' head(result)
#' 
#' @export
#' 



MDSeq <- function(count, X = NULL, U = NULL, contrast = list(mean=NULL,dispersion=NULL),
                  initial.pars = list(beta = NULL, gamma = NULL), offsets = NULL,
                  test.ZI=TRUE, ZI.alpha=0.05, ZI = rep(TRUE, nrow(counts)), 
                  verbose = TRUE, mc.cores = 1, ...)
{
  startTime <- proc.time()
  
  counts <- as.data.frame(count)
  nsample <- ncol(counts)
  
  varName <- extractName(counts, X, U, contrast)
  if(is.na(varName[1])) {
    stop('Please check inputs!')
  } else {
    nout <- length(varName)
  }
  
  if (!is.null(X)) {
    if (ncol(counts) != nrow(X)) stop('Sample size does not match between response and covariates!')
    X <- as.matrix(X)
  }
  if (!is.null(U)) {
    if (ncol(counts) != nrow(U)) stop('Sample size does not match between response and covariates!')
    U <- as.matrix(U)
  }
  
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
  
  if (test.ZI==FALSE) {
    if(nrow(counts) != length(ZI)) stop('Dimension of ZI incorrect!')
    counts <- cbind(counts, ZI)
  }
  
  if (mc.cores > 1) {
    
    require(parallel)
    
    no_cores <- detectCores() - 1
    if(mc.cores > no_cores) { 
      mc.cores <- no_cores
      if (verbose) message("Only ", mc.cores, " threads can use!")
    } else {
      if (verbose) message(mc.cores, " threads are using!")  
    }
    
    countsList <- as.list(as.data.frame(t(counts)))
    
    if (.Platform$OS.type == "unix") {
      
      Dat <- mclapply(countsList, ff, X, U, contrast, offsets, test.ZI, ZI.alpha, nsample, nout=nout, ini.beta.idx, ini.gamma.idx, 
                      mc.cores = mc.cores)
      
    } else {
      
      type <- if (exists("mcfork", mode="function")) "FORK" else "PSOCK"
      cl <- makeCluster(mc.cores, type=type)
      setDefaultCluster(cl)
      clusterEvalQ(NULL, library(MDSeq))
      clusterEvalQ(NULL, library(quadprog))
      Dat <- parLapply(NULL, countsList, ff, X, U, contrast, offsets, test.ZI, ZI.alpha, nsample, nout, ini.beta.idx, ini.gamma.idx)
      stopCluster(cl)
      
    }
    
    Dat <- data.frame(matrix(unlist(Dat), nrow=nrow(counts), byrow=T))
    
  } else {
    
    if (verbose) message("Only 1 thread is using!") 
    
    Dat <- as.data.frame(t(apply(counts, 1, ff, X, U, contrast, offsets, test.ZI, ZI.alpha, nsample, nout, ini.beta.idx, ini.gamma.idx)))
    
  }
  
  rownames(Dat) <- rownames(counts)
  
  varName <- extractName(count, X, U, contrast)
  
  if (!is.na(varName[1])) {
    if (test.ZI==FALSE) {
      colnames(Dat) <- c(varName, 's', 'se.s', 'loglik', 'convergence', 'EM.convergence')
    } else {
      colnames(Dat) <- c(varName, 's', 'se.s', 'loglik', 'convergence', 'EM.convergence', 'ZI', 'ZI.LR', 'ZI.pval')
    }
  }
  
  totalTime <- round((proc.time() - startTime)[3], 2)
  if (verbose) message("Total time elapsed:",totalTime, " seconds")
  
  rslt <- list(Dat=Dat,
               X=X,
               U=U,
               counts=counts,
               contrast=contrast,
               test.ZI=test.ZI,
               mc.cores=mc.cores,
               totalTime=totalTime)
  
  class(rslt) <- "ZIMD"
  
  return(rslt)
  
}







ff <- function(x, X, U, contrast, offsets, test.ZI, ZI.alpha, nsample, nout, ini.beta.idx, ini.gamma.idx, ...)
{
  #out <- rep(NA, nout)
  
  if (test.ZI==FALSE) {
    y <- x[-length(x)]
    ZI <- x[length(x)]
  } else {
    y <- x
    ans <- pval <- dev <- NA
  }
  
  counts <- y[1:nsample]
  out <- rep(NA,nout)
  
  ibeta <- igamma <- NULL
  
  if (ini.beta.idx || ini.gamma.idx) {
    if (ini.beta.idx) {
      ibeta <- y[(nsample+1):(nsample+ncol(X)+1)]
    }
    if (ini.gamma.idx) {
      igamma <- y[(length(y)-ncol(U)+1):length(y)]
    }
  }
  
  if (test.ZI==FALSE) {
    
    if (ZI %in% c(0,1)) {
      
      suppressWarnings(fit <- try(glm.ZIMD(counts, X, U, contrast, initial.pars = list(beta=ibeta, gamma=igamma), offsets, ZI=ZI,...)))
      
      if(!inherits(fit,'try-error')) {
        out <- unlist(c(fit$coefficients,fit$se,fit$coefficients.contrast,fit$s,loglik=fit$loglik,convergence=fit$convergence,EM.convergence=fit$EM.convergence))
      }
      
    }
    if(length(out)!=(nout+5)) { out <- rep(NA, (nout+5))}
    
  } else {
    
    if (all(counts[!is.na(counts)] > 0)) {
      
      ans <- 0; dev = 0; pval = 1;
      
      suppressWarnings(fit <- try(glm.ZIMD(counts, X, U, contrast, initial.pars = list(beta=ibeta, gamma=igamma), offsets,...)))
      
      if (!inherits(fit,'try-error')) {
        out <- c(unlist(c(fit$coefficients,fit$se, fit$coefficients.contrast,fit$s,loglik=fit$loglik,
                          convergence=fit$convergence,EM.convergence=fit$EM.convergence)),
                 ZI=ans, LR=dev, pval=pval)
      } else {
        out <- c(out, ans, dev, pval)
      }
      
    } else {
      
      suppressWarnings(fit <- try(glm.ZIMD(counts, X, U, contrast, initial.pars = list(beta=ibeta, gamma=igamma), offsets, ZI=TRUE,...)))
      
      suppressWarnings(fit0 <- try(glm.ZIMD(counts, X, U, contrast, initial.pars = list(beta=ibeta, gamma=igamma), offsets, ZI=FALSE,...)))
      
      if (inherits(fit,'try-error') || inherits(fit0,'try-error')) {
        
        ans <- 2
        out <- c(out, ans, dev, pval)
        
      } else {
        
        loglik <- fit$loglik
        loglik0 <- fit0$loglik
        
        if (is.na(loglik) || is.na(loglik0)) {
          
          ans <- 2
          out <- c(out, ans, dev, pval)
          
        } else {
          
          dev <- 2*abs(loglik - loglik0)
          pval <- pchisq(dev, df = 1, lower.tail = F)
          if (pval < ZI.alpha) {
            ans <- 1
            fit <- fit
          } else {
            ans <- 0
            fit <- fit0
          }
        }
        
        out <- c(unlist(c(fit$coefficients,fit$se,fit$coefficients.contrast,fit$s,loglik=fit$loglik,
                          convergence=fit$convergence,EM.convergence=fit$EM.convergence)),
                 ZI=ans, LR=dev, pval=pval)
        
      }
    }
    if(length(out)!=(nout+8)) { out <- rep(NA, (nout+8))}
  }
  return(out)
}


extractName <- function(counts, X = NULL, U = NULL, contrast = list(mean=NULL,dispersion=NULL),
                        initial.pars = list(beta = NULL, gamma = NULL))
{
  i = 1
  varName = NA
  while(TRUE){
    suppressWarnings(fitn <- try(glm.ZIMD(counts[i,], X, U, contrast), silent = T))
    if (!inherits(fitn,'try-error')) {
      if (!anyNA(unlist(c(fitn$coefficients,fitn$se)))) {
        varName <- names(unlist(c(fitn$coefficients,fitn$se,fitn$coefficients.contrast)))
        break;
      } else {
        i = i + 1
      }
    } else {
      i = i + 1
    }
    if (i > 100) break;
  }
  return(varName)
}
