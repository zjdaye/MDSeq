#' A ZIMD class
#' 
#' @export
#' 

ZIMD <- function(y, ...) UseMethod("ZIMD")


#' Printing a running summary for a ZIMD object 
#'
#' @export
#' 

print.ZIMD <- function(obj, ...)
{
  if (obj$test.ZI) {
    cat("Inputs:", ncol(obj$counts), "observations with total",nrow(obj$counts),"genes.\n", sep=' ')
    cat("Test Zero-inflation for each gene:", "True.\n")
  } else {
    cat("Inputs:", ncol(obj$counts)-1, "observations with total",nrow(obj$counts),"genes.\n", sep=' ')
    cat("Test Zero-inflation for each gene:", "False.\n")
  }
  
  if (obj$mc.cores == 1){
    cat("Only single thread was used.\n")
  } else {
    cat("Request for multithreads:", obj$mc.cores, "threads were used.\n")
  }
  
  cat("Total time elapsed (seconds):", obj$totalTime,"\n")

}


#' Extracting result from a ZIMD object (two factor test)
#' 
#' This program will extract testing results and can simultaneously perform 
#' a two-sided Wald's test whether log2FC = 0 or an inequality test, in which 
#' the null hypothesis is that "|log2FC| is less or equal than a given threshold"; 
#' wehreas the alternative is "|log2FC| is greater than threshold".  
#' 
#' @param obj a ZIMD object containing test statistics.
#' @param get which one to extract, 'contrast' or a specific covariate.
#' @param compare if get='contrast', a comparison of group A versus group B must be specified. 
#'                A and B have to be two of group names in the contrast matrix.
#' @param log2FC.threshold log2 fold change threshold for inequality testing. Default is 0. 
#'                          That indicates a test that the log2 fold change 
#'                          is equal to 0. If log2FC.threshold is specified not equal to 0, 
#'                          an inequality test will be performed for |log2FC|>threshold.
#' @param p.adj p-value adjustment method, default is 'BY' correction.
#' @param Log2 logical, if True, difference in mean and dispersion are given in log2 scale.

#'                
#' @return a matrix containing testing results.
#' \item{log2FC.mean}{log2 fold change for mean}
#' \item{statistic.mean}{test statistic for mean}
#' \item{pvalue.mean}{p-value for mean}
#' \item{FDR.mean}{FDR after correction for mean}
#' \item{log2FC.dispersion}{log2 fold change for dispersion}
#' \item{statistic.dispersion}{test statistic for dispersion}
#' \item{pvalue.dispersion}{p-value for dispersion}
#' \item{FDR.dispersion}{FDR after correction for dispersion}
#' 
#' @seealso inequality test is done by inequality.test.
#' 
#' 
#' @export
#'

extract.ZIMD <- function(obj, get = 'contrast', compare=list(A=NULL, B=NULL), 
                         log2FC.threshold = 0, p.adj = 'BY', Log2 = TRUE, ...)
{
  Dat <- obj$Dat
  
  result <- NULL
  
  if (is.null(compare$A) || is.null(compare$B)) {
    stop("Please specify a single comparison as 'A vs B'.")
  }
  
  if (get == "contrast") {
    if (is.null(obj$contrast$mean) && is.null(obj$contrast$dispersion)) {stop('This result does not contain \'contrast\' term!')}
    AB <- paste(c(compare$A, compare$B), collapse = 'vs')
    betaAB <- paste(c('beta',AB), collapse = '.')
    gammaAB <- paste(c('gamma',AB), collapse = '.')
    if ( (betaAB %in% colnames(Dat)) || (gammaAB %in% colnames(Dat)) ) {
      
      if (!is.null(obj$contrast$mean) && (betaAB %in% colnames(Dat)) ) {
        beta.mean <- as.vector(Dat[,betaAB])
        se.beta.mean <- as.vector(Dat[,paste(c('beta.se',AB),collapse = '.')])
        mat.mean <- cbind(beta.mean, se.beta.mean)
        if (log2FC.threshold != 0) {
          pval.mean <- t(apply(mat.mean, 1, function(x) inequality.test(x[1], x[2], log2FC.threshold=log2FC.threshold)))
        } else {
          pval.mean <- cbind(beta.mean/se.beta.mean, 2*pnorm(-abs(beta.mean/se.beta.mean)))
        }
        fdr.mean <- p.adjust(pval.mean[,2], method = p.adj)
        #log2
        if (Log2) {
          result <- cbind(result, beta.mean/log(2), pval.mean, fdr.mean)
          colnames(result) <- c(paste(AB,'mean.log2FC',log2FC.threshold,sep='.'), "Statistics.mean", "Pvalue.mean", "FDR.mean")
        } else {
          result <- cbind(result, beta.mean, pval.mean, fdr.mean)
          colnames(result) <- c(paste(AB,'mean.log2FC',log2FC.threshold,sep='.'), "Statistics.mean", "Pvalue.mean", "FDR.mean")
        }
      }
      if (!is.null(obj$contrast$dispersion) && (gammaAB %in% colnames(Dat))) {
        gamma.dispersion <- as.vector(Dat[,gammaAB])
        se.gamma.dispersion <- as.vector(Dat[,paste(c('gamma.se',AB),collapse = '.')])
        mat.dispersion <- cbind(gamma.dispersion, se.gamma.dispersion)
        if (log2FC.threshold != 0) {
          pval.dispersion <- t(apply(mat.dispersion, 1, function(x) inequality.test(x[1], x[2], log2FC.threshold=log2FC.threshold))) 
        } else {
          pval.dispersion <- cbind(gamma.dispersion/se.gamma.dispersion, 2*pnorm(-abs(gamma.dispersion/se.gamma.dispersion)))
        }
        fdr.dispersion <- p.adjust(pval.dispersion[,2], method = p.adj)
        #log2
        if (Log2) {
          result.d <- cbind( gamma.dispersion/log(2), pval.dispersion, fdr.dispersion)
          colnames(result.d) <- c(paste(AB,'dispersion.log2FC',log2FC.threshold,sep='.'), "Statistics.dispersion", "Pvalue.dispersion", "FDR.dispersion")
        }else{
          result.d <- cbind( gamma.dispersion, pval.dispersion, fdr.dispersion)
          colnames(result.d) <- c(paste(AB,'dispersion.log2FC',log2FC.threshold,sep='.'), "Statistics.dispersion", "Pvalue.dispersion", "FDR.dispersion")
        }
        result <- cbind(result, result.d)
      }
    }
  } else {
    beta.cov <- paste("beta", get, sep = '.')
    gamma.cov <- paste("gamma", get, sep = '.')
    if (!(beta.cov %in% colnames(Dat)) && !(gamma.cov %in% colnames(Dat))) {
      stop("Could not find this covariate! Please select \'contrast\' or other valid covariate!")
    } else {
      if (beta.cov %in% colnames(Dat)) {
        beta.mean <- as.vector(Dat[,beta.cov])
        se.beta.mean <- as.vector(Dat[, paste('se',beta.cov,sep = '.')])
        mat.mean <- cbind(beta.mean, se.beta.mean)
        if (log2FC.threshold != 0) {
          pval.mean <- t(apply(mat.mean, 1, function(x) inequality.test(x[1], x[2], log2FC.threshold=log2FC.threshold))) 
        } else {
          pval.mean <- cbind(beta.mean/se.beta.mean, 2*pnorm(-abs(beta.mean/se.beta.mean)))
        }
        fdr.mean <- p.adjust(pval.mean[,2], method = p.adj)
        #log2
        if (Log2) {
          result <- cbind(result, beta.mean/log(2), pval.mean, fdr.mean)
          colnames(result) <- c(paste(get,'mean','log2FC',log2FC.threshold,sep='.'), "Statistics.mean", "Pvalue.mean", "FDR.mean")
        } else {
          result <- cbind(result, beta.mean, pval.mean, fdr.mean)
          colnames(result) <- c(paste(get,'mean','log2FC',log2FC.threshold,sep='.'), "Statistics.mean", "Pvalue.mean", "FDR.mean")
        }
      }
      if (gamma.cov %in% colnames(Dat)) {
        gamma.dispersion <- as.vector(Dat[,gamma.cov])
        se.gamma.dispersion <- as.vector(Dat[, paste('se',gamma.cov,sep = '.')])
        mat.dispersion <- cbind(gamma.dispersion, se.gamma.dispersion)
        if (log2FC.threshold != 0) {
          pval.dispersion <- t(apply(mat.dispersion, 1, function(x) inequality.test(x[1], x[2], log2FC.threshold=log2FC.threshold))) 
        } else {
          pval.dispersion <- cbind(gamma.dispersion/se.gamma.dispersion, 2*pnorm(-abs(gamma.dispersion/se.gamma.dispersion)))
        }
        fdr.dispersion <- p.adjust(pval.dispersion[,2], method = p.adj)
        #log2
        if (Log2) {
          result.d <- cbind( gamma.dispersion/log(2), pval.dispersion, fdr.dispersion)
          colnames(result.d) <- c(paste(get,'dispersion','log2FC',log2FC.threshold,sep='.'), "Statistics.dispersion", "Pvalue.dispersion", "FDR.dispersion")
        } else {
          result.d <- cbind( gamma.dispersion, pval.dispersion, fdr.dispersion)
          colnames(result.d) <- c(paste(get,'dispersion','logFC',log2FC.threshold,sep='.'), "Statistics.dispersion", "Pvalue.dispersion", "FDR.dispersion")
        }
        result <- cbind(result, result.d)
      }
    }
  }
  rownames(result) <- rownames(Dat)
  return(as.data.frame(result))
}
