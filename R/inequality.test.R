#' Inequality test (two factors)
#' 
#' This function will perform an inequality test that evaluates log fold-changes beyond a given threshold level. The null hypothesis 
#' is "|log2FC| is less or equal than given threshold"; and the 
#' alternative is "|log2FC| is greater than given threshold".  
#' 
#' @param est estimated coefficients.
#' @param se standard error of estimated coefficients.
#' @param log2FC.threshold non-negative value which specifies 
#'                          a log2 fold change threshold, default is 1.
#' 
#' @return test statistics and p-value for inequality test
#' 
#' @export
#' 


inequality.test <- function(est, se, log2FC.threshold=1)
{
  # convert log2FC.threshold to natural log scale
  logFC.threshold <- log2FC.threshold * log(2)
  stats <- p.value <- NA
  if (!is.na(est) && !is.na(se)) {
    if (abs(est) <= logFC.threshold) {
      p.value <- 1
      stats <- 0
    } else {
      if (est>0) {
        stats <- (est - logFC.threshold)^2 / se^2
        p.value <- 0.5*pchisq(stats, df=1, lower.tail = F)
      }else{
        stats <- (est + logFC.threshold)^2 / se^2
        p.value <- 0.5*pchisq(stats, df=1, lower.tail = F)
      }   
    } 
  }
  return(c(stats=stats, pval=p.value))
}
