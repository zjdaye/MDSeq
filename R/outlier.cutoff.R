#' Finding cutoff point for outlier checking
#' 
#' The cutoff points of outlier checking will be found
#'  from a variance-gamma distribution.
#' 
#' @param alpha specific significance level. Default is alpha=0.05.
#' 
#' @return a cutoff value corresponding to specific significance level.
#' 
#' @examples 
#' # cutoff point at alpha=0.05
#' cf <- outlier.cutoff(alpha=0.05) 
#' 
#' @export
#' 

outlier.cutoff <- function(alpha=0.05){
  # save several commonly used cutoff values  
  reserved <- c(3.190205, 4.363898, 7.208451)
  ord <- c(0.1, 0.05, 0.01)
  
  if (alpha %in% ord) {
    return(reserved[ord==alpha])
  } else {
    require(VarianceGamma)
    return(2*qvg(1-alpha/2,vgC=0,sigma=1,theta=0,nu=2))
  }
} 
