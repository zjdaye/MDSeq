#' A glmZIMD class
#' 
#' @export
#' 

glmZIMD <- function(y, ...) UseMethod("glmZIMD")


#' @export
#'
print.glmZIMD <- function(x, ...)
{
  cat("Call:\n")
  print(x$call)
  
  beta <- as.matrix(cbind(Estimate=x$coefficients$beta, "Std. Error"=x$se$se.beta, 
                          "z value"=x$zvalue$zvalue.beta, "Pr(>|z|)"=x$pvalue$pvalue.beta))
  gamma <- as.matrix(cbind(Estimate=x$coefficients$gamma, "Std. Error"=x$se$se.gamma, 
                           "z value"=x$zvalue$zvalue.gamma, "Pr(>|z|)"=x$pvalue$pvalue.gamma))
  coefficients <- list(beta=beta, gamma=gamma)
  
  cat("\n--------------\n")
  cat("Coefficients:\n")
  cat("--------------\n")
  cat("Mean:\n")
  print(coefficients$beta, digits=4)
  cat("\nDispersion:\n")
  print(coefficients$gamma, digits=4)
  
  cat('\n')
  
  cat('Zero-inflation parameter:\n')
  s <- cbind(Estimation=x$s$s, Std.Error=x$s$se.s)
  rownames(s) <- '                      '
  print(s, digits = 4)
}