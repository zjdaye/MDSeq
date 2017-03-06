#' Creating design matrix from contrast matrices and setting contrasts if not supplied
#' 
#' @param group vector of factors containing group/treatment information
#' @param contrast.matrix list containing predefined contrast matices for both mean and dispersion
#' @param contrast.type if a contrast matix is given, which type of contrast 
#'                      will be used in setting contrasts. Default is sum to zero contrasts, 'contr.sum'. 
#'                      Also, other contrast types can be used, such as 
#'                      "contr.helmert", "contr.poly", "contr.treatment", and "contr.SAS".
#' 
#' @return a list containing design matrices for both mean and dispersion
#' 
#' @examples 
#' group <- factor(rep(1:3,4))
#' get.model.matrix(group)
#' 
#' #treatment contrast - more tranditional glm forming
#' get.model.matrix(group, contrast.type="contr.treatment")
#' 
#' @export
#'    


get.model.matrix <- function(group, contrast.matrix=list(mean=NULL, dispersion=NULL), contrast.type="contr.sum") {
  if (!is.factor(group)) { stop("group is not a vector of factors.")  }
  factors = levels(group); n.factors = length(factors);
  
  if (!is.element(contrast.type, c("contr.helmert", "contr.poly", "contr.sum", "contr.treatment", "contr.SAS"))) 
    { stop("contrast.type is not available in R!") }
  
  # Compute contrast from constrast.type if not given.
  if (is.null(contrast.matrix$mean)) { contrast.matrix$mean = do.call(contrast.type, list(n.factors))
  } else { 
    if (nrow(contrast.matrix$mean) != n.factors) { stop("Rows of contrast.matrix$mean do not match with number of factors.") }
    if (ncol(contrast.matrix$mean) != (n.factors-1)) { stop("Columns of contrast.matrix$mean do not match with number of factors.") }
  }
  contrasts(group) = contrast.matrix$mean
  design.mean = model.matrix(~group)
  attributes(design.mean)$contrast.matrix = attributes(design.mean)$contrast.matrix$group  # no need for group list
  
  if (is.null(contrast.matrix$dispersion)) { contrast.matrix$dispersion = do.call(contrast.type, list(n.factors))
  } else { 
    if (nrow(contrast.matrix$dispersion) != n.factors) { stop("Rows of contrast.matrix$dispersion do not match with number of factors.") }
    if (ncol(contrast.matrix$dispersion) != (n.factors-1)) { stop("Columns of contrast.matrix$dispersion do not match with number of factors.") }
  }
  
  contrasts(group) = contrast.matrix$dispersion
  design.dispersion = model.matrix(~group)
  attributes(design.dispersion)$contrast.matrix = attributes(design.dispersion)$contrast.matrix$group # no need for group list
  
  return(design.matrix=list(mean=design.mean, dispersion=design.dispersion))
}
