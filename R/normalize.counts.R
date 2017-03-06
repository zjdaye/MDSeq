#' Normalization to count data
#'
#' Normalizing count data to eliminate composition biases between libraries.
#' 
#' @param counts matrix or data frame containing the count data.
#' @param group vector giving the experimental treatment/condition for each sample (optional).
#' @param lib.sizes vector giving library size of each sample. If lib.sizes=NULL, which is 
#'                  default, then library size will be estimated as the column sum of each 
#'                  sample.
#' @param method specific normalization method. The optional methods are 'TMM' (Robinson and Oshlack, 2010), 
#'                "RLE" (Anders and Huber, 2010), "upperquartile" (Bullard et al, 2010), 
#'                "cqn" (Hansen, Irizarry, and Wu, 2012). These normalization methods are applied using 
#'                the edgeR package, except "cqn" which uses the cqn package.
#' @param verbose logical, if True, progress will be reported.
#' @param cqn.annotation optional, gene information used in cqn method, 
#'                       containing GC contents and gene length.
#' @param ... other arguments could be passed to normalization methods.
#' 
#' @return a numeric matrix containing normalized count data.
#' 
#' @examples 
#' library(MDSeq)
#' data(sampleData)
#' 
#' dat <- sample.exprs
#' dim(dat)
#' 
#' # Lowly expressed gene filtered by mean cpm value
#' dat.filtered <- filter.counts(dat, mean.cpm.cutoff = 0.1)
#' dim(dat.filtered)
#' 
#' # normalize by using TMM method   
#' dat.normalized <- normalize.counts(dat.filtered, method="TMM") 
#' 
#' @export
#' 


normalize.counts <- function(counts, group = NULL, lib.sizes = NULL,
                              method = c("TMM", "RLE", "upperquartile", "cqn"), 
                              verbose = TRUE, cqn.annotation = NULL, ...)
{
  method <- match.arg(method)
  
  if (is.null(lib.sizes)) {
    lib.sizes <- colSums(counts)
  }
  
  if (verbose) {
    message('Using ', method, ' normaliztion.')
  }
  
  if (method != 'cqn') {
    require(edgeR)
    
    y <- DGEList(counts=counts, group=group, lib.size=lib.sizes)
    
    if (verbose) {
      message('Calculating normalization(scaling) factors with the ', method, ' method.')
    }    
    y <- suppressWarnings(calcNormFactors(y, method = method))
    
    if (verbose) {
      message('Adjusting counts to rescaling factors.')
    }
    counts.adj <- ceiling(t(t(counts) * (1/y$samples$norm.factors)) - 0.5)
    
  } else {
    require(cqn)
    
    if (!is.matrix(cqn.annotation) && !is.data.frame(cqn.annotation) || is.null(dim(cqn.annotation)) || dim(cqn.annotation)[2] != 2)
      stop("cqn requires a two-column matrix or data frame in the argument 'cqn.annotation'.")
    
    common <- intersect(rownames(counts), rownames(cqn.annotation[!is.na(cqn.annotation[, 1]) & !is.na(cqn.annotation[, 2]), ]))
    if (length(common) == 0)
      stop("No row names in 'counts' match any row name in 'cqn.annotation'.")
    
    cqnNorm <- cqn::cqn(counts[common, ],
                        lengths=as.numeric(cqn.annotation[common, 1]),
                        x=as.numeric(cqn.annotation[common, 2]),
                        sizeFactors=lib.sizes, verbose=verbose)
    
    log2r <- sweep(cqnNorm$y + cqnNorm$offset, 2, log2(lib.sizes/1e6), FUN="+")
    
    counts.adj <- ceiling(2^log2r-0.5)
    
  }
  
  return(counts.adj)
  
}
