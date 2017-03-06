#' Filtering out lowly expressed genes
#' 
#' @param counts matrix of raw counts.
#' @param mean.cpm.cutoff cutoff point for mean counts per million (cpm).
#' @param num.per.group If num.per.group > 0, gene will be filtered by
#'                      number of samples with cpm > threshold per group. If 
#'                      num.per.group = 4, that means gene which's cpm > threshold 
#'                      in at least 4 samples will be kept. Note that  
#'                      mean.cpm.cutoff is unless when num.per.group > 0. The 
#'                      default is num.per.group equals to 0.
#' @param threshold.per.group given threshold for cpm per group. Default is 1.
#' @param group vector of factors containing group/treatment information
#' @param lib.sizes library size for each sample. If not given, it will be 
#'                  estimated by the column sum of each sample.
#' @param remove.zeros logical, if True (default), 
#'                      non-expressed genes will be removed first.
#' 
#' @return filtered counts matrix
#' 
#' @examples 
#' library(MDSeq)
#' data(sampleData)
#' 
#' dat <- sample.exprs
#' dim(dat) 
#' 
#' # group information
#' group <- sample.pheno$group
#'  
#' # lowly expressed gene filtered by mean cpm value across all samples
#' dat.filtered <- filter.counts(dat, mean.cpm.cutoff = 0.1)  
#' dim(dat.filtered)
#' 
#' # lowly expressed gene filtered by the number of samples which are 
#' # equal to or greather than given threshold in each treatment/group
#' dat.filtered <- filter.counts(dat, num.per.group = 4, threshold.per.group = 1, group = group)  
#' dim(dat.filtered)
#' 
#' @export
#' 



filter.counts <- function(counts, mean.cpm.cutoff=0.05,  num.per.group = 0, threshold.per.group = 1, 
                          group = NULL, lib.sizes = NULL, remove.zeros = TRUE)
{
  require(edgeR)
  
  if (remove.zeros) {
    counts <- counts[rowSums(counts) > 0, ]   
  }
  
  if (is.null(lib.sizes)) {
    lib.sizes <- colSums(counts)
  }
  
  y <- edgeR::DGEList(counts = counts, lib.size = lib.sizes)
  cpmy <- cpm(y)
  
  if (num.per.group == 0) {
    isexpr <- rowMeans(cpmy) > mean.cpm.cutoff
  } else {
    if (is.null(group)) {
      stop("Group information is not provided.")
    } else {
      unq <- unique(group)
      idxm <- matrix(NA, nrow=nrow(counts), ncol=length(unq))
      for(i in 1:length(unq)) {
        idxm[,i] <- rowSums(cpmy[,group==unq[i]] > threshold.per.group) >= num.per.group
      }
      isexpr <- rowSums(idxm) == length(unq)
    }
  }
  
  
  fcounts <- y[isexpr, , keep.lib.sizes = FALSE]
  
  return(fcounts$counts)
  
}