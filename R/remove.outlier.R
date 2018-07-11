#' Removing outliers from count data
#' 
#' Outlier influential upon the effects of specific covariate/covariates will be removed 
#' from a matix of count data.
#' 
#' @param counts matrix or data frame containing gene expression counts data
#' @param X vector or matrix containing covariates in the mean GLM.
#' @param U vector or matrix containing covariates in the dispersion GLM,  
#'          which can be different from X.
#' @param contrast list containing contrast matices for mean and dispersion.
#' @param offsets numeric vector of same size as rows of counts giving offsets for use in log-linear models.
#' @param check which covariate(covariates) to check. if check="contrast", outliers 
#'              influential on the effects of contrasts will be checked; if check="other", outliers 
#'              influential on the effects of all other covariates (except contrasts) will be checked; 
#'              if check="all", outliers influential on the effects of all covariates will be checked. 
#'              In addition, outliers influential towards specific covariate(/covariates) can be checked.
#' @param cut.off cutoff point defining the significance level of outlier selections. 
#' @param ZI logical, if True (defualt), the effect of zero-inflation will be considered.  
#' @param conserv logical, if TRUE (default), outlier detection will skip any gene with estimation 
#'                         error for either coefficients or standard errors of coefficients and label it as an error; 
#'                         if FALSE, outlier detection will only skip any gene with estimation error 
#'                         on coefficients. 
#' @param verbose logical, if True, progress will be reported.
#' @param min.sample.size threshold to keep mininum valid sample size. Default is 50.
#' @param mc.cores integer equals to or greater than 1 (default), indicating how many cores 
#'                  will be used in computation. If mc.cores > 0, then parallel processes 
#'                  will be used, based on 'parallel' package.
#' @param ... other arguments could be passed to outlier checking.
#' 
#' @return a list containing updated count data and summary information.
#' \item{count}{count matrix, in which outliers are set to NA's}
#' \item{outliers}{a matrix containing two columns. 1. status, indicating outlier checking status. 
#'                  if 'status' is not equal to 0, it indicates an error occurred in the checking process.
#'                  2. num.outliers, numbers of outliers found at each gene.}
#'                  
#' @details MDSeq requires sufficient samples and variability within samples to provide variance estimation.  Status=1 error will be produced if sample size minus numbers of outliers is insufficient, and Status=2 error will be produced when variance cannot be estimated, such as when all cases or all controls have zero counts.  The MDSeq does not pool variances across genes in order to allow the interpretation of gene expression variability.
#'                  
#' @seealso outlier checking of single gene is done by check.outlier.
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
#' # the first 100 genes were used as an example
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
#' @export
#' 
#' 


remove.outlier <- function(counts, X = NULL, U = NULL, contrast = list(mean=NULL,dispersion=NULL), offsets=NULL,
                           check = 'contrast', ZI = TRUE, verbose = TRUE, 
                           cut.off=outlier.cutoff(0.05), min.sample.size=50, mc.cores = 1, conserv = TRUE, ...)
{
  startTime <- proc.time()
  
  counts <- as.data.frame(counts)
  
  if (!any(check %in% c(colnames(X), colnames(U), 'contrast'))) 
    stop('Could not find the variable! Please choose \'contrast\' or other valid covariate.')
  
  if (mc.cores > 1) {
    
    require(parallel)
    
    no_cores <- detectCores(all.tests=TRUE) - 1
    if(mc.cores > no_cores) { 
      mc.cores <- no_cores
      if (verbose) message("Only ", mc.cores, " threads can be used!")
    } else {
      if (verbose) message(mc.cores, " threads are using!")  
    }
    
    countsList <- as.list(as.data.frame(t(counts)))
    
    if (.Platform$OS.type == "unix") {
      
      Dat <- mclapply(countsList, fout, X, U, contrast, offsets,
                      check, cut.off, min.sample.size, conserv,
                      mc.cores = mc.cores)
      
    } else {
      type <- if (exists("mcfork", mode="function")) "FORK" else "PSOCK"
      cl <- makeCluster(mc.cores, type=type)
      setDefaultCluster(cl)
      clusterEvalQ(NULL, library(MDSeq))
      clusterEvalQ(NULL, library(quadprog))
      Dat <- parLapply(NULL, countsList, fout, X, U, contrast, offsets,
                       check, cut.off, min.sample.size, conserv)
      stopCluster(cl)
      
    }
    
    Dat <- data.frame(matrix(unlist(Dat), nrow=nrow(counts), byrow=T))
    
  } else {
    
    if (verbose) message("Only using 1 thread!") 
    
    Dat <- as.data.frame(t(apply(counts, 1, fout, X, U, contrast, offsets, check, cut.off, min.sample.size, conserv)))
    
  }
  
  count <- Dat[,1:(ncol(Dat)-2)]
  outlier.idx <- Dat[,(ncol(Dat)-1):ncol(Dat)]
  
  rownames(count) <- rownames(outlier.idx) <- rownames(counts)
  colnames(count) <- colnames(counts)
  colnames(outlier.idx) <- c('status', 'num.outliers')
  
  totalTime <- round((proc.time() - startTime)[3], 2)
  if (verbose) message("Total time elapsed:",totalTime, " seconds")
  
  return(list(count=count, outliers=outlier.idx))
  
}


fout <- function(x, X, U, contrast, offsets, check, cut.off, min.sample.size, conserv, ...)
{
  ans <- as.integer(x)
  
  suppressWarnings(ot <- try(check.outlier(ans, X, U, contrast, offsets, check, cut.off, verbose = F, output=F, conserv = conserv, ...),
                             silent = T))
  
  if (!inherits(ot,'try-error')) {
    if ((length(ans) -length(ot)) > min.sample.size) {
      ans[ot] <- NA
      ans <- c(ans, 0, length(ot))
    } else {
      ans <- c(ans, 1, NA)
    }
  } else {
    ans <- c(ans, 2, NA)
  }
  
  return(as.integer(ans))
  
}
