#' RNA-Seq sample data
#' 
#' A RNA-Seq sample data which contains raw counts of expression 
#' and phenotype information. There are total 200 samples, for which one half 
#' of them are in the case group, and the other half are in the control group. 
#' In the raw count matrix, each of 200 samples has RNA-Seq counts measured 
#' from 50000 genes. Two phenotype information were kept 
#' as X1 and X2 for each sample.
#' 
#' @format two data frames contain raw count and phenotype information, respectively.
#' \describe{
#'  \item{sample.exprs}{raw count matrix, 50000X200}
#'  \item{sample.pheno}{phenotype information. covariates, 'X1' and 'X2'; 
#'                       treatment assignment, 'group'}
#' }
#' 
#' @name sampleData
#' @docType data
#' 
NULL