#' Generating random samples (genes) from ZIMD distribution
#' 
#' @param n number of total observations in each simulated sample.
#' @param nsim number of simulated samples.
#' @param contrast list containing contrast matrices 
#'                  for mean and dispersion, respectively.
#' @param beta preset coefficients in the mean GLM
#' @param gamma preset coefficients in the dispersion GLM
#' @param beta.size number of trials (zero or more)
#' @param gamma.size number of trials (zero or more)
#' @param beta.prob probability of success on each trial
#' @param gamma.prob probability of success on each trial
#' @param s probability of technical excess zeros
#' @param tol tolerance for simulation computation.
#' @param max_val maximum value for simulated count.
#' @param seed random seed
#' 
#' @return a list containing random samples and 
#'          corresponding covariates for mean and dispersion.
#' \item{y}{matrix containing random samples.}
#' \item{X}{vector or matrix containing covariates for mean.}
#' \item{U}{vector or matrix containing covariates for dispersion.}
#' 
#' @examples 
#' # simulate 1000 genes, each has 200 observations.
#' n <- 200; nsim <- 1000; seed <- 123
#' 
#' # set covariate coefficient for mean and dispersion
#' beta <- c(5, 1, 0.5, 0.5); gamma <- c(4, 0, 0.5, 0.5)
#' 
#' # set proportion of zero-inflation
#' s <- 0
#' 
#' # set design matrix for mean and dispersion 
#' factors <- c('Control', 'Case')
#' group <- gl(n=2, k=n/2, labels=factors)
#' design <- get.model.matrix(group, contrast.type="contr.treatment")
#' 
#' # simulate data nsim (genes) X n (observations)
#' ysim <- sim.ZIMD(n=n, nsim=nsim, beta=beta, gamma=gamma, contrast = design, s=s, seed=seed)
#' 
#' @export
#'  


sim.ZIMD <- function(n=200, nsim=1000,
                    beta, beta.size=2, beta.prob=c(.2,.5),
                    gamma, gamma.size=2, gamma.prob=c(.2,.5),
                    contrast, s, tol=1e-12, max_val=1e300, seed=NA)
{
  if(!is.na(seed)) { set.seed(seed) }
  
  p <- length(beta)
  q <- length(gamma)
  
  X <- matrix(rbinom(p*n, beta.size, prob=beta.prob), n, p)
  U <- matrix(rbinom(q*n, gamma.size, prob=gamma.prob), n, q)
  
  X[,1:dim(contrast$mean)[2]] <- contrast$mean
  U[,1:dim(contrast$dispersion)[2]] <- contrast$dispersion
  
  mu0 <- pmin( exp(X%*%beta), max_val)
  phi0 <- pmax( exp(U%*%gamma), 1+tol)
  theta <- mu0/(phi0-1)
  
  y <- matrix(NA, nrow=nsim, ncol=n)
  
  for (j in 1:nsim) {
    for(i in 1:n) { y[j, i] <- rnbinom(n=1, size=theta[i], mu=mu0[i]) }
    y[j, sample(1:n, floor(n*s))] <- 0 
  }
  
  rownames(y) <- paste('gene', 1:nsim, sep='')
  colnames(y) <- paste('obs', 1:n, sep='')
  
  return(list(y=y,
              X=X,
              U=U))
}