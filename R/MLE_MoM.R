# functions for getting MLE ---------------------------------------------------

g <- function(N, muHat, sumLog, alpha) {
  #' first derivative of log likelihood
  #'
  #' @param N total sample size
  #' @param muHat sample mean
  #' @param sumLog sum of log(x_i)
  #' @param alpha parameter value
  #' @export

 N * log(alpha / muHat) - N * digamma(alpha) + sumLog
}

h <- function(N, muHat, alpha) {
  #' second derivative of log likelihood
  #'
  #' @param N total sample size
  #' @param muHat sample mean
  #' @param alpha parameter value
  #' @export

  N / alpha - N * trigamma(alpha)
}

gammaFit <- function(z, maxIter = 1000, tol = 1e-4) {
  #' Maximum likelihood and method of moments estimates for iid observations of gamma random variables
  #'
  #' This function computes maximum likelihood and method of moments estimates
  #' for iid observations of gamma random variables.
  #' Maximum likelihood estimatation is via a Newton-Raphson routine.
  #' @param z data
  #' @param maxIter maximum number of iterations (maxIter = 1000 by default)
  #' @param tol tolerance for convergence (|logLik_k - logLik_{k-1}| < tol)
  #' @export
  #' @examples
  #' nx <- 60
  #' ny <- 60
  #' alpha <- 1
  #' lambda <- 4
  #' x <- rgamma(n = nx, shape = alpha, rate = lambda)
  #' y <- rgamma(n = ny, shape = alpha, rate = lambda)
  #' fit <- gammaFit(z = c(x,y))
  #' fit$MLE
  #' fit$MoM

  N <- length(z)
  muHat <- mean(z)
  sigmaHat <- sd(z)

  # MoM estimates
  alpha0 <- muHat^2 / sigmaHat^2
  lambda0 <- muHat / sigmaHat^2

  # Newton-Raphson routine
  sumLog <- sum(log(z))
  alphaHat <- rep(NA, maxIter)
  ll <- rep(NA, maxIter)
  alphaHat[1] <- alpha0
  ll[1] <- -Inf
  conv <- FALSE
  k <- 2

  while(!conv & k <= maxIter) {
    alphaHat[k] <- alphaHat[k-1] - g(N, muHat, sumLog, alphaHat[k-1]) / h(N, muHat, alphaHat[k-1])
    lambdaHat <- alphaHat[k] / muHat
    ll[k] <- N * alphaHat[k] * log(lambdaHat) - N * lgamma(alphaHat[k]) + 
             (alphaHat[k] - 1) * sumLog - lambdaHat * N * muHat
    conv <- (abs(ll[k] - ll[k-1])) < tol
    k <- k + 1
  }

  MLE <- c(alphaHat[k - 1], lambdaHat)
  MoM <- c(alpha0, lambda0)
  names(MLE) <- names(MoM) <- c("alpha", "lambda")

  return(list(MLE = MLE, MoM = MoM, ll = ll[2:(k-1)], alphaHat = alphaHat[2:(k-1)], 
              k = k-1, conv = conv))
}