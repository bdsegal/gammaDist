# True CDF functions ----------------------------------------------------------

cdfSub <- function(v, q, alpha1, alpha2, lambda1, lambda2) {
  #' Utility function
  #'
  #' @param v integration term
  #' @param q quantile
  #' @param alpha1 parameter value
  #' @param alpha2 parameter value
  #' @param lambda1 parameter value
  #' @param lambda2 parameter value
  #' @export

  v^(alpha2 - 1) * exp(-lambda2 * v) * Igamma(alpha1, lambda1*(v + q))
}

pgammaDif <- function(q, nx, ny, alpha, lambda) {
  #' True CDF for difference in means \eqn{\bar{x} - \bar{y}}
  #'
  #' This function computes the true CDF for the difference in means \eqn{\bar{x} - \bar{y}}
  #' of gamma random variables under the null of equality in distributions.
  #' CDF derived by Klar (2015).
  #' @param q scalar or vector of quantiles
  #' @param nx sample size of first group
  #' @param ny samle size of second group
  #' @param alpha common parameter value
  #' @param lambda common parameter value
  #' @examples
  #' nx <- 60
  #' ny <- 60
  #' alpha <- 1
  #' lambda <- 4
  #' q <- seq(-0.2, 0.2, 0.001)
  #' plot(q, pgammaDif(q = q, nx, ny, alpha, lambda),
  #'      type = "l", ylab = "F(q)")
  #' @export

  alpha1 <- nx * alpha
  alpha2 <- ny * alpha
  lambda1 <- nx * lambda
  lambda2 <- ny * lambda

  # integrate with adaptive quadrature
  A <- integrate(cdfSub, lower = max(0, -x), upper = Inf, 
                 q = q, alpha1 = alpha1, alpha2 = alpha2,
                 lambda1 = lambda1, lambda2 = lambda2)$value

  exp(alpha2*log(lambda2) - lgamma(alpha1) - lgamma(alpha2) + log(A))
}

pgammaDif <- Vectorize(FUN = pgammaDif, vectorize.args = "q")