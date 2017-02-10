# saddle point approximation functions ----------------------------------------

K <- function(t, nx, ny, alpha, lambda) {
  #' cumulant generating function (CGF)
  #'
  #' @param t CGF argument
  #' @param nx sample size of first group
  #' @param ny sample size of second group
  #' @param alpha common parameter value
  #' @param lambda common parameter value
  #' @export

  -alpha*(nx*log(1-t/(nx*lambda)) + ny*log(1+t/(ny*lambda)))
}

Kprime <- function(t, nx, ny, alpha, lambda) {
  #' first derivative of cumulant generating function (CGF)
  #'  
  #' @param t CGF argument
  #' @param nx sample size of first group
  #' @param ny sample size of second group
  #' @param alpha common parameter value
  #' @param lambda common parameter value
  #' @export

  alpha * (nx + ny) * t / ((nx * lambda - t) * (ny * lambda + t))
}

Kprime2 <- function(t, nx, ny, alpha, lambda) {
  #' second derivative of cumulant generating function (CGF)
  #'
  #' @param t CGF argument
  #' @param nx sample size of first group
  #' @param ny sample size of second group
  #' @param alpha common parameter value
  #' @param lambda common parameter value
  #' @export

  a <- (nx * lambda - t) * (ny * lambda + t)
  alpha * (nx + ny) *(t^2 + nx*ny*lambda^2) / a^2
}

saddleEqn <- function(t, q, nx, ny, alpha, lambda){
  #' saddlepoint equation: solve to obtain the saddlepoint \eqn{\hat{t}}
  #' 
  #' @param t CGF argument
  #' @param q quantile
  #' @param nx sample size of first group
  #' @param ny sample size of second group
  #' @param alpha common parameter value
  #' @param lambda common parameter value
  #' @export

  Kprime(t, nx, ny, alpha, lambda) - q
}

dgammaDifSaddle <- function(z, nx, ny, alpha, lambda) {
  #' Saddlepoint density for difference in means, \eqn{\bar{x} - \bar{y}}
  #'
  #' This function calculates the saddlepoint density for difference in means,
  #' \eqn{\bar{x} - \bar{y}} of gamma random variables under the null of equal
  #'  distributions.
  #' @param z scalar or vector of quantiles
  #' @param nx sample size of first group
  #' @param ny sample size of second group
  #' @param alpha common parameter value
  #' @param lambda common parameter value
  #' @examples
  #' nx <- 60
  #' ny <- 60
  #' alpha <- 1
  #' lambda <- 4
  #' z <- seq(-0.2, 0.2, 0.001)
  #' plot(z, dgammaDifSaddle(z = z, nx, ny, alpha, lambda),
  #'      type = "l", ylab = "f(q)")
  #' @export

  tHat <- sapply(z, function(z) {uniroot(f = saddleEqn,
                  interval = c(-ny*lambda, nx*lambda),
                  q = z, nx, ny, alpha, lambda)$root})

  1 / sqrt(2 * pi * Kprime2(tHat, nx, ny, alpha, lambda)) * 
           exp(K(tHat, nx, ny, alpha, lambda) - tHat * z)
}

pgammaDifSaddle <- function(q, nx, ny, alpha, lambda, lower.tail = TRUE) {
  #' Saddlepoint distribution for difference in means, \eqn{\bar{x} - \bar{y}}
  #'
  #' This function calculates the saddlepoint cumulative distribution for 
  #' the difference in means, \eqn{\bar{x} - \bar{y}} of gamma
  #' random variables under the null of equal distributions. Only valid
  #' for \eqn{\bar{x} - \bar{y} \ne 0}.
  #' @param q scalar or vector of quantiles
  #' @param nx sample size of first group
  #' @param ny sample size of second group
  #' @param alpha common parameter value
  #' @param lambda common parameter value
  #' @examples
  #' nx <- 60
  #' ny <- 60
  #' alpha <- 1
  #' lambda <- 4
  #' q <- seq(-0.2, 0.2, 0.001)
  #' plot(q, pgammaDifSaddle(q = q, nx, ny, alpha, lambda),
  #'      type = "l", ylab = "F(q)")
  #' @export
  
  tHat <- sapply(q, function(q) {uniroot(f = saddleEqn,
                  interval = c(-ny*lambda, nx*lambda),
                  q = q, nx, ny, alpha, lambda)$root})
  wHat <- sign(tHat) * sqrt(2 * (tHat*q - K(tHat, nx, ny, alpha, lambda)))
  uHat <- tHat * sqrt(Kprime2(tHat, nx, ny, alpha, lambda))

  if (lower.tail) {
    pnorm(wHat, lower.tail = TRUE) + dnorm(wHat)*(1/wHat - 1/uHat)
  } else {
    pnorm(wHat, lower.tail = FALSE) - dnorm(wHat)*(1/wHat - 1/uHat)
  }
}