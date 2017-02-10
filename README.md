# gammaDist

This package implements distribution functions for gamma random variables, including the cumulative distribution function (CDF) for the difference in the means of gamma random variables, a saddlepoint approximation to the CDF, and a Newton-Raphson routine for computing the maximum likelihood estimates for iid gamma random variables. The CDF is derived by Klar (2015).

# Installation

```{r}
# install.packages("devtools")
# library("devtools")
install_github("bdsegal/gammaDist")
```

## Examples

```{r}
library(gammaDist)

# generate data
nx <- 60
ny <- 60
alpha <- 1
lambda <- 4
x <- rgamma(n = nx, shape = alpha, rate = lambda)
y <- rgamma(n = ny, shape = alpha, rate = lambda)

# MLE
fit <- gammaFit(z = c(x,y))
fit$MLE
fit$MoM

# Plot true distribution and saddlepoint approximation
q <- seq(-0.2, 0.2, 0.001)
plot(q, pgammaDif(q = q, nx, ny, alpha, lambda), type = "l", lty = 1, lwd = 2,
     ylab = "F(q)")
lines(q, pgammaDifSaddle(q = q, nx, ny, alpha, lambda), col = "red", 
      lty = 2, lwd = 2)
legend("topleft", legend = c("True", "Saddlepoint"), 
       col = c("black", "red"), lty = 1:2, inset = 0.02)

# Plot saddlepoint approximation of the density
plot(q, dgammaDifSaddle(z = q, nx, ny, alpha, lambda), type = "l",
     ylab = "F(z)", xlab = "z")
```

## References
Klar, B. (2015). A note on gamma difference distributions. Journal of Statistical Computation and Simulation 85, 3708–3715.