# gammaDist

This is an R package that implements gamma difference distributions under the null of equal distributions, including the cumulative distribution function (CDF) derived by Klar (2015), and saddlepoint approximations to the CDF and density. In particular, if $X_i \sim \text{Gamma}(\alpha, \lambda), i=1, \ldots, n_x$ and $Y_j \sim \text{Gamma}(\alpha, \lambda), j=1, \ldots, n_y$, gammDist computes the CDF and density for the random variable $Z = (1/n_x) \sum_i X_i - (1/n_y) \sum_j Y_j$. gammaDist also includes a function for computing the maximum likelihood and method of moments estimates for iid gamma random variables.

# Installation

```{r}
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

# get parameter estimates for pooled data
fit <- gammaFit(z = c(x,y))
fit$MLE
fit$MoM

# Plot true distribution and saddlepoint approximation
q <- seq(-0.2, 0.2, 0.001)
plot(q, pgammaDif(q = q, nx, ny, alpha, lambda),
     type = "l", lty = 1, lwd = 2, ylab = "F(q)")
lines(q, pgammaDifSaddle(q = q, nx, ny, alpha, lambda),
      col = "red", lty = 2, lwd = 2)
legend("topleft", legend = c("True", "Saddlepoint"), 
       col = c("black", "red"), lty = 1:2, inset = 0.02)

# Plot saddlepoint approximation of the density
plot(q, dgammaDifSaddle(z = q, nx, ny, alpha, lambda),
     type = "l", ylab = "f(z)", xlab = "z")
```

## References
Klar, B. (2015). A note on gamma difference distributions. Journal of Statistical Computation and Simulation 85, 3708–3715.