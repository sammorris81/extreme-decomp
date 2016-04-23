library(numDeriv)
set.seed(2000)
beta <- abs(rnorm(5))

ns <- 25
nt <- 4
y <- abs(matrix(rnorm(ns * nt), ns, nt)) + 10

np <- length(beta)
X.mu <- abs(array(rnorm(ns * nt * np), dim = c(ns, nt, np)))
X.mu[, , 1] <- 1

logsig <- matrix(rnorm(ns * nt), ns, nt)
theta <- abs(matrix(rnorm(ns * nt), ns, nt))

alpha <- 0.4
xi <- 0.2

grad_loglike_betamu(beta1 = beta, X.mu = X.mu, y = y, theta = theta,
                    logsig = logsig, xi = xi, thresh = -Inf, alpha = alpha)

grad(func = loglike_mu, x = beta, X.mu = X.mu, y = y, theta = theta,
     logsig = logsig, xi = xi, thresh = -Inf, alpha = alpha)

grad_logpost_betamu(beta1 = beta, beta.mu = 0, beta.sd = 2, X.mu = X.mu,
                    y = y, theta = theta, logsig = logsig, xi = xi,
                    thresh = -Inf, alpha = alpha)

grad(func = logpost_mu, x = beta, beta.mu = 0, beta.sd = 2, X.mu = X.mu,
     y = y, theta = theta, logsig = logsig, xi = xi,
     thresh = -Inf, alpha = alpha)

hess_logpost_betamu(beta1 = beta, beta.mu = 0, beta.sd = 2, X.mu = X.mu,
                    y = y, theta = theta, logsig = logsig, xi = xi,
                    thresh = -Inf, alpha = alpha)

hessian(func = logpost_mu, x = beta, beta.mu = 0, beta.sd = 2, X.mu = X.mu,
        y = y, theta = theta, logsig = logsig, xi = xi,
        thresh = -Inf, alpha = alpha)

library(microbenchmark)
microbenchmark(grad_loglike_betamu(beta1 = beta, X.mu = X.mu, y = y, theta = theta,
                                   logsig = logsig, xi = xi, thresh = -Inf, alpha = alpha),
               grad_logpost_betamu(beta1 = beta, beta.mu = 0, beta.sd = 0, X.mu = X.mu,
                                   y = y, theta = theta, logsig = logsig, xi = xi,
                                   thresh = -Inf, alpha = alpha))

X.sig <- abs(array(rnorm(ns * nt * np), dim = c(ns, nt, np)))
X.sig[, , 1] <- 1

mu <- matrix(rnorm(ns * nt), ns, nt)

grad_loglike_betasig(beta2 = beta, X.sig = X.sig, y = y, theta = theta,
                     mu = mu, xi = xi, thresh = -Inf, alpha = alpha)

grad(func = loglike_sig, x = beta, X.sig = X.sig, y = y, theta = theta,
     mu = mu, xi = xi, thresh = -Inf, alpha = alpha)

grad_logpost_betasig(beta2 = beta, beta.mu = 0, beta.sd = 2, X.sig = X.sig,
                    y = y, theta = theta, mu = mu, xi = xi,
                    thresh = -Inf, alpha = alpha)

grad(func = logpost_sig, x = beta, beta.mu = 0, beta.sd = 2, X.sig = X.sig,
     y = y, theta = theta, mu = mu, xi = xi,
     thresh = -Inf, alpha = alpha)


hess_logpost_betasig(beta2 = beta, beta.mu = 0, beta.sd = 2, X.sig = X.sig,
                     y = y, theta = theta, mu = mu, xi = xi,
                     thresh = -Inf, alpha = alpha)

hessian(func = logpost_sig, x = beta, beta.mu = 0, beta.sd = 2, X.sig = X.sig,
        y = y, theta = theta, mu = mu, xi = xi, thresh = -Inf, alpha = alpha)


hess_logpost_betasig(beta2 = beta, beta.mu = 5, beta.sd = 2, X.sig = X.sig,
                     y = y, theta = theta, mu = mu, xi = xi,
                     thresh = -Inf, alpha = alpha) /
hessian(func = logpost_sig, x = beta, beta.mu = 5, beta.sd = 2, X.sig = X.sig,
        y = y, theta = theta, mu = mu, xi = xi, thresh = -Inf, alpha = alpha)

microbenchmark(
grad_loglike_betasig(beta2 = beta, X.sig = X.sig, y = y, theta = theta,
                     mu = mu, xi = xi, thresh = -Inf, alpha = alpha),
grad_logpost_betasig(beta2 = beta, beta.mu = 0, beta.sd = 1, X.sig = X.sig,
                     y = y, theta = theta, mu = mu, xi = xi,
                     thresh = -Inf, alpha = alpha)
)



temp <- function(beta2, beta.mu, beta.sd, X.sig, y, theta, mu,
                 xi, thresh, alpha) {

  logsig <- 0
  p.sig <- dim(X.sig)[3]
  for (j in 1:p.sig) {
    logsig <- logsig + X.sig[, , j] * beta2[j]
  }

  sigma <- exp(logsig)
  res <- y - mu

  tx <- (1 + xi * res / sigma)

  return(sum(tx^(-1 / (alpha * xi) - 1) / sigma))
}

grad_temp <- function(beta2, beta.mu, beta.sd, X.sig, y, theta, mu,
                      xi, thresh, alpha) {

  logsig <- 0
  p.sig <- dim(X.sig)[3]
  for (j in 1:p.sig) {
    logsig <- logsig + X.sig[, , j] * beta2[j]
  }

  sigma <- exp(logsig)
  res <- y - mu

  tx <- (1 + xi * res / sigma)

  grad_part <- - tx^(-1 / (alpha * xi) - 1) * sigma^(-2) * sigma +
    sigma^(-1) * tx^(-1 / (alpha * xi) - 2) * (1 / (alpha * xi) + 1) *
    (xi * res / sigma)

  grad_part <- tx^(-1 / (alpha * xi) - 1) * (-1 + (alpha * xi + 1) * res * xi / (alpha * xi * tx * sigma)) / sigma

  grad <- rep(0, p.sig)
  for (j in 1:p.sig) {
    grad[j] <- sum(grad_part * X.sig[, , j])
  }

  return(grad)
}

grad(func = temp, x = beta, beta.mu = 0, beta.sd = 2, X.sig = X.sig,
     y = y, theta = theta, mu = mu, xi = xi,
     thresh = -Inf, alpha = alpha)

grad_temp(beta2 = beta, beta.mu = 0, beta.sd = 2, X.sig = X.sig,
          y = y, theta = theta, mu = mu, xi = xi,
          thresh = -Inf, alpha = alpha)