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

grad_logpost_betamu(beta1 = beta, beta.mu = 0, beta.sd = 1, X.mu = X.mu,
                    y = y, theta = theta, logsig = logsig, xi = xi,
                    thresh = -Inf, alpha = alpha)

grad(func = logpost_mu, x = beta, beta.mu = 0, beta.sd = 1, X.mu = X.mu,
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

grad_logpost_betasig(beta2 = beta, beta.mu = 0, beta.sd = 1, X.sig = X.sig,
                    y = y, theta = theta, mu = mu, xi = xi,
                    thresh = -Inf, alpha = alpha)

microbenchmark(
grad_loglike_betasig(beta2 = beta, X.sig = X.sig, y = y, theta = theta,
                     mu = mu, xi = xi, thresh = -Inf, alpha = alpha),
grad_logpost_betasig(beta2 = beta, beta.mu = 0, beta.sd = 1, X.sig = X.sig,
                     y = y, theta = theta, mu = mu, xi = xi,
                     thresh = -Inf, alpha = alpha)
)

grad(func = logpost_sig, x = beta, beta.mu = 0, beta.sd = 1, X.sig = X.sig,
     y = y, theta = theta, mu = mu, xi = xi,
     thresh = -Inf, alpha = alpha)

