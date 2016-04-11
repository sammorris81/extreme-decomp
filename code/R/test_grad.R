set.seed(1234)
beta1 <- beta2 <- runif(3)
X.mu <- array(1, dim = c(5, 2, 3))
alpha <- 0.4
X.mu[, , 1] <- 1
X.mu[, , 2] <- matrix(c(rep(-0.5, 5), rep(0.5, 5)), 5, 2)
X.mu[, , 3] <- rnorm(10)

X.sig <- X.mu

Y <- abs(matrix(rnorm(10), 5, 2))
theta <- abs(matrix(rnorm(10), 5, 2))
xi <- 0.2
logsig <- matrix(rnorm(10), 5, 2)

loglike_mu(beta1 = beta1, X.mu = X.mu, y = Y, theta = theta, logsig = logsig,
           xi = xi, thresh = -Inf, alpha = alpha)

grad(func = loglike_mu, x = beta1, X.mu = X.mu, y = Y, theta = theta, logsig = logsig,
     xi = xi, thresh = -Inf, alpha = alpha)

mu <- 0
for (i in 1:3) {
  mu <- mu + X.mu[, , i] * beta1[i]
}

grad(func = loglike_mu, x = beta1, X.mu = X.mu, y = Y, theta = theta, logsig = logsig,
     xi = xi, thresh = -Inf, alpha = alpha)
grad_loglike_betamu(beta1 = beta1, X.mu = X.mu, y = Y, theta = theta,
                    logsig = logsig, xi = xi, thresh = -Inf, alpha = alpha)

library(microbenchmark)
microbenchmark(grad(func = loglike_mu, x = beta1, X.mu = X.mu, y = Y, theta = theta, logsig = logsig,
                    xi = xi, thresh = -Inf, alpha = alpha),
               grad_loglike_betamu(beta1 = beta1, X.mu = X.mu, y = Y, theta = theta,
                                   logsig = logsig, xi = xi, thresh = -Inf, alpha = alpha))

sigma <- exp(logsig)
mu_star  <- mu + sigma * ((theta^xi) - 1) / xi
sig_star <- alpha * sigma * (theta^xi)
xi_star  <- alpha * xi

t1 <- (1 + xi_star * (Y - mu_star) / sig_star)
t2 <- (1 + alpha * xi * (Y - (mu + sigma * (theta^xi - 1) / xi))/(alpha * sigma * theta^xi))
t3 <- 1 + xi * (Y - mu - sigma * theta^xi / xi + sigma / xi)/(sigma * theta^xi)
t4 <- xi * (Y - mu)/(sigma * theta^xi) + 1 / theta^xi

loglike_sig(beta2 = beta2, X.sig = X.sig, y = Y, theta = theta, mu = mu,
            xi = xi, thresh = Inf, alpha = alpha)
logsig <- 0
for (i in 1:3) {
  logsig <- logsig + X.sig[, , i] * beta2[i]
}

sum(loglike(y = Y, theta = theta, mu = mu, logsig = logsig, xi = xi,
            thresh = Inf, alpha = alpha))

grad(func = loglike_sig, x = beta2, X.sig = X.sig, y = Y, theta = theta,
     mu = mu, xi = xi, thresh = -Inf, alpha = alpha)
grad_loglike_betasig(beta2 = beta2, X.sig = X.sig, y = Y, theta = theta,
                     mu = mu, xi = xi, thresh = -Inf, alpha = alpha)

microbenchmark(grad(func = loglike_sig, x = beta2, X.sig = X.sig, y = Y, theta = theta,
                    mu = mu, xi = xi, thresh = -Inf, alpha = alpha),
               grad_loglike_betasig(beta2 = beta2, X.sig = X.sig, y = Y, theta = theta,
                                    mu = mu, xi = xi, thresh = -Inf, alpha = alpha))
