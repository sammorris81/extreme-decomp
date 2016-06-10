rm(list=ls())
library(fields)
library(Rcpp)
library(emulator)
library(microbenchmark)
library(SpatialExtremes)
library(numDeriv)
library(fields)

#### Verify gradients - on beta ####
rm(list=ls())
source("../../../usefulR/usefulfunctions.R", chdir = TRUE)
openblas.set.num.threads(3)
source("auxfunctions.R", chdir = TRUE)
source("updatemodel.R", chdir = TRUE)

set.seed(2000)
ns <- 10
nt <- 3

phi.t <- 0.2
s <- cbind(runif(ns), runif(ns))
d <- rdist(s)
Sigma.t <- exp(-d / phi.t)
Qb.t    <- chol2inv(chol(Sigma.t))

tau1.int.t  <- rgamma(1, 1, 1)
tau1.time.t <- rgamma(1, 1, 1)
tau2.int.t  <- rgamma(1, 1, 1)
tau2.time.t <- rgamma(1, 1, 1)
beta1.int.t <- beta1.time.t <- matrix(0, ns, nt)
beta2.int.t <- beta2.time.t <- matrix(0, ns, nt)
for (t in 1:nt) {
  beta1.int.t[, t]  <- t(chol(Sigma.t)) %*% rnorm(ns) / sqrt(tau1.int.t[1])
  beta1.time.t[, t] <- t(chol(Sigma.t)) %*% rnorm(ns) / sqrt(tau1.time.t[1])
  beta2.int.t[, t]  <- t(chol(Sigma.t)) %*% rnorm(ns) / sqrt(tau2.int.t[1])
  beta2.time.t[, t] <- t(chol(Sigma.t)) %*% rnorm(ns) / sqrt(tau2.time.t[1])
}
beta1.int.mn.t <- beta1.time.mn.t <- 0
beta2.int.mn.t <- beta2.time.mn.t <- 0

mu.t <- matrix(0, ns, nt)
ls.t <- matrix(0, ns, nt)
time <- (1:nt - nt / 2) / nt
for (t in 1:nt) {
  mu.t[, t] <- beta1.int.t[, t] + beta1.time.t[, t] * time[t]
  ls.t[, t] <- beta2.int.t[, t] + beta2.time.t[, t] * time[t]
}
xi.t <- 0
y.t <- rgev(n = ns * nt, loc = mu.t, exp(ls.t), xi.t)

nknots <- 4
theta.t <- matrix(abs(rnorm(ns * nt)), ns, nt)
alpha.t <- 0.4
thresh.t <- matrix(median(y.t), ns, nt)

xi.t <- 0
lp.beta1.int <- logpost.beta1.int(beta.int = beta1.int.t[, t],
                                  beta.int.mn = beta1.int.mn.t,
                                  tau = tau1.int.t, Qb = Qb.t,
                                  beta.time = beta1.time.t[, t], time = time[t],
                                  y = y.t[, t], ls = ls.t[, t], xi = xi.t,
                                  theta = theta.t[, t], thresh = thresh.t[, t],
                                  alpha = alpha.t)

mean(grad(func = logpost.beta1.int, x = beta1.int.t[, t],
          beta.int.mn = beta1.int.mn.t, tau = tau1.int.t, Qb = Qb.t,
          beta.time = beta1.time.t[, t], time = time[t],
          y = y.t[, t], ls = ls.t[, t], xi = xi.t,
          theta = theta.t[, t], thresh = thresh.t[, t],
          alpha = alpha.t) /
       logpost.beta1.int.grad(beta.int = beta1.int.t[, t],
                              beta.int.mn = beta1.int.mn.t,
                              tau = tau1.int.t, Qb = Qb.t,
                              beta.time = beta1.time.t[, t], time = time[t],
                              y = y.t[, t], ls = ls.t[, t], xi = xi.t,
                              theta = theta.t[, t], thresh = thresh.t[, t],
                              alpha = alpha.t))

sd(grad(func = logpost.beta1.int, x = beta1.int.t[, t],
        beta.int.mn = beta1.int.mn.t, tau = tau1.int.t, Qb = Qb.t,
        beta.time = beta1.time.t[, t], time = time[t],
        y = y.t[, t], ls = ls.t[, t], xi = xi.t,
        theta = theta.t[, t], thresh = thresh.t[, t],
        alpha = alpha.t) /
     logpost.beta1.int.grad(beta.int = beta1.int.t[, t],
                            beta.int.mn = beta1.int.mn.t,
                            tau = tau1.int.t, Qb = Qb.t,
                            beta.time = beta1.time.t[, t], time = time[t],
                            y = y.t[, t], ls = ls.t[, t], xi = xi.t,
                            theta = theta.t[, t], thresh = thresh.t[, t],
                            alpha = alpha.t))

lp.beta1.time <- logpost.beta1.time(beta.time = beta1.time.t[, t],
                                    beta.time.mn = beta1.time.mn.t,
                                    time = time[t], tau = tau1.int.t, Qb = Qb.t,
                                    beta.int = beta1.int.t[, t],
                                    y = y.t[, t], ls = ls.t[, t], xi = xi.t,
                                    theta = theta.t[, t],
                                    thresh = thresh.t[, t], alpha = alpha.t)

mean(grad(func = logpost.beta1.time, x = beta1.time.t[, t],
          beta.time.mn = beta1.time.mn.t, time = time[t], tau = tau1.int.t,
          Qb = Qb.t, beta.int = beta1.int.t[, t],
          y = y.t[, t], ls = ls.t[, t], xi = xi.t,
          theta = theta.t[, t], thresh = thresh.t[, t],
          alpha = alpha.t) /
       logpost.beta1.time.grad(beta.time = beta1.time.t[, t],
                               beta.time.mn = beta1.time.mn.t, time = time[t],
                               tau = tau1.int.t, Qb = Qb.t,
                               beta.int = beta1.int.t[, t],
                               y = y.t[, t], ls = ls.t[, t], xi = xi.t,
                               theta = theta.t[, t],
                               thresh = thresh.t[, t], alpha = alpha.t))

sd(grad(func = logpost.beta1.time, x = beta1.time.t[, t],
        beta.time.mn = beta1.time.mn.t, time = time[t], tau = tau1.int.t,
        Qb = Qb.t, beta.int = beta1.int.t[, t],
        y = y.t[, t], ls = ls.t[, t], xi = xi.t,
        theta = theta.t[, t], thresh = thresh.t[, t],
        alpha = alpha.t) /
     logpost.beta1.time.grad(beta.time = beta1.time.t[, t],
                             beta.time.mn = beta1.time.mn.t, time = time[t],
                             tau = tau1.int.t, Qb = Qb.t,
                             beta.int = beta1.int.t[, t],
                             y = y.t[, t], ls = ls.t[, t], xi = xi.t,
                             theta = theta.t[, t],
                             thresh = thresh.t[, t], alpha = alpha.t))

this.grad <- logpost.beta1.grad(beta.int = beta1.int.t[, t],
                                beta.int.mn = beta1.int.mn.t,
                                beta.time = beta1.time.t[, t],
                                beta.time.mn = beta1.time.mn.t,
                                time = time[t],
                                tau = tau1.int.t, Qb = Qb.t,
                                y = y.t[, t], ls = ls.t[, t], xi = xi.t,
                                theta = theta.t[, t], thresh = thresh.t[, t],
                                alpha = alpha.t)

logpost.beta1.int.grad(beta.int = beta1.int.t[, t],
                       beta.int.mn = beta1.int.mn.t,
                       tau = tau1.int.t, Qb = Qb.t,
                       beta.time = beta1.time.t[, t], time = time[t],
                       y = y.t[, t], ls = ls.t[, t], xi = xi.t,
                       theta = theta.t[, t], thresh = thresh.t[, t],
                       alpha = alpha.t)

logpost.beta1.time.grad(beta.time = beta1.time.t[, t],
                        beta.time.mn = beta1.time.mn.t, time = time[t],
                        tau = tau1.int.t, Qb = Qb.t,
                        beta.int = beta1.int.t[, t],
                        y = y.t[, t], ls = ls.t[, t], xi = xi.t,
                        theta = theta.t[, t],
                        thresh = thresh.t[, t], alpha = alpha.t)

microbenchmark(logpost.beta1.grad(beta.int = beta1.int.t[, t],
                                  beta.int.mn = beta1.int.mn.t,
                                  beta.time = beta1.time.t[, t],
                                  beta.time.mn = beta1.time.mn.t,
                                  time = time[t],
                                  tau = tau1.int.t, Qb = Qb.t,
                                  y = y.t[, t], ls = ls.t[, t], xi = xi.t,
                                  theta = theta.t[, t], thresh = thresh.t[, t],
                                  alpha = alpha.t),
               logpost.beta1.int.grad(beta.int = beta1.int.t[, t],
                                      beta.int.mn = beta1.int.mn.t,
                                      tau = tau1.int.t, Qb = Qb.t,
                                      beta.time = beta1.time.t[, t], time = time[t],
                                      y = y.t[, t], ls = ls.t[, t], xi = xi.t,
                                      theta = theta.t[, t], thresh = thresh.t[, t],
                                      alpha = alpha.t),
               logpost.beta1.time.grad(beta.time = beta1.time.t[, t],
                                       beta.time.mn = beta1.time.mn.t, time = time[t],
                                       tau = tau1.int.t, Qb = Qb.t,
                                       beta.int = beta1.int.t[, t],
                                       y = y.t[, t], ls = ls.t[, t], xi = xi.t,
                                       theta = theta.t[, t],
                                       thresh = thresh.t[, t], alpha = alpha.t))

lp.beta2.int <- logpost.beta2.int(beta.int = beta2.int.t[, t],
                                  beta.int.mn = beta2.int.mn.t,
                                  tau = tau2.int.t, Qb = Qb.t,
                                  beta.time = beta1.time.t[, t], time = time[t],
                                  y = y.t[, t], mu = mu.t[, t], xi = xi.t,
                                  theta = theta.t[, t], thresh = thresh.t[, t],
                                  alpha = alpha.t)

mean(grad(func = logpost.beta2.int, x = beta2.int.t[, t],
          beta.int.mn = beta2.int.mn.t, tau = tau2.int.t, Qb = Qb.t,
          beta.time = beta2.time.t[, t], time = time[t],
          y = y.t[, t], mu = mu.t[, t], xi = xi.t,
          theta = theta.t[, t], thresh = thresh.t[, t],
          alpha = alpha.t) /
       logpost.beta2.int.grad(beta.int = beta2.int.t[, t],
                              beta.int.mn = beta2.int.mn.t,
                              tau = tau2.int.t, Qb = Qb.t,
                              beta.time = beta2.time.t[, t], time = time[t],
                              y = y.t[, t], mu = mu.t[, t], xi = xi.t,
                              theta = theta.t[, t], thresh = thresh.t[, t],
                              alpha = alpha.t))

sd(grad(func = logpost.beta2.int, x = beta2.int.t[, t],
        beta.int.mn = beta2.int.mn.t, tau = tau2.int.t, Qb = Qb.t,
        beta.time = beta2.time.t[, t], time = time[t],
        y = y.t[, t], mu = mu.t[, t], xi = xi.t,
        theta = theta.t[, t], thresh = thresh.t[, t],
        alpha = alpha.t) /
     logpost.beta2.int.grad(beta.int = beta2.int.t[, t],
                            beta.int.mn = beta2.int.mn.t,
                            tau = tau2.int.t, Qb = Qb.t,
                            beta.time = beta2.time.t[, t], time = time[t],
                            y = y.t[, t], mu = mu.t[, t], xi = xi.t,
                            theta = theta.t[, t], thresh = thresh.t[, t],
                            alpha = alpha.t))

lp.beta2.time <- logpost.beta2.time(beta.time = beta2.time.t[, t],
                                    beta.time.mn = beta2.time.mn.t,
                                    time = time[t], tau = tau2.int.t, Qb = Qb.t,
                                    beta.int = beta1.int.t[, t],
                                    y = y.t[, t], mu = mu.t[, t], xi = xi.t,
                                    theta = theta.t[, t], thresh = thresh.t[, t],
                                    alpha = alpha.t)

mean(grad(func = logpost.beta2.time, x = beta2.time.t[, t],
          beta.time.mn = beta2.time.mn.t, time = time[t], tau = tau2.int.t,
          Qb = Qb.t, beta.int = beta2.int.t[, t],
          y = y.t[, t], mu = mu.t[, t], xi = xi.t,
          theta = theta.t[, t], thresh = thresh.t[, t],
          alpha = alpha.t) /
       logpost.beta2.time.grad(beta.time = beta2.time.t[, t],
                               beta.time.mn = beta2.time.mn.t,
                               time = time[t], tau = tau2.int.t, Qb = Qb.t,
                               beta.int = beta2.int.t[, t],
                               y = y.t[, t], mu = mu.t[, t], xi = xi.t,
                               theta = theta.t[, t], thresh = thresh.t[, t],
                               alpha = alpha.t))

sd(grad(func = logpost.beta2.time, x = beta2.time.t[, t],
        beta.time.mn = beta2.time.mn.t, time = time[t], tau = tau2.int.t,
        Qb = Qb.t, beta.int = beta2.int.t[, t],
        y = y.t[, t], mu = mu.t[, t], xi = xi.t,
        theta = theta.t[, t], thresh = thresh.t[, t],
        alpha = alpha.t) /
     logpost.beta2.time.grad(beta.time = beta2.time.t[, t],
                             beta.time.mn = beta2.time.mn.t,
                             time = time[t], tau = tau2.int.t, Qb = Qb.t,
                             beta.int = beta2.int.t[, t],
                             y = y.t[, t], mu = mu.t[, t], xi = xi.t,
                             theta = theta.t[, t], thresh = thresh.t[, t],
                             alpha = alpha.t))

#### testing means ####
rm(list=ls())
source("../../../usefulR/usefulfunctions.R", chdir = TRUE)
openblas.set.num.threads(3)
source("auxfunctions.R", chdir = TRUE)
source("updatemodel.R", chdir = TRUE)

set.seed(2000)
ns <- 600
nt <- 30

phi.t <- 0.2
s <- cbind(runif(ns), runif(ns))
d <- rdist(s)
Sigma.t <- exp(-d / phi.t)
Qb.t    <- chol2inv(chol(Sigma.t))

tau1.int.t  <- rgamma(1, 10, 10)
tau1.time.t <- rgamma(1, 10, 10)
beta1.int.t <- beta1.time.t <- matrix(0, ns, nt)
beta1.int.mn.t <- 100
beta1.time.mn.t <- 0
for (t in 1:nt) {
  beta1.int.t[, t] <- beta1.int.mn.t +
    t(chol(Sigma.t)) %*% rnorm(ns) / sqrt(tau1.int.t)
  beta1.time.t[, t] <- beta1.time.mn.t +
    t(chol(Sigma.t)) %*% rnorm(ns) / sqrt(tau1.time.t)
}
# beta1.int.mn.t <- beta1.time.mn.t <- 0

SS1.int.t <- SS1.time.t <- rep(0, nt)
mu.t <- matrix(0, ns, nt)
time <- (1:nt - nt / 2) / nt
for (t in 1:nt) {
  mu.t[, t] <- beta1.int.t[, t] + beta1.time.t[, t] * time[t]
  SS1.int.t[t]  <- quad.form(Qb.t, beta1.int.t[, t] - beta1.int.mn.t)
  SS1.time.t[t] <- quad.form(Qb.t, beta1.time.t[, t] - beta1.time.mn.t)
}

ls.t <- matrix(0, ns, nt)
xi.t <- 0
y.t <- rgev(n = ns * nt, loc = mu.t, exp(ls.t), xi.t)

thresh.t <- matrix(-Inf, ns, nt)
theta.t <- matrix(1, ns, nt)
theta.xi.t <- theta.t^xi.t
alpha.t <- 1

# initialize values
beta1.int.mn  <- 0
beta1.time.mn <- 0
beta1.pri.sd  <- 100

curll <- loglike(y = y.t, mu = mu, ls = ls.t, xi = xi.t,
                 theta = theta.t, theta.xi = theta.xi.t,
                 thresh = thresh.t, alpha = alpha.t)

niters <- 10000
burn   <- 8000
keep.beta1.int.mn <- keep.beta1.time.mn <- rep(0, niters)

set.seed(3366)  # demo
for (iter in 1:niters) {
  this.update <- updateGPMean(beta.sd = beta1.pri.sd, Qb = Qb.t,
                              beta.int = beta1.int.t, tau.int = tau1.int.t,
                              beta.time = beta1.time.t, tau.time = tau1.time.t)
  beta1.int.mn  <- this.update$beta.int.mn
  beta1.time.mn <- this.update$beta.time.mn

  keep.beta1.int.mn[iter] <- beta1.int.mn
  keep.beta1.time.mn[iter] <- beta1.time.mn

  if (iter %% 1000 == 0) {
    par(mfrow = c(1, 2))

    plot(keep.beta1.int.mn[1:iter], type = "l",
         main = paste("Beta 1 intercept mean: ", beta1.int.mn.t, sep = ""))
    plot(keep.beta1.time.mn[1:iter], type = "l",
         main = paste("Beta 1 time mean: ", beta1.time.mn.t, sep = ""))

  }
}

#### testing taus ####
rm(list=ls())
source("../../../usefulR/usefulfunctions.R", chdir = TRUE)
openblas.set.num.threads(3)
source("auxfunctions.R", chdir = TRUE)
source("updatemodel.R", chdir = TRUE)

set.seed(2000)
ns <- 600
nt <- 30

phi.t <- 0.2
s <- cbind(runif(ns), runif(ns))
d <- rdist(s)
Sigma.t <- exp(-d / phi.t)
Qb.t    <- chol2inv(chol(Sigma.t))

tau1.int.t  <- rgamma(1, 1, 1)
tau1.time.t <- rgamma(1, 1, 1)
beta1.int.t <- beta1.time.t <- matrix(0, ns, nt)
beta1.int.mn.t <- 5
beta1.time.mn.t <- 7
for (t in 1:nt) {
  beta1.int.t[, t] <- beta1.int.mn.t +
    t(chol(Sigma.t)) %*% rnorm(ns) / sqrt(tau1.int.t)
  beta1.time.t[, t] <- beta1.time.mn.t +
    t(chol(Sigma.t)) %*% rnorm(ns) / sqrt(tau1.time.t)
}

SS1.int.t <- SS1.time.t <- rep(0, nt)
mu.t <- matrix(0, ns, nt)
time <- (1:nt - nt / 2) / nt
for (t in 1:nt) {
  mu.t[, t] <- beta1.int.t[, t] + beta1.time.t[, t] * time[t]
  SS1.int.t[t]  <- quad.form(Qb.t, beta1.int.t[, t] - beta1.int.mn.t)
  SS1.time.t[t] <- quad.form(Qb.t, beta1.time.t[, t] - beta1.time.mn.t)
}

ls.t <- matrix(0, ns, nt)
xi.t <- 0
y.t <- rgev(n = ns * nt, loc = mu.t, exp(ls.t), xi.t)

thresh.t <- matrix(-Inf, ns, nt)
theta.t <- matrix(1, ns, nt)
theta.xi.t <- theta.t^xi.t
alpha.t <- 1

# initialize values
beta1.int.mn  <- 0
beta1.time.mn <- 0
beta1.pri.sd  <- 100

curll <- loglike(y = y.t, mu = mu, ls = ls.t, xi = xi.t,
                 theta = theta.t, theta.xi = theta.xi.t,
                 thresh = thresh.t, alpha = alpha.t)

niters <- 10000
burn   <- 8000
keep.tau1.int <- keep.tau1.time <- rep(0, niters)

set.seed(3366)  # demo
for (iter in 1:niters) {
  this.update <- updateGPTau(SS.int = SS1.int.t, SS.time = SS1.time.t,
                             tau.a = 0.1, tau.b = 0.1, ns = ns)
  tau1.int  <- this.update$tau.int
  tau1.time <- this.update$tau.time

  keep.tau1.int[iter]  <- tau1.int
  keep.tau1.time[iter] <- tau1.time

  if (iter %% 1000 == 0) {
    par(mfrow = c(1, 2))

    plot(keep.tau1.int[1:iter], type = "l",
         main = paste("Tau 1 intercept: ", tau1.int.t, sep = ""))
    plot(keep.tau1.time[1:iter], type = "l",
         main = paste("Tau 1 time: ", tau1.time.t, sep = ""))

  }
}

#### testing means and taus ####
rm(list=ls())
source("../../../usefulR/usefulfunctions.R", chdir = TRUE)
openblas.set.num.threads(3)
source("auxfunctions.R", chdir = TRUE)
source("updatemodel.R", chdir = TRUE)

set.seed(2000)
ns <- 600
nt <- 30

phi.t <- 0.2
s <- cbind(runif(ns), runif(ns))
d <- rdist(s)
Sigma.t <- exp(-d / phi.t)
Qb.t    <- chol2inv(chol(Sigma.t))

tau1.int.t  <- rgamma(1, 1, 1)
tau1.time.t <- rgamma(1, 1, 1)
beta1.int.t <- beta1.time.t <- matrix(0, ns, nt)
beta1.int.mn.t <- 5
beta1.time.mn.t <- 7
for (t in 1:nt) {
  beta1.int.t[, t] <- beta1.int.mn.t +
    t(chol(Sigma.t)) %*% rnorm(ns) / sqrt(tau1.int.t)
  beta1.time.t[, t] <- beta1.time.mn.t +
    t(chol(Sigma.t)) %*% rnorm(ns) / sqrt(tau1.time.t)
}

SS1.int.t <- SS1.time.t <- rep(0, nt)
mu.t <- matrix(0, ns, nt)
time <- (1:nt - nt / 2) / nt
for (t in 1:nt) {
  mu.t[, t] <- beta1.int.t[, t] + beta1.time.t[, t] * time[t]
  SS1.int.t[t]  <- quad.form(Qb.t, beta1.int.t[, t] - beta1.int.mn.t)
  SS1.time.t[t] <- quad.form(Qb.t, beta1.time.t[, t] - beta1.time.mn.t)
}

ls.t <- matrix(0, ns, nt)
xi.t <- 0
y.t <- rgev(n = ns * nt, loc = mu.t, exp(ls.t), xi.t)

thresh.t <- matrix(-Inf, ns, nt)
theta.t <- matrix(1, ns, nt)
theta.xi.t <- theta.t^xi.t
alpha.t <- 1

# initialize values
beta1.int.mn  <- 0
beta1.time.mn <- 0
beta1.pri.sd  <- 100
tau1.int      <- 0.1
tau1.time     <- 0.1

SS1.int <- SS1.time <- rep(0, nt)
for (t in 1:nt) {
  SS1.int[t]  <- quad.form(Qb.t, beta1.int.t[, t] - beta1.int.mn)
  SS1.time[t] <- quad.form(Qb.t, beta1.time.t[, t] - beta1.time.mn)
}

curll <- loglike(y = y.t, mu = mu, ls = ls.t, xi = xi.t,
                 theta = theta.t, theta.xi = theta.xi.t,
                 thresh = thresh.t, alpha = alpha.t)

niters <- 10000
burn   <- 8000
keep.beta1.int.mn <- keep.beta1.time.mn <- rep(0, niters)
keep.tau1.int <- keep.tau1.time <- rep(0, niters)

set.seed(3366)  # demo
for (iter in 1:niters) {
  this.update <- updateGPMean(beta.sd = beta1.pri.sd, Qb = Qb.t,
                              beta.int = beta1.int.t, tau.int = tau1.int,
                              beta.time = beta1.time.t, tau.time = tau1.time)
  beta1.int.mn  <- this.update$beta.int.mn
  SS1.int       <- this.update$SS.int
  beta1.time.mn <- this.update$beta.time.mn
  SS1.time      <- this.update$SS.time

  keep.beta1.int.mn[iter] <- beta1.int.mn
  keep.beta1.time.mn[iter] <- beta1.time.mn

  this.update <- updateGPTau(SS.int = SS1.int, SS.time = SS1.time,
                             tau.a = 0.1, tau.b = 0.1, ns = ns)
  tau1.int  <- this.update$tau.int
  tau1.time <- this.update$tau.time

  keep.tau1.int[iter]  <- tau1.int
  keep.tau1.time[iter] <- tau1.time

  if (iter %% 2000 == 0) {
    par(mfrow = c(2, 2))

    plot(keep.beta1.int.mn[1:iter], type = "l",
         main = paste("Beta 1 intercept mean: ", beta1.int.mn.t, sep = ""))
    plot(keep.beta1.time.mn[1:iter], type = "l",
         main = paste("Beta 1 time mean: ", beta1.time.mn.t, sep = ""))

    plot(keep.tau1.int[1:iter], type = "l",
         main = paste("Tau 1 intercept: ", tau1.int.t, sep = ""))
    plot(keep.tau1.time[1:iter], type = "l",
         main = paste("Tau 1 time: ", tau1.time.t, sep = ""))

  }
}

#### testing beta1 ####
rm(list=ls())
source("../../../usefulR/usefulfunctions.R", chdir = TRUE)
openblas.set.num.threads(3)
source("auxfunctions.R", chdir = TRUE)
source("updatemodel.R", chdir = TRUE)

set.seed(2000)
ns <- 600
nt <- 30

phi.t <- 0.2
s <- cbind(runif(ns), runif(ns))
d <- rdist(s)
Sigma.t <- exp(-d / phi.t)
Qb.t    <- chol2inv(chol(Sigma.t))

tau1.int.t  <- rgamma(1, 1, 1)
tau1.time.t <- rgamma(1, 1, 1)
beta1.int.t <- beta1.time.t <- matrix(0, ns, nt)
for (t in 1:nt) {
  beta1.int.t[, t] <- t(chol(Sigma.t)) %*% rnorm(ns) / sqrt(tau1.int.t[1])
  beta1.time.t[, t] <- t(chol(Sigma.t)) %*% rnorm(ns) / sqrt(tau1.time.t[1])
}
beta1.int.mn.t <- beta1.time.mn.t <- 0

mu.t <- matrix(0, ns, nt)
time <- (1:nt - nt / 2) / nt
for (t in 1:nt) {
  mu.t[, t] <- beta1.int.t[, t] + beta1.time.t [, t] * time[t]
}
ls.t <- matrix(0, ns, nt)
xi.t <- 0
y.t <- rgev(n = ns * nt, loc = mu.t, exp(ls.t), xi.t)

# initialize values
beta1.int  <- beta1.int.t + rnorm(ns * nt)
beta1.time <- beta1.time.t + rnorm(ns * nt)
mu <- matrix(0, ns, nt)
SS1.int <- SS1.time <- rep(0, nt)
for (t in 1:nt) {
  mu[, t] <- beta1.int[, t] + beta1.time[, t] * time[t]
  SS1.int[t]  <- quad.form(Qb.t, beta1.int[, t] - beta1.int.mn.t)
  SS1.time[t] <- quad.form(Qb.t, beta1.time[, t] - beta1.time.mn.t)
}
thresh.t <- matrix(-Inf, ns, nt)
theta.t <- matrix(1, ns, nt)
theta.xi.t <- theta.t^xi.t
alpha.t <- 1
curll <- loglike(y = y.t, mu = mu, ls = ls.t, xi = xi.t,
                 theta = theta.t, theta.xi = theta.xi.t,
                 thresh = thresh.t, alpha = alpha.t)

niters <- 10000
burn   <- 8000
keep.beta1.int <- keep.beta1.time <- array(0, dim = c(niters, ns, nt))
acc.beta1  <- att.beta1 <- MH.beta1 <- matrix(0.01, ns, nt)

set.seed(3366)  # demo
Rprof(filename = "Rprof.out", line.profiling = TRUE)
for (iter in 1:niters) {
  this.update <- updateBeta1(beta.int = beta1.int, beta.int.mn = beta1.int.mn.t,
                             SS.int = SS1.int, tau.int = tau1.int.t,
                             beta.time = beta1.time,
                             beta.time.mn = beta1.time.mn.t,
                             SS.time = SS1.time, tau.time = tau1.time.t,
                             mu = mu, time = time,
                             y = y.t, theta = theta.t, theta.xi = theta.xi.t,
                             ls = ls.t, xi = xi.t, thresh = thresh.t,
                             alpha = alpha.t, Qb = Qb.t, curll = curll,
                             acc = acc.beta1, att = att.beta1, MH = MH.beta1)

  beta1.int  <- this.update$beta.int
  SS1.int    <- this.update$SS.int
  beta1.time <- this.update$beta.time
  SS1.time   <- this.update$SS.time
  mu         <- this.update$mu
  curll      <- this.update$curll
  acc.beta1  <- this.update$acc
  att.beta1  <- this.update$att

  keep.beta1.int[iter, , ]  <- beta1.int
  keep.beta1.time[iter, , ] <- beta1.time

  if (iter < burn / 2) {
    this.update <- mhUpdate(acc = acc.beta1, att = att.beta1, MH = MH.beta1,
                            target.min = 0.5, target.max = 0.8,
                            nattempts = 50)
    acc.beta1 <- this.update$acc
    att.beta1 <- this.update$att
    MH.beta1  <- this.update$MH
  }

  if (iter %% 100 == 0) {
    par(mfrow = c(4, 4))
    acc.rate <- round(acc.beta1 / att.beta1, 3)
    for (i in 1:4) {
      plot(keep.beta1.int[1:iter, 1, i], type = "l",
           main = paste("intercept site 1, day ", i, " (",
                        round(beta1.int.t[1, i], 3), ")", sep = ""),
           ylab = paste("MH: ", MH.beta1[1, i], sep = ""),
           xlab = acc.rate[1, i])
      plot(keep.beta1.time[1:iter, 1, i], type = "l",
           main = paste("time site 1, day ", i, " (",
                        round(beta1.time.t[1, i], 3), ")", sep = ""),
           ylab = paste("MH: ", MH.beta1[1, i], sep = ""),
           xlab = acc.rate[1, i])
      plot(keep.beta1.int[1:iter, 2, i], type = "l",
           main = paste("intercept site 1, day ", i, " (",
                        round(beta1.int.t[2, i], 3), ")", sep = ""),
           ylab = paste("MH: ", MH.beta1[2, i], sep = ""),
           xlab = acc.rate[2, i])
      plot(keep.beta1.time[1:iter, 2, i], type = "l",
           main = paste("time site 1, day ", i, " (",
                        round(beta1.time.t[2, i], 3), ")", sep = ""),
           ylab = paste("MH: ", MH.beta1[2, i], sep = ""),
           xlab = acc.rate[2, i])
    }

  }
}
Rprof(NULL)
summaryRprof(filename = "Rprof.out", lines = "show")