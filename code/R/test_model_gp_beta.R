rm(list=ls())
library(fields)
library(Rcpp)
library(emulator)
library(microbenchmark)
library(SpatialExtremes)
library(MASS)
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

tau.int.t  <- c(rgamma(1, 1, 1), rgamma(1, 1, 1))
tau.time.t <- c(rgamma(1, 1, 1), rgamma(1, 1, 1))
beta.int.mn.t  <- rnorm(2)
beta.time.mn.t <- rnorm(2)
beta.int.t <- beta.time.t <- matrix(0, ns, 2)
for (p in 1:2) {
  beta.int.t[, p] <- beta.int.mn.t[p] +
    t(chol(Sigma.t)) %*% rnorm(ns, 0, 1 / sqrt(tau.int.t[p]))
  beta.time.t[, p] <- beta.time.mn.t[p] +
    t(chol(Sigma.t)) %*% rnorm(ns, 0, 1 / sqrt(tau.time.t[p]))
}

mu.t <- ls.t <- matrix(0, ns, nt)
time <- (1:nt - nt / 2) / nt
for (t in 1:nt) {
  mu.t[, t] <- beta.int.t[, 1] + beta.time.t[, 1] * time[t]
  ls.t[, t] <- beta.int.t[, 2] + beta.time.t[, 2] * time[t]
}

SS.int.t <- SS.time.t <- rep(0, 2)
for (p in 1:2) {
  SS.int.t[p]  <- quad.form(Qb.t, beta.int.t[, p] - beta.int.mn.t[p])
  SS.time.t[p] <- quad.form(Qb.t, beta.time.t[, p] - beta.time.mn.t[p])
}

xi.t <- 0
y.t <- rgev(n = ns * nt, loc = mu.t, exp(ls.t), xi.t)

thresh.t <- matrix(-Inf, ns, nt)
theta.t <- matrix(1, ns, nt)
theta.xi.t <- theta.t^xi.t
alpha.t <- 1

lp.beta1.int <- logpost.beta1.int(beta.int = beta.int.t[, 1],
                                  beta.mn = beta.int.mn.t[1],
                                  tau = tau.int.t[1], Qb = Qb.t,
                                  beta.time = beta.time.t[, 1], time = time,
                                  y = y.t, ls = ls.t, xi = xi.t,
                                  theta = theta.t, theta.xi = theta.xi.t,
                                  thresh = thresh.t, alpha = alpha.t)

mean(as.vector(jacobian(func = logpost.beta1.int, x = beta.int.t[, 1],
          beta.mn = beta.int.mn.t[1], tau = tau.int.t[1], Qb = Qb.t,
          beta.time = beta.time.t[, 1], time = time,
          y = y.t, ls = ls.t, xi = xi.t,
          theta = theta.t, theta.xi = theta.xi.t,
          thresh = thresh.t,
          alpha = alpha.t)) /
       logpost.beta1.int.grad(beta.int = beta.int.t[, 1],
                         beta.mn = beta.int.mn.t[1],
                         tau = tau.int.t[1], Qb = Qb.t,
                         beta.time = beta.time.t[, 1], time = time,
                         y = y.t, ls = ls.t, xi = xi.t,
                         theta = theta.t, theta.xi = theta.xi.t,
                         thresh = thresh.t, alpha = alpha.t))

sd(as.vector(jacobian(func = logpost.beta1.int, x = beta.int.t[, 1],
                      beta.mn = beta.int.mn.t[1], tau = tau.int.t[1], Qb = Qb.t,
                      beta.time = beta.time.t[, 1], time = time,
                      y = y.t, ls = ls.t, xi = xi.t,
                      theta = theta.t, theta.xi = theta.xi.t,
                      thresh = thresh.t,
                      alpha = alpha.t)) /
     logpost.beta1.int.grad(beta.int = beta.int.t[, 1],
                            beta.mn = beta.int.mn.t[1],
                            tau = tau.int.t[1], Qb = Qb.t,
                            beta.time = beta.time.t[, 1], time = time,
                            y = y.t, ls = ls.t, xi = xi.t,
                            theta = theta.t, theta.xi = theta.xi.t,
                            thresh = thresh.t, alpha = alpha.t))


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

tau.int.t  <- c(rgamma(1, 1, 1), rgamma(1, 1, 1))
tau.time.t <- c(rgamma(1, 1, 1), rgamma(1, 1, 1))
beta.int.mn.t  <- rnorm(2)
beta.time.mn.t <- rnorm(2)
beta.int.t <- beta.time.t <- matrix(0, ns, 2)
for (p in 1:2) {
  beta.int.t[, p] <- beta.int.mn.t[p] +
    t(chol(Sigma.t)) %*% rnorm(ns, 0, 1 / sqrt(tau.int.t[p]))
  beta.time.t[, p] <- beta.time.mn.t[p] +
    t(chol(Sigma.t)) %*% rnorm(ns, 0, 1 / sqrt(tau.time.t[p]))
}

mu.t <- ls.t <- matrix(0, ns, nt)
time <- (1:nt - nt / 2) / nt
for (t in 1:nt) {
  mu.t[, t] <- beta.int.t[, 1] + beta.time.t[, 1] * time[t]
  ls.t[, t] <- beta.int.t[, 2] + beta.time.t[, 2] * time[t]
}

SS.int.t <- SS.time.t <- rep(0, 2)
for (p in 1:2) {
  SS.int.t[p]  <- quad.form(Qb.t, beta.int.t[, p] - beta.int.mn.t[p])
  SS.time.t[p] <- quad.form(Qb.t, beta.time.t[, p] - beta.time.mn.t[p])
}

xi.t <- 0
y.t <- rgev(n = ns * nt, loc = mu.t, exp(ls.t), xi.t)

thresh.t <- matrix(-Inf, ns, nt)
theta.t <- matrix(1, ns, nt)
theta.xi.t <- theta.t^xi.t
alpha.t <- 1

# initialize values
beta.int.mn <- beta.time.mn <- rep(0, 2)
SS.int <- SS.time <- rep(0, 2)
beta.pri.sd  <- c(100, 10)

curll <- loglike(y = y.t, mu = mu.t, ls = ls.t, xi = xi.t,
                 theta = theta.t, theta.xi = theta.xi.t,
                 thresh = thresh.t, alpha = alpha.t)

niters <- 10000
burn   <- 8000
keep.beta.int.mn <- keep.beta.time.mn <- matrix(0, niters, 2)

set.seed(3366)  # demo
for (iter in 1:niters) {
  for (p in 1:2) {
    this.update <- updateGPMean(beta.sd = beta.pri.sd[p], Qb = Qb.t,
                                beta.int = beta.int.t[, p],
                                tau.int = tau.int.t[p],
                                beta.time = beta.time.t[, p],
                                tau.time = tau.time.t[p])
    beta.int.mn[p]  <- this.update$beta.int.mn
    SS.int[p] <- this.update$SS.int
    beta.time.mn[p] <- this.update$beta.time.mn
    SS.time[p] <- this.update$SS.time
  }

  keep.beta.int.mn[iter, ] <- beta.int.mn
  keep.beta.time.mn[iter, ] <- beta.time.mn

  if (iter %% 1000 == 0) {
    par(mfrow = c(2, 2))
    for (p in 1:2) {
      plot(keep.beta.int.mn[1:iter, p], type = "l",
           main = paste("Beta ", p, " intercept mean: ",
                        beta.int.mn.t[p], sep = ""))
      plot(keep.beta.time.mn[1:iter, p], type = "l",
           main = paste("Beta ", p, " time mean: ",
                        beta.time.mn.t[p], sep = ""))
    }
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

tau.int.t  <- c(rgamma(1, 1, 1), rgamma(1, 1, 1))
tau.time.t <- c(rgamma(1, 1, 1), rgamma(1, 1, 1))
beta.int.mn.t  <- rnorm(2)
beta.time.mn.t <- rnorm(2)
beta.int.t <- beta.time.t <- matrix(0, ns, 2)
for (p in 1:2) {
  beta.int.t[, p] <- beta.int.mn.t[p] +
    t(chol(Sigma.t)) %*% rnorm(ns, 0, 1 / sqrt(tau.int.t[p]))
  beta.time.t[, p] <- beta.time.mn.t[p] +
    t(chol(Sigma.t)) %*% rnorm(ns, 0, 1 / sqrt(tau.time.t[p]))
}

mu.t <- ls.t <- matrix(0, ns, nt)
time <- (1:nt - nt / 2) / nt
for (t in 1:nt) {
  mu.t[, t] <- beta.int.t[, 1] + beta.time.t[, 1] * time[t]
  ls.t[, t] <- beta.int.t[, 2] + beta.time.t[, 2] * time[t]
}

SS.int.t <- SS.time.t <- rep(0, 2)
for (p in 1:2) {
  SS.int.t[p]  <- quad.form(Qb.t, beta.int.t[, p] - beta.int.mn.t[p])
  SS.time.t[p] <- quad.form(Qb.t, beta.time.t[, p] - beta.time.mn.t[p])
}

xi.t <- 0
y.t <- rgev(n = ns * nt, loc = mu.t, exp(ls.t), xi.t)

thresh.t <- matrix(-Inf, ns, nt)
theta.t <- matrix(1, ns, nt)
theta.xi.t <- theta.t^xi.t
alpha.t <- 1

# initialize values
tau.int <- tau.time <- rep(0, 2)

curll <- loglike(y = y.t, mu = mu.t, ls = ls.t, xi = xi.t,
                 theta = theta.t, theta.xi = theta.xi.t,
                 thresh = thresh.t, alpha = alpha.t)

niters <- 10000
burn   <- 8000
keep.tau.int <- keep.tau.time <- matrix(0, niters, 2)

set.seed(3366)  # demo
for (iter in 1:niters) {
  for (p in 1:2) {
    this.update <- updateGPTau(SS.int = SS.int.t[p],
                               SS.time = SS.time.t[p],
                               tau.a = 0.1, tau.b = 0.1, ns = ns)
    tau.int[p]  <- this.update$tau.int
    tau.time[p] <- this.update$tau.time
  }

  keep.tau.int[iter, ]  <- tau.int
  keep.tau.time[iter, ] <- tau.time

  if (iter %% 1000 == 0) {
    par(mfrow = c(2, 2))
    for (p in 1:2) {
      plot(keep.tau.int[1:iter, p], type = "l",
           main = paste("Tau int ", p, " : ", tau.int.t[p], sep = ""))
      plot(keep.tau.time[1:iter, p], type = "l",
           main = paste("Tau time ", p, " : ", tau.time.t[p], sep = ""))
    }
  }
}

#### testing beta.mu ####
rm(list=ls())
source("../../../usefulR/usefulfunctions.R", chdir = TRUE)
openblas.set.num.threads(3)
source("auxfunctions.R", chdir = TRUE)
source("updatemodel.R", chdir = TRUE)

set.seed(2000)
ns <- 30
nt <- 5

phi.t <- 0.2
s <- cbind(runif(ns), runif(ns))
d <- rdist(s)
Sigma.t <- exp(-d / phi.t)
Qb.t    <- chol2inv(chol(Sigma.t))

tau.int.t  <- c(0.25, 0.25)
tau.time.t <- c(0.25, 0.25)
beta.int.mn.t  <- rep(0, 2)
beta.time.mn.t <- rep(0, 2)
beta.int.t <- beta.time.t <- matrix(0, ns, 2)
for (p in 1:2) {
  beta.int.t[, p] <- mvrnorm(1, mu = rep(beta.int.mn.t[p], ns),
                             Sigma = Sigma.t / tau.int.t[p])
  beta.time.t[, p] <- mvrnorm(1, mu = rep(beta.time.mn.t[p], ns),
                              Sigma = Sigma.t / tau.time.t[p])
}

mu.t <- ls.t <- matrix(0, ns, nt)
time <- (1:nt - nt / 2) / nt
for (t in 1:nt) {
  mu.t[, t] <- beta.int.t[, 1] + beta.time.t[, 1] * time[t]
  ls.t[, t] <- beta.int.t[, 2] + beta.time.t[, 2] * time[t]
}

SS.int.t <- SS.time.t <- rep(0, 2)
for (p in 1:2) {
  SS.int.t[p]  <- quad.form(Qb.t, beta.int.t[, p] - beta.int.mn.t[p])
  SS.time.t[p] <- quad.form(Qb.t, beta.time.t[, p] - beta.time.mn.t[p])
}

xi.t <- 0.01
y.t <- matrix(0, ns, nt)
for (t in 1:nt) {
  y.t[, t] <- rgev(n = ns, loc = mu.t[, t], exp(ls.t[, t]), xi.t)
}

thresh.t <- matrix(-Inf, ns, nt)
theta.t <- matrix(1, ns, nt)
theta.xi.t <- theta.t^xi.t
alpha.t <- 1

# initialize values
beta.int  <- beta.int.t + rnorm(ns, 0, 0.5)
beta.time <- beta.time.t + rnorm(ns, 0, 0.5)
SS.int <- SS.time <- rep(0, 2)
for (p in 1:2) {
  SS.int[p]  <- quad.form(Qb.t, beta.int[, p] - beta.int.mn.t[p])
  SS.time[p] <- quad.form(Qb.t, beta.time[, p] - beta.time.mn.t[p])
}

mu <- ls <- matrix(0, ns, nt)
for (t in 1:nt) {
  mu[, t] <- beta.int[, 1] + beta.time[, 1] * time[t]
  ls[, t] <- beta.int[, 2] + beta.time[, 2] * time[t]
}

curll <- loglike(y = y.t, mu = mu.t, ls = ls, xi = xi.t,
                 theta = theta.t, theta.xi = theta.xi.t,
                 thresh = thresh.t, alpha = alpha.t)

niters <- 30000
burn   <- 25000
keep.beta.int <- keep.beta.time <- array(0, dim = c(niters, ns, 2))
keep.mu <- array(0, dim = c(niters, ns, nt))
acc.beta.int  <- att.beta.int  <- MH.beta.int  <- matrix(0.1, ns, 2)
acc.beta.time <- att.beta.time <- MH.beta.time <- matrix(0.1, ns, 2)
MH.beta.int[, 2] <- 0.1
MH.beta.time[, 2] <- 0.1

beta.int <- beta.int.t
mu <- mu.t
set.seed(3366)  # demo
for (iter in 1:niters) {
  this.update <- updateBeta1.int(beta.int = beta.int[, 1],
                                 beta.mn = beta.int.mn.t[1],
                                 SS = SS.int[1], tau = tau.int.t[1],
                                 beta.time = beta.time[, 1], time = time,
                                 y = y.t, theta = theta.t,
                                 theta.xi = theta.xi.t, mu = mu, ls = ls.t,
                                 xi = xi.t, thresh = thresh.t, alpha = alpha.t,
                                 Qb = Qb.t, curll = curll,
                                 acc = acc.beta.int[, 1],
                                 att = att.beta.int[, 1],
                                 MH = MH.beta.int[, 1])
  beta.int[, 1] <- this.update$beta.int
  SS.int[1]     <- this.update$SS
  mu            <- this.update$mu
  curll         <- this.update$curll
  acc.beta.int[, 1] <- this.update$acc
  att.beta.int[, 1] <- this.update$att

  this.update <- updateBeta1.time(beta.time = beta.time[, 1],
                                  beta.mn = beta.time.mn.t[1],
                                  SS = SS.time[1], tau = tau.time.t[1],
                                  beta.int = beta.int[, 1], time = time,
                                  y = y.t, theta = theta.t,
                                  theta.xi = theta.xi.t, mu = mu, ls = ls.t,
                                  xi = xi.t, thresh = thresh.t, alpha = alpha.t,
                                  Qb = Qb.t, curll = curll,
                                  acc = acc.beta.time[, 1],
                                  att = att.beta.time[, 1],
                                  MH = MH.beta.time[, 1])
  beta.time[, 1] <- this.update$beta.time
  SS.time[1]     <- this.update$SS
  mu             <- this.update$mu
  curll          <- this.update$curll
  acc.beta.time[, 1] <- this.update$acc
  att.beta.time[, 1] <- this.update$att

  keep.beta.int[iter, , ] <- beta.int
  keep.beta.time[iter, , ] <- beta.time
  keep.mu[iter, , ] <- mu

  if (iter < burn / 2) {
    this.update <- mhUpdate(acc = acc.beta.int[, 1], att = att.beta.int[, 1],
                            MH = MH.beta.int[, 1], target.min = 0.3,
                            target.max = 0.6, nattempts = 400)
    acc.beta.int[, 1] <- this.update$acc
    att.beta.int[, 1] <- this.update$att
    MH.beta.int[, 1]  <- this.update$MH

    this.update <- mhUpdate(acc = acc.beta.time[, 1], att = att.beta.time[, 1],
                            MH = MH.beta.time[, 1], target.min = 0.3,
                            target.max = 0.6, nattempts = 400)
    acc.beta.time[, 1] <- this.update$acc
    att.beta.time[, 1] <- this.update$att
    MH.beta.time[, 1]  <- this.update$MH
  }

  if (iter %% 1000 == 0) {
    if (iter < burn) {
      start <- max(1, iter - 5000)
    } else {
      start <- burn + 1
    }
    par(mfrow = c(4, 4))
    acc.rate.int <- round(acc.beta.int / att.beta.int, 3)
    acc.rate.time <- round(acc.beta.time / att.beta.time, 3)
    for (i in 1:4) {
      plot(keep.beta.int[start:iter, i, 1], type = "l",
           main = paste("beta int ", i, " : ", round(beta.int.t[i], 3),
                        sep = ""),
           ylab = paste("MH: ", round(MH.beta.int[i, 1], 3)),
           xlab = acc.rate.int[i, 1])
      plot(keep.beta.time[start:iter, i, 1], type = "l",
           main = paste("beta time ", i, " : ", round(beta.time.t[i], 3),
                        sep = ""),
           ylab = paste("MH: ", round(MH.beta.time[i, 1], 3)),
           xlab = acc.rate.time[i, 1])
    }
    for (i in 1:2) {
      for (j in 1:4) {
        plot(keep.mu[start:iter, i, j], type = "l",
             main = paste("mu ", i, ", ", j, ": ", round(mu.t[i, j], 3),
                          sep = ""))
      }
    }
  }
}

#### testing entire mu update ####
rm(list=ls())
source("../../../usefulR/usefulfunctions.R", chdir = TRUE)
openblas.set.num.threads(3)
source("auxfunctions.R", chdir = TRUE)
source("updatemodel.R", chdir = TRUE)

set.seed(2000)
ns <- 30
nt <- 5

phi.t <- 0.2
s <- cbind(runif(ns), runif(ns))
d <- rdist(s)
Sigma.t <- exp(-d / phi.t)
Qb.t    <- chol2inv(chol(Sigma.t))

tau.int.t  <- c(0.25, 0.25)
tau.time.t <- c(0.25, 0.25)
beta.int.mn.t  <- rep(0, 2)
beta.time.mn.t <- rep(0, 2)
beta.int.t <- beta.time.t <- matrix(0, ns, 2)
for (p in 1:2) {
  beta.int.t[, p] <- mvrnorm(1, mu = rep(beta.int.mn.t[p], ns),
                             Sigma = Sigma.t / tau.int.t[p])
  beta.time.t[, p] <- mvrnorm(1, mu = rep(beta.time.mn.t[p], ns),
                              Sigma = Sigma.t / tau.time.t[p])
}

mu.t <- ls.t <- matrix(0, ns, nt)
time <- (1:nt - nt / 2) / nt
for (t in 1:nt) {
  mu.t[, t] <- beta.int.t[, 1] + beta.time.t[, 1] * time[t]
  ls.t[, t] <- beta.int.t[, 2] + beta.time.t[, 2] * time[t]
}

SS.int.t <- SS.time.t <- rep(0, 2)
for (p in 1:2) {
  SS.int.t[p]  <- quad.form(Qb.t, beta.int.t[, p] - beta.int.mn.t[p])
  SS.time.t[p] <- quad.form(Qb.t, beta.time.t[, p] - beta.time.mn.t[p])
}

xi.t <- 0.01
y.t <- matrix(0, ns, nt)
for (t in 1:nt) {
  y.t[, t] <- rgev(n = ns, loc = mu.t[, t], exp(ls.t[, t]), xi.t)
}

thresh.t <- matrix(-Inf, ns, nt)
theta.t <- matrix(1, ns, nt)
theta.xi.t <- theta.t^xi.t
alpha.t <- 1

# initialize values
beta.int  <- beta.int.t + rnorm(ns, 0, 0.5)
beta.time <- beta.time.t + rnorm(ns, 0, 0.5)
beta.int.mn <- beta.time.mn <- rep(0, 2)
tau.int <- tau.time <- rep(1, 2)
SS.int <- SS.time <- rep(0, 2)
for (p in 1:2) {
  SS.int[p]  <- quad.form(Qb.t, beta.int[, p] - beta.int.mn[p])
  SS.time[p] <- quad.form(Qb.t, beta.time[, p] - beta.time.mn[p])
}
beta.pri.sd  <- c(100, 10)

mu <- ls <- matrix(0, ns, nt)
for (t in 1:nt) {
  mu[, t] <- beta.int[, 1] + beta.time[, 1] * time[t]
  ls[, t] <- beta.int[, 2] + beta.time[, 2] * time[t]
}

curll <- loglike(y = y.t, mu = mu, ls = ls.t, xi = xi.t,
                 theta = theta.t, theta.xi = theta.xi.t,
                 thresh = thresh.t, alpha = alpha.t)

niters <- 30000
burn   <- 25000
keep.beta.int <- keep.beta.time <- array(0, dim = c(niters, ns, 2))
keep.mu <- array(0, dim = c(niters, ns, nt))
keep.beta.int.mn <- keep.beta.time.mn <- matrix(0, niters, 2)
keep.tau.int <- keep.tau.time <- matrix(0, niters, 2)
acc.beta.int  <- att.beta.int  <- MH.beta.int  <- matrix(0.1, ns, 2)
acc.beta.time <- att.beta.time <- MH.beta.time <- matrix(0.1, ns, 2)

beta.int <- beta.int.t
mu <- mu.t
set.seed(3366)  # demo
for (iter in 1:niters) {
  this.update <- updateBeta1.int(beta.int = beta.int[, 1],
                                 beta.mn = beta.int.mn[1],
                                 SS = SS.int[1], tau = tau.int[1],
                                 beta.time = beta.time[, 1], time = time,
                                 y = y.t, theta = theta.t,
                                 theta.xi = theta.xi.t, mu = mu, ls = ls.t,
                                 xi = xi.t, thresh = thresh.t, alpha = alpha.t,
                                 Qb = Qb.t, curll = curll,
                                 acc = acc.beta.int[, 1],
                                 att = att.beta.int[, 1],
                                 MH = MH.beta.int[, 1])
  beta.int[, 1] <- this.update$beta.int
  SS.int[1]     <- this.update$SS
  mu            <- this.update$mu
  curll         <- this.update$curll
  acc.beta.int[, 1] <- this.update$acc
  att.beta.int[, 1] <- this.update$att

  this.update <- updateBeta1.time(beta.time = beta.time[, 1],
                                  beta.mn = beta.time.mn[1],
                                  SS = SS.time[1], tau = tau.time[1],
                                  beta.int = beta.int[, 1], time = time,
                                  y = y.t, theta = theta.t,
                                  theta.xi = theta.xi.t, mu = mu, ls = ls.t,
                                  xi = xi.t, thresh = thresh.t, alpha = alpha.t,
                                  Qb = Qb.t, curll = curll,
                                  acc = acc.beta.time[, 1],
                                  att = att.beta.time[, 1],
                                  MH = MH.beta.time[, 1])
  beta.time[, 1] <- this.update$beta.time
  SS.time[1]     <- this.update$SS
  mu             <- this.update$mu
  curll          <- this.update$curll
  acc.beta.time[, 1] <- this.update$acc
  att.beta.time[, 1] <- this.update$att

  keep.beta.int[iter, , ] <- beta.int
  keep.beta.time[iter, , ] <- beta.time
  keep.mu[iter, , ] <- mu

  for (p in 1:2) {
    this.update <- updateGPMean(beta.sd = beta.pri.sd[p], Qb = Qb.t,
                                beta.int = beta.int[, p],
                                tau.int = tau.int[p],
                                beta.time = beta.time[, p],
                                tau.time = tau.time[p])
    beta.int.mn[p]  <- this.update$beta.int.mn
    SS.int[p] <- this.update$SS.int
    beta.time.mn[p] <- this.update$beta.time.mn
    SS.time[p] <- this.update$SS.time
  }

  keep.beta.int.mn[iter, ] <- beta.int.mn
  keep.beta.time.mn[iter, ] <- beta.time.mn

  for (p in 1:2) {
    this.update <- updateGPTau(SS.int = SS.int[p],
                               SS.time = SS.time[p],
                               tau.a = 0.1, tau.b = 0.1, ns = ns)
    tau.int[p]  <- this.update$tau.int
    tau.time[p] <- this.update$tau.time
  }

  keep.tau.int[iter, ]  <- tau.int
  keep.tau.time[iter, ] <- tau.time


  if (iter < burn / 2) {
    this.update <- mhUpdate(acc = acc.beta.int[, 1], att = att.beta.int[, 1],
                            MH = MH.beta.int[, 1], target.min = 0.3,
                            target.max = 0.6, nattempts = 400)
    acc.beta.int[, 1] <- this.update$acc
    att.beta.int[, 1] <- this.update$att
    MH.beta.int[, 1]  <- this.update$MH

    this.update <- mhUpdate(acc = acc.beta.time[, 1], att = att.beta.time[, 1],
                            MH = MH.beta.time[, 1], target.min = 0.3,
                            target.max = 0.6, nattempts = 400)
    acc.beta.time[, 1] <- this.update$acc
    att.beta.time[, 1] <- this.update$att
    MH.beta.time[, 1]  <- this.update$MH
  }

  if (iter %% 1000 == 0) {
    if (iter < burn) {
      start <- max(1, iter - 5000)
    } else {
      start <- burn + 1
    }
    par(mfrow = c(4, 4))
    acc.rate.int <- round(acc.beta.int / att.beta.int, 3)
    acc.rate.time <- round(acc.beta.time / att.beta.time, 3)
    for (i in 1:4) {
      plot(keep.beta.int[start:iter, i, 1], type = "l",
           main = paste("beta int ", i, " : ", round(beta.int.t[i], 3),
                        sep = ""),
           ylab = paste("MH: ", round(MH.beta.int[i, 1], 3)),
           xlab = acc.rate.int[i, 1])
      plot(keep.beta.time[start:iter, i, 1], type = "l",
           main = paste("beta time ", i, " : ", round(beta.time.t[i], 3),
                        sep = ""),
           ylab = paste("MH: ", round(MH.beta.time[i, 1], 3)),
           xlab = acc.rate.time[i, 1])
    }
    for (i in 1) {
      for (j in 1:4) {
        plot(keep.mu[start:iter, i, j], type = "l",
             main = paste("mu ", i, ", ", j, ": ", round(mu.t[i, j], 3),
                          sep = ""))
      }
    }
    plot(keep.tau.int[start:iter, 1], type = "l",
         main = paste("tau int ", round(tau.int.t[1], 3), sep = ""))
    plot(keep.tau.time[start:iter, 1], type = "l",
         main = paste("tau time ", round(tau.time.t[1], 3), sep = ""))
    plot(keep.beta.int.mn[start:iter, 1], type = "l",
         main = paste("beta int mean ", round(beta.int.mn.t[1], 3), sep = ""))
    plot(keep.beta.time.mn[start:iter, 1], type = "l",
         main = paste("beta time mean ", round(beta.time.mn.t[1], 3), sep = ""))
  }
}

#### testing entire ls update ####
rm(list=ls())
source("../../../usefulR/usefulfunctions.R", chdir = TRUE)
openblas.set.num.threads(3)
source("auxfunctions.R", chdir = TRUE)
source("updatemodel.R", chdir = TRUE)

set.seed(2000)
ns <- 30
nt <- 5

phi.t <- 0.2
s <- cbind(runif(ns), runif(ns))
d <- rdist(s)
Sigma.t <- exp(-d / phi.t)
Qb.t    <- chol2inv(chol(Sigma.t))

tau.int.t  <- c(0.25, 0.25)
tau.time.t <- c(0.25, 0.25)
beta.int.mn.t  <- rep(0, 2)
beta.time.mn.t <- rep(0, 2)
beta.int.t <- beta.time.t <- matrix(0, ns, 2)
for (p in 1:2) {
  beta.int.t[, p] <- mvrnorm(1, mu = rep(beta.int.mn.t[p], ns),
                             Sigma = Sigma.t / tau.int.t[p])
  beta.time.t[, p] <- mvrnorm(1, mu = rep(beta.time.mn.t[p], ns),
                              Sigma = Sigma.t / tau.time.t[p])
}

mu.t <- ls.t <- matrix(0, ns, nt)
time <- (1:nt - nt / 2) / nt
for (t in 1:nt) {
  mu.t[, t] <- beta.int.t[, 1] + beta.time.t[, 1] * time[t]
  ls.t[, t] <- beta.int.t[, 2] + beta.time.t[, 2] * time[t]
}

SS.int.t <- SS.time.t <- rep(0, 2)
for (p in 1:2) {
  SS.int.t[p]  <- quad.form(Qb.t, beta.int.t[, p] - beta.int.mn.t[p])
  SS.time.t[p] <- quad.form(Qb.t, beta.time.t[, p] - beta.time.mn.t[p])
}

xi.t <- 0.01
y.t <- matrix(0, ns, nt)
for (t in 1:nt) {
  y.t[, t] <- rgev(n = ns, loc = mu.t[, t], exp(ls.t[, t]), xi.t)
}

thresh.t <- matrix(-Inf, ns, nt)
theta.t <- matrix(1, ns, nt)
theta.xi.t <- theta.t^xi.t
alpha.t <- 1

# initialize values
beta.int  <- beta.int.t + rnorm(ns, 0, 0.5)
beta.time <- beta.time.t + rnorm(ns, 0, 0.5)
beta.int.mn <- beta.time.mn <- rep(0, 2)
tau.int <- tau.time <- rep(1, 2)
SS.int <- SS.time <- rep(0, 2)
for (p in 1:2) {
  SS.int[p]  <- quad.form(Qb.t, beta.int[, p] - beta.int.mn[p])
  SS.time[p] <- quad.form(Qb.t, beta.time[, p] - beta.time.mn[p])
}
beta.pri.sd  <- c(100, 10)

mu <- ls <- matrix(0, ns, nt)
for (t in 1:nt) {
  mu[, t] <- beta.int[, 1] + beta.time[, 1] * time[t]
  ls[, t] <- beta.int[, 2] + beta.time[, 2] * time[t]
}

curll <- loglike(y = y.t, mu = mu, ls = ls.t, xi = xi.t,
                 theta = theta.t, theta.xi = theta.xi.t,
                 thresh = thresh.t, alpha = alpha.t)

niters <- 30000
burn   <- 25000
keep.beta.int <- keep.beta.time <- array(0, dim = c(niters, ns, 2))
keep.ls <- array(0, dim = c(niters, ns, nt))
keep.beta.int.mn <- keep.beta.time.mn <- matrix(0, niters, 2)
keep.tau.int <- keep.tau.time <- matrix(0, niters, 2)
acc.beta.int  <- att.beta.int  <- MH.beta.int  <- matrix(0.1, ns, 2)
acc.beta.time <- att.beta.time <- MH.beta.time <- matrix(0.1, ns, 2)

beta.int <- beta.int.t
mu <- mu.t
set.seed(3366)  # demo
for (iter in 1:niters) {
  p <- 2

  this.update <- updateBeta2.int(beta.int = beta.int[, p],
                                 beta.mn = beta.int.mn[p],
                                 SS = SS.int[p], tau = tau.int[p],
                                 beta.time = beta.time[, p], time = time,
                                 y = y.t, theta = theta.t,
                                 theta.xi = theta.xi.t, mu = mu.t, ls = ls,
                                 xi = xi.t, thresh = thresh.t, alpha = alpha.t,
                                 Qb = Qb.t, curll = curll,
                                 acc = acc.beta.int[, p],
                                 att = att.beta.int[, p],
                                 MH = MH.beta.int[, p])
  beta.int[, p] <- this.update$beta.int
  SS.int[p]     <- this.update$SS
  if (p == 1) {
    mu <- this.update$mu
  } else {
    ls <- this.update$ls
  }
  curll         <- this.update$curll
  acc.beta.int[, p] <- this.update$acc
  att.beta.int[, p] <- this.update$att

  this.update <- updateBeta2.time(beta.time = beta.time[, p],
                                  beta.mn = beta.time.mn[p],
                                  SS = SS.time[p], tau = tau.time[p],
                                  beta.int = beta.int[, p], time = time,
                                  y = y.t, theta = theta.t,
                                  theta.xi = theta.xi.t, mu = mu.t, ls = ls,
                                  xi = xi.t, thresh = thresh.t, alpha = alpha.t,
                                  Qb = Qb.t, curll = curll,
                                  acc = acc.beta.time[, p],
                                  att = att.beta.time[, p],
                                  MH = MH.beta.time[, p])
  beta.time[, p] <- this.update$beta.time
  SS.time[p]     <- this.update$SS
  if (p == 1) {
    mu <- this.update$mu
  } else {
    ls <- this.update$ls
  }
  curll          <- this.update$curll
  acc.beta.time[, p] <- this.update$acc
  att.beta.time[, p] <- this.update$att

  keep.beta.int[iter, , ] <- beta.int
  keep.beta.time[iter, , ] <- beta.time
  keep.ls[iter, , ] <- ls

  if (iter < burn / 2) {
    this.update <- mhUpdate(acc = acc.beta.int[, p], att = att.beta.int[, p],
                            MH = MH.beta.int[, 1], target.min = 0.3,
                            target.max = 0.6, nattempts = 400)
    acc.beta.int[, p] <- this.update$acc
    att.beta.int[, p] <- this.update$att
    MH.beta.int[, p]  <- this.update$MH

    this.update <- mhUpdate(acc = acc.beta.time[, p], att = att.beta.time[, p],
                            MH = MH.beta.time[, p], target.min = 0.3,
                            target.max = 0.6, nattempts = 400)
    acc.beta.time[, p] <- this.update$acc
    att.beta.time[, p] <- this.update$att
    MH.beta.time[, p]  <- this.update$MH
  }

  for (p in 1:2) {
    this.update <- updateGPMean(beta.sd = beta.pri.sd[p], Qb = Qb.t,
                                beta.int = beta.int[, p],
                                tau.int = tau.int[p],
                                beta.time = beta.time[, p],
                                tau.time = tau.time[p])
    beta.int.mn[p]  <- this.update$beta.int.mn
    SS.int[p] <- this.update$SS.int
    beta.time.mn[p] <- this.update$beta.time.mn
    SS.time[p] <- this.update$SS.time
  }

  keep.beta.int.mn[iter, ] <- beta.int.mn
  keep.beta.time.mn[iter, ] <- beta.time.mn

  for (p in 1:2) {
    this.update <- updateGPTau(SS.int = SS.int[p],
                               SS.time = SS.time[p],
                               tau.a = 0.1, tau.b = 0.1, ns = ns)
    tau.int[p]  <- this.update$tau.int
    tau.time[p] <- this.update$tau.time
  }

  keep.tau.int[iter, ]  <- tau.int
  keep.tau.time[iter, ] <- tau.time


  if (iter %% 1000 == 0) {
    if (iter < burn) {
      start <- max(1, iter - 5000)
    } else {
      start <- burn + 1
    }
    par(mfrow = c(4, 4))
    acc.rate.int <- round(acc.beta.int / att.beta.int, 3)
    acc.rate.time <- round(acc.beta.time / att.beta.time, 3)
    for (i in 1:4) {
      plot(keep.beta.int[start:iter, i, 2], type = "l",
           main = paste("beta int ", i, " : ", round(beta.int.t[i, 2], 3),
                        sep = ""),
           ylab = paste("MH: ", round(MH.beta.int[i, 2], 3)),
           xlab = acc.rate.int[i, 2])
      plot(keep.beta.time[start:iter, i, 2], type = "l",
           main = paste("beta time ", i, " : ", round(beta.time.t[i, 2], 3),
                        sep = ""),
           ylab = paste("MH: ", round(MH.beta.time[i, 2], 3)),
           xlab = acc.rate.time[i, 2])
    }
    for (i in 1) {
      for (j in 1:4) {
        plot(keep.ls[start:iter, i, j], type = "l",
             main = paste("ls ", i, ", ", j, ": ", round(ls.t[i, j], 3),
                          sep = ""))
      }
    }
    plot(keep.tau.int[start:iter, 2], type = "l",
         main = paste("tau int ", round(tau.int.t[2], 3), sep = ""))
    plot(keep.tau.time[start:iter, 2], type = "l",
         main = paste("tau time ", round(tau.time.t[2], 3), sep = ""))
    plot(keep.beta.int.mn[start:iter, 2], type = "l",
         main = paste("beta int mean ", round(beta.int.mn.t[2], 3), sep = ""))
    plot(keep.beta.time.mn[start:iter, 2], type = "l",
         main = paste("beta time mean ", round(beta.time.mn.t[2], 3), sep = ""))
  }
}

#### testing both mu and ls update ####
rm(list=ls())
source("../../../usefulR/usefulfunctions.R", chdir = TRUE)
openblas.set.num.threads(3)
source("auxfunctions.R", chdir = TRUE)
source("updatemodel.R", chdir = TRUE)

set.seed(2000)
ns <- 30
nt <- 5

phi.t <- 0.2
s <- cbind(runif(ns), runif(ns))
d <- rdist(s)
Sigma.t <- exp(-d / phi.t)
Qb.t    <- chol2inv(chol(Sigma.t))

tau.int.t  <- c(0.25, 0.25)
tau.time.t <- c(0.25, 0.25)
beta.int.mn.t  <- rep(0, 2)
beta.time.mn.t <- rep(0, 2)
beta.int.t <- beta.time.t <- matrix(0, ns, 2)
for (p in 1:2) {
  beta.int.t[, p] <- mvrnorm(1, mu = rep(beta.int.mn.t[p], ns),
                             Sigma = Sigma.t / tau.int.t[p])
  beta.time.t[, p] <- mvrnorm(1, mu = rep(beta.time.mn.t[p], ns),
                              Sigma = Sigma.t / tau.time.t[p])
}

mu.t <- ls.t <- matrix(0, ns, nt)
time <- (1:nt - nt / 2) / nt
for (t in 1:nt) {
  mu.t[, t] <- beta.int.t[, 1] + beta.time.t[, 1] * time[t]
  ls.t[, t] <- beta.int.t[, 2] + beta.time.t[, 2] * time[t]
}

SS.int.t <- SS.time.t <- rep(0, 2)
for (p in 1:2) {
  SS.int.t[p]  <- quad.form(Qb.t, beta.int.t[, p] - beta.int.mn.t[p])
  SS.time.t[p] <- quad.form(Qb.t, beta.time.t[, p] - beta.time.mn.t[p])
}

xi.t <- 0.01
y.t <- matrix(0, ns, nt)
for (t in 1:nt) {
  y.t[, t] <- rgev(n = ns, loc = mu.t[, t], exp(ls.t[, t]), xi.t)
}

thresh.t <- matrix(-Inf, ns, nt)
theta.t <- matrix(1, ns, nt)
theta.xi.t <- theta.t^xi.t
alpha.t <- 1

# initialize values
beta.int  <- beta.int.t + rnorm(ns, 0, 0.5)
beta.time <- beta.time.t + rnorm(ns, 0, 0.5)
beta.int.mn <- beta.time.mn <- rep(0, 2)
tau.int <- tau.time <- rep(1, 2)
SS.int <- SS.time <- rep(0, 2)
for (p in 1:2) {
  SS.int[p]  <- quad.form(Qb.t, beta.int[, p] - beta.int.mn[p])
  SS.time[p] <- quad.form(Qb.t, beta.time[, p] - beta.time.mn[p])
}
beta.pri.sd  <- c(100, 10)

mu <- ls <- matrix(0, ns, nt)
for (t in 1:nt) {
  mu[, t] <- beta.int[, 1] + beta.time[, 1] * time[t]
  ls[, t] <- beta.int[, 2] + beta.time[, 2] * time[t]
}

curll <- loglike(y = y.t, mu = mu, ls = ls.t, xi = xi.t,
                 theta = theta.t, theta.xi = theta.xi.t,
                 thresh = thresh.t, alpha = alpha.t)

niters <- 30000
burn   <- 25000
keep.beta.int <- keep.beta.time <- array(0, dim = c(niters, ns, 2))
keep.mu <- keep.ls <- array(0, dim = c(niters, ns, nt))
keep.beta.int.mn <- keep.beta.time.mn <- matrix(0, niters, 2)
keep.tau.int <- keep.tau.time <- matrix(0, niters, 2)
acc.beta.int  <- att.beta.int  <- MH.beta.int  <- matrix(0.1, ns, 2)
acc.beta.time <- att.beta.time <- MH.beta.time <- matrix(0.1, ns, 2)

beta.int <- beta.int.t
mu <- mu.t
set.seed(3366)  # demo
for (iter in 1:niters) {

  this.update <- updateBeta1.int(beta.int = beta.int[, 1],
                                 beta.mn = beta.int.mn[1],
                                 SS = SS.int[1], tau = tau.int[1],
                                 beta.time = beta.time[, 1], time = time,
                                 y = y.t, theta = theta.t,
                                 theta.xi = theta.xi.t, mu = mu, ls = ls,
                                 xi = xi.t, thresh = thresh.t, alpha = alpha.t,
                                 Qb = Qb.t, curll = curll,
                                 acc = acc.beta.int[, 1],
                                 att = att.beta.int[, 1],
                                 MH = MH.beta.int[, 1])
  beta.int[, 1] <- this.update$beta.int
  SS.int[1]     <- this.update$SS
  mu <- this.update$mu
  curll         <- this.update$curll
  acc.beta.int[, 1] <- this.update$acc
  att.beta.int[, 1] <- this.update$att

  this.update <- updateBeta1.time(beta.time = beta.time[, 1],
                                  beta.mn = beta.time.mn[1],
                                  SS = SS.time[1], tau = tau.time[1],
                                  beta.int = beta.int[, 1], time = time,
                                  y = y.t, theta = theta.t,
                                  theta.xi = theta.xi.t, mu = mu, ls = ls,
                                  xi = xi.t, thresh = thresh.t, alpha = alpha.t,
                                  Qb = Qb.t, curll = curll,
                                  acc = acc.beta.time[, 1],
                                  att = att.beta.time[, 1],
                                  MH = MH.beta.time[, 1])
  beta.time[, 1] <- this.update$beta.time
  SS.time[1]     <- this.update$SS
  mu             <- this.update$mu
  curll          <- this.update$curll
  acc.beta.time[, 1] <- this.update$acc
  att.beta.time[, 1] <- this.update$att

  this.update <- updateBeta2.int(beta.int = beta.int[, 2],
                                 beta.mn = beta.int.mn[2],
                                 SS = SS.int[2], tau = tau.int[2],
                                 beta.time = beta.time[, 2], time = time,
                                 y = y.t, theta = theta.t,
                                 theta.xi = theta.xi.t, mu = mu, ls = ls,
                                 xi = xi.t, thresh = thresh.t, alpha = alpha.t,
                                 Qb = Qb.t, curll = curll,
                                 acc = acc.beta.int[, 2],
                                 att = att.beta.int[, 2],
                                 MH = MH.beta.int[, 2])
  beta.int[, 2] <- this.update$beta.int
  SS.int[2]     <- this.update$SS
  ls            <- this.update$ls
  curll         <- this.update$curll
  acc.beta.int[, 2] <- this.update$acc
  att.beta.int[, 2] <- this.update$att

  this.update <- updateBeta2.time(beta.time = beta.time[, 2],
                                  beta.mn = beta.time.mn[2],
                                  SS = SS.time[2], tau = tau.time[2],
                                  beta.int = beta.int[, 2], time = time,
                                  y = y.t, theta = theta.t,
                                  theta.xi = theta.xi.t, mu = mu, ls = ls,
                                  xi = xi.t, thresh = thresh.t, alpha = alpha.t,
                                  Qb = Qb.t, curll = curll,
                                  acc = acc.beta.time[, 2],
                                  att = att.beta.time[, 2],
                                  MH = MH.beta.time[, 2])
  beta.time[, 2] <- this.update$beta.time
  SS.time[2]     <- this.update$SS
  ls             <- this.update$ls
  curll          <- this.update$curll
  acc.beta.time[, 2] <- this.update$acc
  att.beta.time[, 2] <- this.update$att

  keep.beta.int[iter, , ] <- beta.int
  keep.beta.time[iter, , ] <- beta.time
  keep.mu[iter, , ] <- mu
  keep.ls[iter, , ] <- ls

  if (iter < burn / 2) {
    for (p in 1:2) {
      this.update <- mhUpdate(acc = acc.beta.int[, p], att = att.beta.int[, p],
                              MH = MH.beta.int[, 1], target.min = 0.3,
                              target.max = 0.6, nattempts = 400)
      acc.beta.int[, p] <- this.update$acc
      att.beta.int[, p] <- this.update$att
      MH.beta.int[, p]  <- this.update$MH

      this.update <- mhUpdate(acc = acc.beta.time[, p],
                              att = att.beta.time[, p],
                              MH = MH.beta.time[, p], target.min = 0.3,
                              target.max = 0.6, nattempts = 400)
      acc.beta.time[, p] <- this.update$acc
      att.beta.time[, p] <- this.update$att
      MH.beta.time[, p]  <- this.update$MH
    }
  }

  for (p in 1:2) {
    this.update <- updateGPMean(beta.sd = beta.pri.sd[p], Qb = Qb.t,
                                beta.int = beta.int[, p],
                                tau.int = tau.int[p],
                                beta.time = beta.time[, p],
                                tau.time = tau.time[p])
    beta.int.mn[p]  <- this.update$beta.int.mn
    SS.int[p] <- this.update$SS.int
    beta.time.mn[p] <- this.update$beta.time.mn
    SS.time[p] <- this.update$SS.time
  }

  keep.beta.int.mn[iter, ] <- beta.int.mn
  keep.beta.time.mn[iter, ] <- beta.time.mn

  for (p in 1:2) {
    this.update <- updateGPTau(SS.int = SS.int[p],
                               SS.time = SS.time[p],
                               tau.a = 0.1, tau.b = 0.1, ns = ns)
    tau.int[p]  <- this.update$tau.int
    tau.time[p] <- this.update$tau.time
  }

  keep.tau.int[iter, ]  <- tau.int
  keep.tau.time[iter, ] <- tau.time


  if (iter %% 1000 == 0) {
    if (iter < burn) {
      start <- max(1, iter - 5000)
    } else {
      start <- burn + 1
    }
    params <- c("mu", "ls")
    par(mfrow = c(5, 4))
    acc.rate.int <- round(acc.beta.int / att.beta.int, 3)
    acc.rate.time <- round(acc.beta.time / att.beta.time, 3)
    for (i in 1:2) {
      for (p in 1:2) {
        plot(keep.beta.int[start:iter, i, p], type = "l",
             main = paste(params[p], " beta int ", i, " : ",
                          round(beta.int.t[i, p], 3), sep = ""),
             ylab = paste("MH: ", round(MH.beta.int[i, p], 3)),
             xlab = acc.rate.int[i, p])
        plot(keep.beta.time[start:iter, i, p], type = "l",
             main = paste(params[p], " beta time ", i, " : ",
                          round(beta.time.t[i, p], 3), sep = ""),
             ylab = paste("MH: ", round(MH.beta.time[i, p], 3)),
             xlab = acc.rate.time[i, p])
      }
    }
    for (i in 1) {
      for (j in 1:2) {
        plot(keep.mu[start:iter, i, j], type = "l",
             main = paste("mu ", i, ", ", j, ": ", round(mu.t[i, j], 3),
                          sep = ""))
      }
      for (j in 1:2) {
        plot(keep.ls[start:iter, i, j], type = "l",
             main = paste("ls ", i, ", ", j, ": ", round(ls.t[i, j], 3),
                          sep = ""))
      }
    }

    for (p in 1:2) {
      plot(keep.tau.int[start:iter, p], type = "l",
           main = paste(params[p], " tau int ", round(tau.int.t[p], 3),
                        sep = ""))
      plot(keep.tau.time[start:iter, p], type = "l",
           main = paste(params[p], " tau time ", round(tau.time.t[p], 3),
                        sep = ""))
    }
    for (p in 1:2) {
      plot(keep.beta.int.mn[start:iter, p], type = "l",
           main = paste(params[p], " beta int mean ",
                        round(beta.int.mn.t[p], 3), sep = ""))
      plot(keep.beta.time.mn[start:iter, p], type = "l",
           main = paste(params[p], " beta time mean ",
                        round(beta.time.mn.t[p], 3), sep = ""))
    }
  }
}

#### testing mu, ls, and xi update ####
rm(list=ls())
source("../../../usefulR/usefulfunctions.R", chdir = TRUE)
openblas.set.num.threads(3)
source("auxfunctions.R", chdir = TRUE)
source("updatemodel.R", chdir = TRUE)

set.seed(2000)
ns <- 30
nt <- 5

phi.t <- 0.2
s <- cbind(runif(ns), runif(ns))
d <- rdist(s)
Sigma.t <- exp(-d / phi.t)
Qb.t    <- chol2inv(chol(Sigma.t))

tau.int.t  <- c(0.25, 0.25)
tau.time.t <- c(0.25, 0.25)
beta.int.mn.t  <- rep(0, 2)
beta.time.mn.t <- rep(0, 2)
beta.int.t <- beta.time.t <- matrix(0, ns, 2)
for (p in 1:2) {
  beta.int.t[, p] <- mvrnorm(1, mu = rep(beta.int.mn.t[p], ns),
                             Sigma = Sigma.t / tau.int.t[p])
  beta.time.t[, p] <- mvrnorm(1, mu = rep(beta.time.mn.t[p], ns),
                              Sigma = Sigma.t / tau.time.t[p])
}

mu.t <- ls.t <- matrix(0, ns, nt)
time <- (1:nt - nt / 2) / nt
for (t in 1:nt) {
  mu.t[, t] <- beta.int.t[, 1] + beta.time.t[, 1] * time[t]
  ls.t[, t] <- beta.int.t[, 2] + beta.time.t[, 2] * time[t]
}

SS.int.t <- SS.time.t <- rep(0, 2)
for (p in 1:2) {
  SS.int.t[p]  <- quad.form(Qb.t, beta.int.t[, p] - beta.int.mn.t[p])
  SS.time.t[p] <- quad.form(Qb.t, beta.time.t[, p] - beta.time.mn.t[p])
}

xi.t <- -0.3
y.t <- matrix(0, ns, nt)
for (t in 1:nt) {
  y.t[, t] <- rgev(n = ns, loc = mu.t[, t], exp(ls.t[, t]), xi.t)
}

thresh.t <- matrix(-Inf, ns, nt)
theta.t <- matrix(1, ns, nt)
theta.xi.t <- theta.t^xi.t
alpha.t <- 1

# initialize values
beta.int  <- beta.int.t + rnorm(ns, 0, 0.5)
beta.time <- beta.time.t + rnorm(ns, 0, 0.5)
beta.int.mn <- beta.time.mn <- rep(0, 2)
tau.int <- tau.time <- rep(1, 2)
SS.int <- SS.time <- rep(0, 2)
for (p in 1:2) {
  SS.int[p]  <- quad.form(Qb.t, beta.int[, p] - beta.int.mn[p])
  SS.time[p] <- quad.form(Qb.t, beta.time[, p] - beta.time.mn[p])
}
xi <- 0.1
theta.xi <- theta.t^xi
beta.pri.sd  <- c(100, 10)

mu <- ls <- matrix(0, ns, nt)
for (t in 1:nt) {
  mu[, t] <- beta.int[, 1] + beta.time[, 1] * time[t]
  ls[, t] <- beta.int[, 2] + beta.time[, 2] * time[t]
}

curll <- loglike(y = y.t, mu = mu, ls = ls, xi = xi,
                 theta = theta.t, theta.xi = theta.xi,
                 thresh = thresh.t, alpha = alpha.t)

niters <- 30000
burn   <- 25000
keep.beta.int <- keep.beta.time <- array(0, dim = c(niters, ns, 2))
keep.mu <- keep.ls <- array(0, dim = c(niters, ns, nt))
keep.beta.int.mn <- keep.beta.time.mn <- matrix(0, niters, 2)
keep.tau.int <- keep.tau.time <- matrix(0, niters, 2)
keep.xi      <- rep(0, niters)
acc.beta.int  <- att.beta.int  <- MH.beta.int  <- matrix(0.1, ns, 2)
acc.beta.time <- att.beta.time <- MH.beta.time <- matrix(0.1, ns, 2)
acc.xi        <- att.xi        <- MH.xi        <- 0.1

set.seed(3366)  # demo
for (iter in 1:niters) {
  this.update <- updateBeta1.int(beta.int = beta.int[, 1],
                                 beta.mn = beta.int.mn[1],
                                 SS = SS.int[1], tau = tau.int[1],
                                 beta.time = beta.time[, 1], time = time,
                                 y = y.t, theta = theta.t,
                                 theta.xi = theta.xi, mu = mu, ls = ls,
                                 xi = xi, thresh = thresh.t, alpha = alpha.t,
                                 Qb = Qb.t, curll = curll,
                                 acc = acc.beta.int[, 1],
                                 att = att.beta.int[, 1],
                                 MH = MH.beta.int[, 1])
  beta.int[, 1] <- this.update$beta.int
  SS.int[1]     <- this.update$SS
  mu <- this.update$mu
  curll         <- this.update$curll
  acc.beta.int[, 1] <- this.update$acc
  att.beta.int[, 1] <- this.update$att

  this.update <- updateBeta1.time(beta.time = beta.time[, 1],
                                  beta.mn = beta.time.mn[1],
                                  SS = SS.time[1], tau = tau.time[1],
                                  beta.int = beta.int[, 1], time = time,
                                  y = y.t, theta = theta.t,
                                  theta.xi = theta.xi, mu = mu, ls = ls,
                                  xi = xi, thresh = thresh.t, alpha = alpha.t,
                                  Qb = Qb.t, curll = curll,
                                  acc = acc.beta.time[, 1],
                                  att = att.beta.time[, 1],
                                  MH = MH.beta.time[, 1])
  beta.time[, 1] <- this.update$beta.time
  SS.time[1]     <- this.update$SS
  mu             <- this.update$mu
  curll          <- this.update$curll
  acc.beta.time[, 1] <- this.update$acc
  att.beta.time[, 1] <- this.update$att

  this.update <- updateBeta2.int(beta.int = beta.int[, 2],
                                 beta.mn = beta.int.mn[2],
                                 SS = SS.int[2], tau = tau.int[2],
                                 beta.time = beta.time[, 2], time = time,
                                 y = y.t, theta = theta.t,
                                 theta.xi = theta.xi, mu = mu, ls = ls,
                                 xi = xi, thresh = thresh.t, alpha = alpha.t,
                                 Qb = Qb.t, curll = curll,
                                 acc = acc.beta.int[, 2],
                                 att = att.beta.int[, 2],
                                 MH = MH.beta.int[, 2])
  beta.int[, 2] <- this.update$beta.int
  SS.int[2]     <- this.update$SS
  ls            <- this.update$ls
  curll         <- this.update$curll
  acc.beta.int[, 2] <- this.update$acc
  att.beta.int[, 2] <- this.update$att

  this.update <- updateBeta2.time(beta.time = beta.time[, 2],
                                  beta.mn = beta.time.mn[2],
                                  SS = SS.time[2], tau = tau.time[2],
                                  beta.int = beta.int[, 2], time = time,
                                  y = y.t, theta = theta.t,
                                  theta.xi = theta.xi, mu = mu, ls = ls,
                                  xi = xi, thresh = thresh.t, alpha = alpha.t,
                                  Qb = Qb.t, curll = curll,
                                  acc = acc.beta.time[, 2],
                                  att = att.beta.time[, 2],
                                  MH = MH.beta.time[, 2])
  beta.time[, 2] <- this.update$beta.time
  SS.time[2]     <- this.update$SS
  ls             <- this.update$ls
  curll          <- this.update$curll
  acc.beta.time[, 2] <- this.update$acc
  att.beta.time[, 2] <- this.update$att

  this.update <- updateXi(xi = xi, xi.min = -0.5, xi.max = 0.5,
                          xi.mn = 0, xi.sd = 0.3, y = y.t, mu = mu, ls = ls,
                          curll = curll, theta = theta.t,
                          theta.xi = theta.xi, thresh = thresh.t,
                          alpha = alpha.t,
                          acc = acc.xi, att = att.xi, MH = MH.xi)
  xi       <- this.update$xi
  theta.xi <- this.update$theta.xi
  curll    <- this.update$curll
  acc.xi   <- this.update$acc
  att.xi   <- this.update$att

  keep.beta.int[iter, , ] <- beta.int
  keep.beta.time[iter, , ] <- beta.time
  keep.mu[iter, , ] <- mu
  keep.ls[iter, , ] <- ls
  keep.xi[iter]     <- xi

  if (iter < burn / 2) {
    for (p in 1:2) {
      this.update <- mhUpdate(acc = acc.beta.int[, p], att = att.beta.int[, p],
                              MH = MH.beta.int[, 1], target.min = 0.3,
                              target.max = 0.6, nattempts = 200)
      acc.beta.int[, p] <- this.update$acc
      att.beta.int[, p] <- this.update$att
      MH.beta.int[, p]  <- this.update$MH

      this.update <- mhUpdate(acc = acc.beta.time[, p],
                              att = att.beta.time[, p],
                              MH = MH.beta.time[, p], target.min = 0.3,
                              target.max = 0.6, nattempts = 200)
      acc.beta.time[, p] <- this.update$acc
      att.beta.time[, p] <- this.update$att
      MH.beta.time[, p]  <- this.update$MH

      this.update <- mhUpdate(acc = acc.xi, att = att.xi, MH = MH.xi,
                              target.min = 0.3, target.max = 0.6,
                              nattempts = 50)
    }
  }

  for (p in 1:2) {
    this.update <- updateGPMean(beta.sd = beta.pri.sd[p], Qb = Qb.t,
                                beta.int = beta.int[, p],
                                tau.int = tau.int[p],
                                beta.time = beta.time[, p],
                                tau.time = tau.time[p])
    beta.int.mn[p]  <- this.update$beta.int.mn
    SS.int[p] <- this.update$SS.int
    beta.time.mn[p] <- this.update$beta.time.mn
    SS.time[p] <- this.update$SS.time
  }

  keep.beta.int.mn[iter, ] <- beta.int.mn
  keep.beta.time.mn[iter, ] <- beta.time.mn

  for (p in 1:2) {
    this.update <- updateGPTau(SS.int = SS.int[p],
                               SS.time = SS.time[p],
                               tau.a = 0.1, tau.b = 0.1, ns = ns)
    tau.int[p]  <- this.update$tau.int
    tau.time[p] <- this.update$tau.time
  }

  keep.tau.int[iter, ]  <- tau.int
  keep.tau.time[iter, ] <- tau.time


  if (iter %% 1000 == 0) {
    if (iter < burn) {
      start <- max(1, iter - 3000)
    } else {
      start <- burn + 1
    }
    params <- c("mu", "ls")
    par(mfrow = c(5, 4))
    acc.rate.int  <- round(acc.beta.int / att.beta.int, 3)
    acc.rate.time <- round(acc.beta.time / att.beta.time, 3)
    acc.rate.xi   <- round(acc.xi / att.xi, 3)
    for (i in 1:2) {
      for (p in 1:2) {
        plot(keep.beta.int[start:iter, i, p], type = "l",
             main = paste(params[p], " beta int ", i, " : ",
                          round(beta.int.t[i, p], 3), sep = ""),
             ylab = paste("MH: ", round(MH.beta.int[i, p], 3)),
             xlab = acc.rate.int[i, p])
        plot(keep.beta.time[start:iter, i, p], type = "l",
             main = paste(params[p], " beta time ", i, " : ",
                          round(beta.time.t[i, p], 3), sep = ""),
             ylab = paste("MH: ", round(MH.beta.time[i, p], 3)),
             xlab = acc.rate.time[i, p])
      }
    }
    for (i in 1) {
      for (j in 1:2) {
        plot(keep.mu[start:iter, i, j], type = "l",
             main = paste("mu ", i, ", ", j, ": ", round(mu.t[i, j], 3),
                          sep = ""))
      }
      for (j in 1) {
        plot(keep.ls[start:iter, i, j], type = "l",
             main = paste("ls ", i, ", ", j, ": ", round(ls.t[i, j], 3),
                          sep = ""))
      }
      plot(keep.xi[start:iter], type = "l",
           main = paste("xi: ", xi.t, sep = ""),
           ylab = paste("MH: ", round(MH.xi, 3)),
           xlab = acc.rate.xi)
    }

    for (p in 1:2) {
      plot(keep.tau.int[start:iter, p], type = "l",
           main = paste(params[p], " tau int ", round(tau.int.t[p], 3),
                        sep = ""))
      plot(keep.tau.time[start:iter, p], type = "l",
           main = paste(params[p], " tau time ", round(tau.time.t[p], 3),
                        sep = ""))
    }
    for (p in 1:2) {
      plot(keep.beta.int.mn[start:iter, p], type = "l",
           main = paste(params[p], " beta int mean ",
                        round(beta.int.mn.t[p], 3), sep = ""))
      plot(keep.beta.time.mn[start:iter, p], type = "l",
           main = paste(params[p], " beta time mean ",
                        round(beta.time.mn.t[p], 3), sep = ""))
    }
  }
}

#### testing bandwidth update ####
rm(list=ls())
source("../../../usefulR/usefulfunctions.R", chdir = TRUE)
openblas.set.num.threads(3)
source("auxfunctions.R", chdir = TRUE)
source("updatemodel.R", chdir = TRUE)

set.seed(2000)
ns <- 200
nt <- 10

bw.t <- 0.2
s <- cbind(runif(ns), runif(ns))
d <- rdist(s)
Sigma.t <- exp(-d / bw.t)
Qb.t    <- chol2inv(chol(Sigma.t))

tau.int.t  <- c(0.25, 0.25)
tau.time.t <- c(0.25, 0.25)
beta.int.mn.t  <- rep(0, 2)
beta.time.mn.t <- rep(0, 2)
beta.int.t <- beta.time.t <- matrix(0, ns, 2)
for (p in 1:2) {
  beta.int.t[, p] <- mvrnorm(1, mu = rep(beta.int.mn.t[p], ns),
                             Sigma = Sigma.t / tau.int.t[p])
  beta.time.t[, p] <- mvrnorm(1, mu = rep(beta.time.mn.t[p], ns),
                              Sigma = Sigma.t / tau.time.t[p])
}

mu.t <- ls.t <- matrix(0, ns, nt)
time <- (1:nt - nt / 2) / nt
for (t in 1:nt) {
  mu.t[, t] <- beta.int.t[, 1] + beta.time.t[, 1] * time[t]
  ls.t[, t] <- beta.int.t[, 2] + beta.time.t[, 2] * time[t]
}

SS.int.t <- SS.time.t <- rep(0, 2)
for (p in 1:2) {
  SS.int.t[p]  <- quad.form(Qb.t, beta.int.t[, p] - beta.int.mn.t[p])
  SS.time.t[p] <- quad.form(Qb.t, beta.time.t[, p] - beta.time.mn.t[p])
}

xi.t <- -0.3
y.t <- matrix(0, ns, nt)
for (t in 1:nt) {
  y.t[, t] <- rgev(n = ns, loc = mu.t[, t], exp(ls.t[, t]), xi.t)
}

thresh.t <- matrix(-Inf, ns, nt)
theta.t <- matrix(1, ns, nt)
theta.xi.t <- theta.t^xi.t
alpha.t <- 1

# initialize values
bw <- 0.5
Sigma <- exp(-d / bw)
Qb    <- chol2inv(chol(Sigma))
logdetQb <- logdet(Qb)

SS.int <- SS.time <- rep(0, 2)
for (i in 1:2) {
  SS.int[i]  <- quad.form(Qb, beta.int.t[, i] - beta.int.mn.t[i])
  SS.time[i] <- quad.form(Qb, beta.time.t[, i] - beta.time.mn.t[i])
}

niters  <- 15000
burn    <- 10000
keep.bw <- rep(0, niters)
acc.bw  <- att.bw <- MH.bw <- 0.5

set.seed(3366)  # demo
Rprof(filename = "Rprof.out", line.profiling = TRUE)
for (iter in 1:niters) {
  this.update <- updateBW(bw = bw, bw.min = 0.001, bw.max = max(d), Qb = Qb,
                          logdetQb = logdetQb, d = d, beta.int = beta.int.t,
                          beta.int.mn = beta.int.mn.t, tau.int = tau.int.t,
                          SS.int = SS.int, beta.time = beta.time.t,
                          tau.time = tau.time.t, SS.time = SS.time,
                          beta.time.mn = beta.time.mn.t,
                          acc = acc.bw, att = att.bw, MH = MH.bw)

  bw       <- this.update$bw
  Qb       <- this.update$Qb
  logdetQb <- this.update$logdetQb
  SS.int   <- this.update$SS.int
  SS.time  <- this.update$SS.time
  acc.bw   <- this.update$acc
  att.bw   <- this.update$att

  keep.bw[iter] <- bw

  if (iter < burn / 2) {
    this.update <- mhUpdate(acc = acc.bw, att = att.bw, MH = MH.bw,
                            nattempts = 50)
    acc.bw <- this.update$acc
    att.bw <- this.update$att
    MH.bw  <- this.update$MH
  }

  if (iter %% 1000 == 0) {
    if (iter <= burn) {
      start <- max(1, iter - 3000)
    } else {
      start <- burn + 1
    }
    acc.rate <- round(acc.bw / att.bw, 3)
    plot(keep.bw[start:iter], type = "l", main = "bandwidth",
         ylab = paste("MH: ", round(MH.bw, 3)),
         xlab = acc.rate)
  }
}
Rprof(filename = NULL)
summaryRprof(filename = "Rprof.out", lines = "show")