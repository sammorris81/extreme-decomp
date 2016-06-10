rm(list=ls())
library(fields)
library(Rcpp)
library(emulator)
library(microbenchmark)
library(SpatialExtremes)
library(numDeriv)
library(fields)

#### testing beta ####
source("../../../usefulR/usefulfunctions.R", chdir = TRUE)
openblas.set.num.threads(3)
source("auxfunctions.R", chdir = TRUE)
source("updatemodel.R", chdir = TRUE)

set.seed(2000)
ns <- 400
nt <- 2
np <- 6
X <- rX(ns, nt, np)
beta.t <- rnorm(np, 0, 10)

phi <- 0.2
s <- cbind(runif(ns), runif(ns))
d <- rdist(s)
Sigma <- exp(-d / phi)
tau   <- rgamma(nt, 0.5, 0.5)
Qb <- chol2inv(chol(Sigma))
Xb <- getXBeta(X = X, beta = beta.t)

if (nt == 1) {
  Xb <- matrix(Xb, ns, nt)
}

mu <- matrix(0, ns, nt)
for (t in 1:nt) {
  mu[, t] <- Xb[, t] + t(chol(Sigma)) %*% rnorm(ns) / sqrt(tau[t])
}

# initialize values
SS <- diag(quad.form(Qb, mu - Xb))

niters <- 10000
beta.keep <- matrix(0, niters, np)
beta <- rep(0, np)
set.seed(3366)  # demo
for (iter in 1:niters) {
  this.update <- updateGPBeta(beta.sd = 100, Qb = Qb,
                              param = mu, X = X, SS = SS, tau = tau)
  beta <- this.update$beta
  Xb   <- this.update$Xb
  SS   <- this.update$SS
  beta.keep[iter, ] <- beta

  if (iter %% 500 == 0) {
    start <- max(1, iter - 2000)
    par(mfrow = c(2, np / 2))
    for (i in 1:6) {
      plot(beta.keep[start:iter, i], type = "l",
           main = paste("Beta = ", round(beta.t[i], 3)))
    }
  }
}

#### testing phi ####
rm(list=ls())
source("../../../usefulR/usefulfunctions.R", chdir = TRUE)
openblas.set.num.threads(3)
source("auxfunctions.R", chdir = TRUE)
source("updatemodel.R", chdir = TRUE)

set.seed(2000)
ns <- 400
nt <- 10
np <- 6
X <- rX(ns, nt, np)
beta.t <- rnorm(np, 0, 10)

phi.t <- 0.2
s <- cbind(runif(ns), runif(ns))
d <- rdist(s)
Sigma.t <- exp(-d / phi.t)
tau   <- rgamma(nt, 0.5, 0.5)
Qb <- chol2inv(chol(Sigma.t))
Xb <- getXBeta(X = X, beta = beta.t)

if (nt == 1) {
  Xb <- matrix(Xb, ns, nt)
}

mu <- matrix(0, ns, nt)
for (t in 1:nt) {
  mu[, t] <- Xb[, t] + t(chol(Sigma.t)) %*% rnorm(ns) / sqrt(tau[t])
}

# initialize values
phi <- 0.05
Qb <- chol2inv(chol(exp(-d / phi)))
SS <- diag(quad.form(Qb, mu - Xb))
phi <- 0.05

niters <- 10000
phi.keep <- rep(0, niters)
acc.phi <- att.phi <- MH.phi <- 0.1
set.seed(3366)  # demo
for (iter in 1:niters) {
  this.update <- updateGPBW(bw = phi, bw.min = 0.01, bw.max = 1.2,
                            bw.mn = 0, bw.sd = 1, Qb = Qb, d = d,
                            mu = mu, Xb1 = Xb, tau1 = tau, SS1 = SS,
                            ls = mu, Xb2 = Xb, tau2 = tau, SS2 = SS,
                            acc = acc.phi, att = att.phi, MH = MH.phi)
  phi <- this.update$bw
  Qb <- this.update$Qb
  SS <- this.update$SS1
  acc.phi <- this.update$acc
  att.phi <- this.update$att

  this.update <- mhUpdate(acc = acc.phi, att = att.phi, MH = MH.phi)
  acc.phi <- this.update$acc
  att.phi <- this.update$att
  MH.phi  <- this.update$MH

  phi.keep[iter] <- phi
  if (iter %% 500 == 0) {
    start <- max(1, iter - 2000)
    plot(phi.keep[start:iter], type = "l")
  }
}

#### testing tau ####
rm(list=ls())
source("../../../usefulR/usefulfunctions.R", chdir = TRUE)
openblas.set.num.threads(3)
source("auxfunctions.R", chdir = TRUE)
source("updatemodel.R", chdir = TRUE)

set.seed(2000)
ns <- 400
nt <- 12
np <- 6
X <- rX(ns, nt, np)
beta.t <- rnorm(np, 0, 10)

phi.t <- 0.2
s <- cbind(runif(ns), runif(ns))
d <- rdist(s)
Sigma.t <- exp(-d / phi.t)
tau.t   <- rgamma(nt, 0.5, 0.5)
Qb <- chol2inv(chol(Sigma.t))
Xb <- getXBeta(X = X, beta = beta.t)

if (nt == 1) {
  Xb <- matrix(Xb, ns, nt)
}

mu <- matrix(0, ns, nt)
for (t in 1:nt) {
  mu[, t] <- Xb[, t] + t(chol(Sigma.t)) %*% rnorm(ns) / sqrt(tau.t[t])
}

SS <- getGPSS(Qb = Qb, param = mu, Xb = Xb)

# initialize values
niters <- 10000
tau <- rep(1, nt)
tau.keep <- matrix(0, niters, nt)
set.seed(3366)  # demo
for (iter in 1:niters) {
  this.update <- updateGPTau(SS = SS, tau.a = 0.1, tau.b = 0.1,
                             ns = ns)
  tau <- this.update$tau

  tau.keep[iter, ] <- tau
  if (iter %% 500 == 0) {
    start <- max(1, iter - 2000)
    par(mfrow = c(4, 3))
    for (t in 1:nt) {
      plot(tau.keep[start:iter, t], type = "l",
           main = paste("tau = ", round(tau.t[t], 3)))
    }
  }
}

#### testing tau, phi, and beta ####
rm(list=ls())
source("../../../usefulR/usefulfunctions.R", chdir = TRUE)
openblas.set.num.threads(3)
source("auxfunctions.R", chdir = TRUE)
source("updatemodel.R", chdir = TRUE)

set.seed(2000)
ns <- 400
nt <- 12
np <- 6
X <- rX(ns, nt, np)
beta.t <- rnorm(np, 0, 10)

phi.t <- 0.2
s <- cbind(runif(ns), runif(ns))
d <- rdist(s)
Sigma.t <- exp(-d / phi.t)
tau.t   <- rgamma(nt, 1, 1)
Qb.t    <- chol2inv(chol(Sigma.t))
Xb.t    <- getXBeta(X = X, beta = beta.t)


mu <- matrix(0, ns, nt)
for (t in 1:nt) {
  mu[, t] <- Xb.t[, t] + t(chol(Sigma.t)) %*% rnorm(ns) / sqrt(tau.t[t])
}

# initialize values
beta <- rep(0, np)
Xb   <- getXBeta(X = X, beta = beta)
tau <- rep(1, nt)
phi <- 0.05
Qb <- chol2inv(chol(exp(-d / phi)))
SS <- getGPSS(Qb = Qb, param = mu, Xb = Xb)

niters <- 2000
burn   <- 1500
beta.sd <- 100
beta <- rep(0, np)
beta.keep <- matrix(0, niters, np)
beta.sd.keep <- rep(0, niters)
tau <- rep(1, nt)
tau.keep <- matrix(0, niters, nt)
phi.keep <- rep(0, niters)
acc.phi <- att.phi <- MH.phi <- 0.1
set.seed(3366)  # demo
for (iter in 1:niters) {
  this.update <- updateGPBeta(beta.sd = beta.sd, Qb = Qb,
                              param = mu, X = X, SS = SS, tau = tau)
  beta <- this.update$beta
  Xb   <- this.update$Xb
  SS   <- this.update$SS
  beta.keep[iter, ] <- beta

  this.update <- updateGPBetaSD(beta = beta, tau.a = 0.1, tau.b = 0.1)
  beta.sd <- this.update$beta.sd

  beta.sd.keep[iter] <- beta.sd

  this.update <- updateGPTau(SS = SS, tau.a = 0.1, tau.b = 0.1,
                             ns = ns)
  tau <- this.update$tau

  tau.keep[iter, ] <- tau

  this.update <- updateGPBW(bw = phi, bw.min = 0.01, bw.max = 1.2,
                            bw.mn = 0, bw.sd = 1, Qb = Qb, d = d,
                            mu = mu, Xb1 = Xb, tau1 = tau, SS1 = SS,
                            ls = mu, Xb2 = Xb, tau2 = tau, SS2 = SS,
                            acc = acc.phi, att = att.phi, MH = MH.phi)
  phi     <- this.update$bw
  Qb      <- this.update$Qb
  SS      <- this.update$SS1
  acc.phi <- this.update$acc
  att.phi <- this.update$att
  phi.keep[iter] <- phi

  if (iter < burn / 2) {
    this.update <- mhUpdate(acc = acc.phi, att = att.phi, MH = MH.phi)
    acc.phi <- this.update$acc
    att.phi <- this.update$att
    MH.phi  <- this.update$MH
  }

  if (iter %% 100 == 0) {
    par(mfrow = c(5, 3))
    for (i in 1:np) {
      plot(beta.keep[1:iter, i], type = "l",
           main = paste("beta: ", round(beta.t[i], 3)))
    }
    plot(beta.sd.keep[1:iter], type = "l", main = "beta sd")
    plot(phi.keep[1:iter], type = "l", main = paste("phi: ", phi.t))
    for(i in 1:7) {
      plot(tau.keep[1:iter, i], type = "l",
           main = paste("tau: ", round(tau.t[i], 3)))
    }
  }
}

#### testing mu ####
rm(list=ls())
source("../../../usefulR/usefulfunctions.R", chdir = TRUE)
openblas.set.num.threads(3)
source("auxfunctions.R", chdir = TRUE)
source("updatemodel.R", chdir = TRUE)

set.seed(2000)
ns <- 400
nt <- 12
np <- 6
X <- rX(ns, nt, np)
beta.t <- rnorm(np, 0, 10)

phi.t <- 0.2
s <- cbind(runif(ns), runif(ns))
d <- rdist(s)
Sigma.t <- exp(-d / phi.t)
tau.t   <- rgamma(nt, 1, 1)
Qb.t    <- chol2inv(chol(Sigma.t))
Xb.t    <- getXBeta(X = X, beta = beta.t)

mu.t <- matrix(0, ns, nt)
for (t in 1:nt) {
  mu.t[, t] <- Xb.t[, t] + t(chol(Sigma.t)) %*% rnorm(ns) / sqrt(tau.t[t])
}

ls.t <- matrix(0, ns, nt)
xi.t <- 0.01
y.t <- rgev(n = ns * nt, loc = mu.t, exp(ls.t), xi.t)

Sigma <- solve(Qb.t * tau.t[t])

# initialize values
mu <- matrix(mu.t + rnorm(ns * nt), ns, nt)
SS <- getGPSS(Qb = Qb.t, param = mu, Xb = Xb.t)
curll <- matrix(0, ns, nt)
for (t in 1:nt) {
  curll[, t] <- dgev(x = y.t[, t], loc = mu[, t], exp(ls.t[, t]), xi.t,
                     log = TRUE)
}

niters <- 10000
burn   <- 8000
mu.keep <- array(0, dim = c(niters, ns, nt))
acc.mu <- att.mu <- MH.mu <- matrix(0.2, ns, nt)

set.seed(3366)  # demo
for (iter in 1:niters) {
  this.update <- updateMuTest(mu = mu, Qb = Qb.t, tau = tau.t, Xb = Xb.t,
                              y = y.t, ls = ls.t, xi = xi.t,
                              SS = SS, curll = curll, acc = acc.mu,
                              att = att.mu, MH = MH.mu)
  mu <- this.update$mu
  SS <- this.update$SS
  curll <- this.update$curll
  acc.mu <- this.update$acc
  att.mu <- this.update$att
  mu.keep[iter, , ] <- mu

  if (iter < burn / 2) {
    this.update <- mhUpdate(acc = acc.mu, att = att.mu, MH = MH.mu,
                            target.min = 0.4, target.max = 0.7,
                            nattempts = 400)
    acc.mu <- this.update$acc
    att.mu <- this.update$att
    MH.mu  <- this.update$MH
  }

  if (iter %% 500 == 0) {
    par(mfrow = c(3, 3))
    start <- max(1, iter - 20000)
    for (i in 1:3) {
      for (j in 1:3) {
        plot(mu.keep[start:iter, i, j], type = "l",
             main = paste("mu: ", round(mu.t[i, j], 3)),
             ylab = round(acc.mu[i, j] / att.mu[i, j], 3),
             xlab = MH.mu[i, j])
      }
    }
  }
}


#### testing mu and tau ####
rm(list=ls())
source("../../../usefulR/usefulfunctions.R", chdir = TRUE)
openblas.set.num.threads(3)
source("auxfunctions.R", chdir = TRUE)
source("updatemodel.R", chdir = TRUE)

set.seed(2000)
ns <- 400
nt <- 12
np <- 6
X <- rX(ns, nt, np)
beta.t <- rnorm(np, 0, 10)

phi.t <- 0.2
s <- cbind(runif(ns), runif(ns))
d <- rdist(s)
Sigma.t <- exp(-d / phi.t)
tau.t   <- rgamma(nt, 1, 1)
Qb.t    <- chol2inv(chol(Sigma.t))
Xb.t    <- getXBeta(X = X, beta = beta.t)

mu.t <- matrix(0, ns, nt)
for (t in 1:nt) {
  mu.t[, t] <- Xb.t[, t] + t(chol(Sigma.t)) %*% rnorm(ns) / sqrt(tau.t[t])
}

ls.t <- matrix(0, ns, nt)
xi.t <- 0.01
y.t <- rgev(n = ns * nt, loc = mu.t, exp(ls.t), xi.t)

# initialize values
mu <- matrix(mu.t + rnorm(ns * nt), ns, nt)
SS <- getGPSS(Qb = Qb.t, param = mu, Xb = Xb.t)
curll <- matrix(0, ns, nt)
for (t in 1:nt) {
  curll[, t] <- dgev(x = y.t[, t], loc = mu[, t], exp(ls.t[, t]), xi.t,
                     log = TRUE)
}

niters <- 60000
burn   <- 50000
mu.keep <- array(0, dim = c(niters, ns, nt))
tau <- rep(1, nt)
tau.keep <- matrix(0, niters, nt)
acc.mu <- att.mu <- MH.mu <- matrix(0.1, ns, nt)

set.seed(3366)  # demo
for (iter in 1:niters) {
  this.update <- updateMuTest(mu = mu, Qb = Qb.t, tau = tau, Xb = Xb.t,
                              y = y.t, ls = ls.t, xi = xi.t,
                              SS = SS, curll = curll, acc = acc.mu,
                              att = att.mu, MH = MH.mu)
  mu <- this.update$mu
  SS <- this.update$SS
  curll <- this.update$curll
  acc.mu <- this.update$acc
  att.mu <- this.update$att
  mu.keep[iter, , ] <- mu

  this.update <- updateGPTau(SS = SS, tau.a = 0.1, tau.b = 0.1,
                             ns = ns)
  tau <- this.update$tau
  tau.keep[iter, ] <- tau

  if (iter < burn / 2) {
    this.update <- mhUpdate(acc = acc.mu, att = att.mu, MH = MH.mu,
                            target.min = 0.5, target.max = 0.7,
                            nattempts = 200)
    acc.mu <- this.update$acc
    att.mu <- this.update$att
    MH.mu  <- this.update$MH
  }

  if (iter %% 1000 == 0) {
    par(mfrow = c(4, 3))
    start <- max(1, iter - 20000)
    for (i in 1:3) {
      for (j in 1:3) {
        plot(mu.keep[start:iter, i, j], type = "l",
             main = paste("mu: ", round(mu.t[i, j], 3)),
             ylab = round(acc.mu[i, j] / att.mu[i, j], 3),
             xlab = MH.mu[i, j])
      }
    }
    for (i in 1:3) {
      plot(tau.keep[start:iter, i], type = "l",
           main = paste("tau: ", round(tau.t[i], 3)))
    }
  }
}

#### testing mu, tau, and beta ####
rm(list=ls())
source("../../../usefulR/usefulfunctions.R", chdir = TRUE)
openblas.set.num.threads(3)
source("auxfunctions.R", chdir = TRUE)
source("updatemodel.R", chdir = TRUE)

set.seed(2000)
ns <- 400
nt <- 12
np <- 6
X <- rX(ns, nt, np)
beta.t <- rnorm(np, 0, 10)

phi.t <- 0.2
s <- cbind(runif(ns), runif(ns))
d <- rdist(s)
Sigma.t <- exp(-d / phi.t)
tau.t   <- rgamma(nt, 1, 1)
Qb.t    <- chol2inv(chol(Sigma.t))
Xb.t    <- getXBeta(X = X, beta = beta.t)

mu.t <- matrix(0, ns, nt)
for (t in 1:nt) {
  mu.t[, t] <- Xb.t[, t] + t(chol(Sigma.t)) %*% rnorm(ns) / sqrt(tau.t[t])
}

ls.t <- matrix(0, ns, nt)
xi.t <- 0.01
y.t <- rgev(n = ns * nt, loc = mu.t, exp(ls.t), xi.t)

# initialize values
mu <- matrix(mu.t + rnorm(ns * nt), ns, nt)
SS <- getGPSS(Qb = Qb.t, param = mu, Xb = Xb.t)
curll <- matrix(0, ns, nt)
for (t in 1:nt) {
  curll[, t] <- dgev(x = y.t[, t], loc = mu[, t], exp(ls.t[, t]), xi.t,
                     log = TRUE)
}

niters <- 20000
burn   <- 15000
beta.sd <- 100
beta <- rep(0, np)
beta.keep <- matrix(0, niters, np)
beta.sd.keep <- rep(0, niters)
mu.keep <- array(0, dim = c(niters, ns, nt))
tau <- rep(1, nt)
tau.keep <- matrix(0, niters, nt)
acc.mu <- att.mu <- MH.mu <- matrix(0.1, ns, nt)

set.seed(3366)  # demo
for (iter in 1:niters) {
  this.update <- updateGPBeta(beta.sd = beta.sd, Qb = Qb.t,
                              param = mu, X = X, SS = SS, tau = tau)
  beta <- this.update$beta
  Xb   <- this.update$Xb
  SS   <- this.update$SS
  beta.keep[iter, ] <- beta

  this.update <- updateGPBetaSD(beta = beta, tau.a = 0.1, tau.b = 1)
  beta.sd <- this.update$beta.sd

  beta.sd.keep[iter] <- beta.sd

  this.update <- updateMuTest(mu = mu, Qb = Qb.t, tau = tau.t, Xb = Xb,
                              y = y.t, ls = ls.t, xi = xi.t,
                              SS = SS, curll = curll, acc = acc.mu,
                              att = att.mu, MH = MH.mu)
  mu <- this.update$mu
  SS <- this.update$SS
  curll <- this.update$curll
  acc.mu <- this.update$acc
  att.mu <- this.update$att
  mu.keep[iter, , ] <- mu

  this.update <- updateGPTau(SS = SS, tau.a = 0.1, tau.b = 0.1,
                             ns = ns)
  tau <- this.update$tau
  tau.keep[iter, ] <- tau

  if (iter < burn / 2) {
    this.update <- mhUpdate(acc = acc.mu, att = att.mu, MH = MH.mu,
                            target.min = 0.5, target.max = 0.7,
                            nattempts = 200)
    acc.mu <- this.update$acc
    att.mu <- this.update$att
    MH.mu  <- this.update$MH
  }

  if (iter %% 1000 == 0) {
    par(mfrow = c(4, 3))
    start <- max(1, iter - 20000)
    for (i in 1:2) {
      for (j in 1:3) {
        plot(mu.keep[start:iter, i, j], type = "l",
             main = paste("mu: ", round(mu.t[i, j], 3)),
             ylab = round(acc.mu[i, j] / att.mu[i, j], 3),
             xlab = MH.mu[i, j])
      }
    }

    for (i in 1:3) {
      plot(beta.keep[start:iter, i], type = "l",
           main = paste("beta: ", round(beta.t[i], 3)))
    }

    for (i in 1:3) {
      plot(tau.keep[start:iter, i], type = "l",
           main = paste("tau: ", round(tau.t[i], 3)))
    }
  }
}

#### Verify gradients - no residual dependence ####
rm(list=ls())
source("../../../usefulR/usefulfunctions.R", chdir = TRUE)
openblas.set.num.threads(3)
source("auxfunctions.R", chdir = TRUE)
source("updatemodel.R", chdir = TRUE)

set.seed(2000)
ns <- 10
nt <- 3
np <- 6
X1 <- rX(ns, nt, np)
X2 <- rX(ns, nt, np)
beta1.t <- rnorm(np, 0, 1)
beta2.t <- rnorm(np, 0, 0.1)

phi.t <- 0.2
s <- cbind(runif(ns), runif(ns))
d <- rdist(s)
Sigma.t <- exp(-d / phi.t)
tau.t   <- rgamma(nt, 1, 1)
Qb.t    <- chol2inv(chol(Sigma.t))

Xb1.t <- getXBeta(X = X1, beta = beta1.t)
Xb2.t <- getXBeta(X = X2, beta = beta2.t)

mu.t <- ls.t <- matrix(0, ns, nt)
for (t in 1:nt) {
  mu.t[, t] <- Xb1.t[, t] + t(chol(Sigma.t)) %*% rnorm(ns) / sqrt(tau.t[t])
  ls.t[, t] <- Xb2.t[, t] + t(chol(Sigma.t)) %*% rnorm(ns) / sqrt(tau.t[t])
}

xi.t <- 0.1
y.t <- rgev(n = ns * nt, loc = mu.t, scale = exp(ls.t), xi.t)

lp.mu <- logpost.mu.test(mu = mu.t[, t], Xb = Xb1.t[, t], tau = tau.t[t],
                         Qb = Qb.t, y = y.t[, t], ls = ls.t[, t], xi = xi.t)

mean(grad(func = logpost.mu.test, x = mu.t[, t], Xb = Xb1.t[, t],
          tau = tau.t[t], Qb = Qb.t, y = y.t[, t], ls = ls.t[, t], xi = xi.t) /
       logpost.mu.grad.test(mu = mu.t[, t], Xb = Xb1.t[, t], tau = tau.t[t],
                            Qb = Qb.t, y = y.t[, t], ls = ls.t[, t],
                            xi = xi.t))

sd(grad(func = logpost.mu.test, x = mu.t[, t], Xb = Xb1.t[, t], tau = tau.t[t],
        Qb = Qb.t, y = y.t[, t], ls = ls.t[, t], xi = xi.t) /
     logpost.mu.grad.test(mu = mu.t[, t], Xb = Xb1.t[, t], tau = tau.t[t],
                          Qb = Qb.t, y = y.t[, t], ls = ls.t[, t],
                          xi = xi.t))

lp.logsig <- logpost.logsig.test(ls = ls.t[, t], Xb = Xb2.t[, t],
                                 tau = tau.t[t], Qb = Qb.t, y = y.t[, t],
                                 mu = mu.t[, t], xi = xi.t)

mean(grad(func = logpost.logsig.test, x = ls.t[, t], Xb = Xb2.t[, t],
          tau = tau.t[t], Qb = Qb.t, y = y.t[, t], mu = mu.t[, t], xi = xi.t) /
       logpost.logsig.grad.test(ls = ls.t[, t], Xb = Xb2.t[, t],
                                tau = tau.t[t], Qb = Qb.t, y = y.t[, t],
                                mu = mu.t[, t], xi = xi.t))

sd(grad(func = logpost.logsig.test, x = ls.t[, t], Xb = Xb2.t[, t],
        tau = tau.t[t], Qb = Qb.t, y = y.t[, t], mu = mu.t[, t], xi = xi.t) /
     logpost.logsig.grad.test(ls = ls.t[, t], Xb = Xb2.t[, t],
                              tau = tau.t[t], Qb = Qb.t, y = y.t[, t],
                              mu = mu.t[, t], xi = xi.t))

#### testing logsig ####
rm(list=ls())
source("../../../usefulR/usefulfunctions.R", chdir = TRUE)
openblas.set.num.threads(3)
source("auxfunctions.R", chdir = TRUE)
source("updatemodel.R", chdir = TRUE)

set.seed(2000)
ns <- 400
nt <- 12
np <- 6
X1.t <- rX(ns, nt, np)
X2.t <- rX(ns, nt, np)
beta1.t <- rnorm(np, 0, 10)
beta2.t <- rnorm(np, 0, 5)

phi.t <- 0.2
s <- cbind(runif(ns), runif(ns))
d <- rdist(s)
Sigma.t <- exp(-d / phi.t)
tau.t   <- rgamma(nt, 1, 1)
Qb.t    <- chol2inv(chol(Sigma.t))
Xb1.t   <- getXBeta(X = X1.t, beta = beta1.t)
Xb2.t   <- getXBeta(X = X2.t, beta = beta2.t)

mu.t <- ls.t <- matrix(0, ns, nt)
for (t in 1:nt) {
  mu.t[, t] <- Xb1.t[, t] + t(chol(Sigma.t)) %*% rnorm(ns) / sqrt(tau.t[t])
  ls.t[, t] <- Xb2.t[, t] + t(chol(Sigma.t)) %*% rnorm(ns) / sqrt(tau.t[t])
}

xi.t <- 0.01
y.t <- rgev(n = ns * nt, loc = mu.t, exp(ls.t), xi.t)

# initialize values
ls <- matrix(ls.t + rnorm(ns * nt, 0, 0.1), ns, nt)
SS <- getGPSS(Qb = Qb.t, param = ls, Xb = Xb2.t)
curll <- matrix(0, ns, nt)
for (t in 1:nt) {
  curll[, t] <- dgev(x = y.t[, t], loc = mu.t[, t], exp(ls[, t]), xi.t,
                     log = TRUE)
}

niters <- 10000
burn   <- 8000
ls.keep <- array(0, dim = c(niters, ns, nt))
acc.ls <- att.ls <- MH.ls <- matrix(0.1, ns, nt)

set.seed(3366)  # demo
for (iter in 1:niters) {
  this.update <- updateLSTest(ls = ls, tau = tau.t, Xb = Xb2.t, SS = SS,
                              y = y.t, mu = mu.t, xi = xi.t,
                              Qb = Qb.t, curll = curll, acc = acc.ls,
                              att = att.ls, MH = MH.ls)
  ls <- this.update$ls
  SS <- this.update$SS
  curll <- this.update$curll
  acc.ls <- this.update$acc
  att.ls <- this.update$att
  ls.keep[iter, , ] <- ls

  if (iter < burn / 2) {
    this.update <- mhUpdate(acc = acc.ls, att = att.ls, MH = MH.ls,
                            target.min = 0.4, target.max = 0.7,
                            nattempts = 400)
    acc.ls <- this.update$acc
    att.ls <- this.update$att
    MH.ls  <- this.update$MH
  }

  if (iter %% 500 == 0) {
    par(mfrow = c(3, 3))
    start <- max(1, iter - 20000)
    for (i in 1:3) {
      for (j in 1:3) {
        plot(ls.keep[start:iter, i, j], type = "l",
             main = paste("logsig: ", round(ls.t[i, j], 3)),
             ylab = round(acc.ls[i, j] / att.ls[i, j], 3),
             xlab = MH.ls[i, j])
      }
    }
  }
}

#### testing logsig, tau, and beta ####
rm(list=ls())
source("../../../usefulR/usefulfunctions.R", chdir = TRUE)
openblas.set.num.threads(3)
source("auxfunctions.R", chdir = TRUE)
source("updatemodel.R", chdir = TRUE)

set.seed(2000)
ns <- 400
nt <- 12
np <- 6
X1.t <- rX(ns, nt, np)
X2.t <- rX(ns, nt, np)
beta1.t <- rnorm(np, 0, 10)
beta2.t <- rnorm(np, 0, 5)

phi.t <- 0.2
s <- cbind(runif(ns), runif(ns))
d <- rdist(s)
Sigma.t <- exp(-d / phi.t)
tau.t   <- rgamma(nt, 1, 1)
Qb.t    <- chol2inv(chol(Sigma.t))
Xb1.t   <- getXBeta(X = X1.t, beta = beta1.t)
Xb2.t   <- getXBeta(X = X2.t, beta = beta2.t)

mu.t <- ls.t <- matrix(0, ns, nt)
for (t in 1:nt) {
  mu.t[, t] <- Xb1.t[, t] + t(chol(Sigma.t)) %*% rnorm(ns) / sqrt(tau.t[t])
  ls.t[, t] <- Xb2.t[, t] + t(chol(Sigma.t)) %*% rnorm(ns) / sqrt(tau.t[t])
}

xi.t <- 0.01
y.t <- rgev(n = ns * nt, loc = mu.t, exp(ls.t), xi.t)

# initialize values
ls <- matrix(ls.t + rnorm(ns * nt, 0, 0.1), ns, nt)
SS <- getGPSS(Qb = Qb.t, param = ls, Xb = Xb2.t)
curll <- matrix(0, ns, nt)
for (t in 1:nt) {
  curll[, t] <- dgev(x = y.t[, t], loc = mu.t[, t], exp(ls[, t]), xi.t,
                     log = TRUE)
}

niters <- 30000
burn   <- 25000
beta.sd <- 100
beta <- rep(0, np)
beta.keep <- matrix(0, niters, np)
beta.sd.keep <- rep(0, niters)
ls.keep <- array(0, dim = c(niters, ns, nt))
tau <- rep(1, nt)
tau.keep <- matrix(0, niters, nt)
acc.ls <- att.ls <- MH.ls <- matrix(0.1, ns, nt)

set.seed(3366)  # demo
for (iter in 1:niters) {
  this.update <- updateGPBeta(beta.sd = beta.sd, Qb = Qb.t,
                              param = ls, X = X2.t, SS = SS, tau = tau)
  beta <- this.update$beta
  Xb2  <- this.update$Xb
  SS   <- this.update$SS
  beta.keep[iter, ] <- beta

  this.update <- updateGPBetaSD(beta = beta, tau.a = 0.1, tau.b = 1)
  beta.sd <- this.update$beta.sd

  beta.sd.keep[iter] <- beta.sd

  this.update <- updateLSTest(ls = ls, tau = tau.t, Xb = Xb2, SS = SS,
                              y = y.t, mu = mu.t, xi = xi.t,
                              Qb = Qb.t, curll = curll,
                              acc = acc.ls, att = att.ls, MH = MH.ls)
  ls    <- this.update$ls
  SS    <- this.update$SS
  curll <- this.update$curll
  acc.ls <- this.update$acc
  att.ls <- this.update$att
  ls.keep[iter, , ] <- ls

  this.update <- updateGPTau(SS = SS, tau.a = 0.1, tau.b = 0.1,
                             ns = ns)
  tau <- this.update$tau
  tau.keep[iter, ] <- tau

  if (iter < burn / 2) {
    this.update <- mhUpdate(acc = acc.ls, att = att.ls, MH = MH.ls,
                            target.min = 0.4, target.max = 0.7,
                            nattempts = 400)
    acc.ls <- this.update$acc
    att.ls <- this.update$att
    MH.ls  <- this.update$MH
  }

  if (iter %% 1000 == 0) {
    par(mfrow = c(4, 3))
    start <- max(1, iter - 20000)
    for (i in 1:2) {
      for (j in 1:3) {
        plot(ls.keep[start:iter, i, j], type = "l",
             main = paste("logsig: ", round(ls.t[i, j], 3)),
             ylab = round(acc.ls[i, j] / att.ls[i, j], 3),
             xlab = MH.ls[i, j])
      }
    }

    for (i in 1:3) {
      plot(beta.keep[start:iter, i], type = "l",
           main = paste("beta: ", round(beta2.t[i], 3)))
    }

    for (i in 1:3) {
      plot(tau.keep[start:iter, i], type = "l",
           main = paste("tau: ", round(tau.t[i], 3)))
    }
  }
}

#### testing basis bandwidth update ####
rm(list=ls())
source("../../../usefulR/usefulfunctions.R", chdir = TRUE)
openblas.set.num.threads(3)
source("auxfunctions.R", chdir = TRUE)
source("updatemodel.R", chdir = TRUE)

set.seed(2000)
# setting np later after X is created
ns <- 400
nt <- 12
nknots <- 5
time.interact <- TRUE

s <- cbind(runif(ns), runif(ns))
knots <- as.matrix(cover.design(R = s, nd = nknots)$design)
d <- rdist(s)
dw2 <- rdist(s, knots)^2
dw2[dw2 < 1e-4] <- 0

# create the matrix of covariates
X1.t <- X2.t <- array(1, dim = c(ns, nt, 2))
for (t in 1:nt) {
  time <- (t - nt / 2) / nt
  X1.t[, t, 2] <- X2.t[, t, 2] <- time
}

bw.basis.t <- 0.2
B.t <- makeW(dw2 = dw2, rho = bw.basis.t)
X1.t <- add.basis.X(X1.t, B.t, time.interact = time.interact)
X2.t <- add.basis.X(X2.t, B.t, time.interact = time.interact)

np <- dim(X1.t)[3]

beta1.t <- rnorm(np, 0, 10)
beta2.t <- rnorm(np, 0, 5)

bw.gp.t <- 0.2
Sigma.t <- exp(-d / bw.gp.t)
tau1.t <- tau2.t <- rgamma(nt, 1, 1)
Qb.t    <- chol2inv(chol(Sigma.t))

Xb1.t <- getXBeta(X = X1.t, beta = beta1.t)
Xb2.t <- getXBeta(X = X2.t, beta = beta2.t)

mu.t <- ls.t <- matrix(0, ns, nt)
for (t in 1:nt) {
  mu.t[, t] <- Xb1.t[, t] + t(chol(Sigma.t)) %*% rnorm(ns) / sqrt(tau1.t[t])
  ls.t[, t] <- Xb2.t[, t] + t(chol(Sigma.t)) %*% rnorm(ns) / sqrt(tau2.t[t])
}

xi.t <- 0.01
y.t <- rgev(n = ns * nt, loc = mu.t, scale = exp(ls.t), shape = xi.t)

# initialize values
bw.basis <- 0.4
bw.basis.min <- quantile(dw2, 0.01)
bw.basis.max <- quantile(dw2, 0.99)
B <- makeW(dw2 = dw2, rho = bw.basis)
X1 <- rep.basis.X(X = X1.t, newB = B, time.interact = time.interact)
X2 <- rep.basis.X(X = X2.t, newB = B, time.interact = time.interact)
Xb1 <- getXBeta(X = X1, beta = beta1.t)
Xb2 <- getXBeta(X = X2, beta = beta1.t)
SS1 <- getGPSS(Qb = Qb.t, param = mu.t, Xb = Xb1)
SS2 <- getGPSS(Qb = Qb.t, param = ls.t, Xb = Xb2)

niters <- 30000
burn   <- 25000

# storage
bw.basis.keep <- rep(0, niters)
Xb1.keep     <- array(0, dim = c(niters, ns, nt))
acc.bw.basis <- att.bw.basis <- MH.bw.basis <- 0.1

set.seed(3366)  # demo
for (iter in 1:niters) {
  this.update <- updateXBasisBW(bw = bw.basis, bw.min = bw.basis.min,
                                bw.max = bw.basis.max,
                                X1 = X1, beta1 = beta1.t, Xb1 = Xb1,
                                mu = mu.t, tau1 = tau1.t, SS1 = SS1,
                                X2 = X2, beta2 = beta2.t, Xb2 = Xb2,
                                ls = ls.t, tau2 = tau2.t, SS2 = SS2,
                                Qb = Qb.t, dw2 = dw2,
                                time.interact = time.interact,
                                acc = acc.bw.basis, att = att.bw.basis,
                                MH = MH.bw.basis)
  bw.basis <- this.update$bw
  X1  <- this.update$X1
  Xb1 <- this.update$Xb1
  SS1 <- this.update$SS1
  X2  <- this.update$X2
  Xb2 <- this.update$Xb2
  SS2 <- this.update$SS2
  acc.bw.basis <- this.update$acc
  att.bw.basis <- this.update$att

  bw.basis.keep[iter] <- bw.basis
  Xb1.keep[iter, , ] <- Xb1

  if (iter < burn / 2) {
    this.update <- mhUpdate(acc = acc.bw.basis, att = att.bw.basis,
                            MH = MH.bw.basis,
                            target.min = 0.3, target.max = 0.6,
                            lower = 0.8, higher = 1.2)
    acc.bw.basis <- this.update$acc
    att.bw.basis <- this.update$att
    MH.bw.basis  <- this.update$MH
  }

  if (iter %% 1000 == 0) {
    par(mfrow = c(2, 5))
    plot(bw.basis.keep[100:iter], type = "l",
         main = paste("BW basis: ", bw.basis.t, sep = ""),
         ylab = round(acc.bw.basis / att.bw.basis, 3),
         xlab = MH.bw.basis)
    for(i in 1:3) { for (j in 1:3) {
      plot(Xb1.keep[100:iter, i, j], type = "l",
           main = paste("Xb1: ", Xb1.t[i, j], sep = ""))
    }}
  }
}

#### testing basis bandwidth, beta1, beta2 ####
rm(list=ls())
source("../../../usefulR/usefulfunctions.R", chdir = TRUE)
openblas.set.num.threads(3)
source("auxfunctions.R", chdir = TRUE)
source("updatemodel.R", chdir = TRUE)

set.seed(2000)
# setting np later after X is created
ns <- 400
nt <- 12
nknots <- 5
time.interact <- TRUE

s <- cbind(runif(ns), runif(ns))
knots <- as.matrix(cover.design(R = s, nd = nknots)$design)
d <- rdist(s)
dw2 <- rdist(s, knots)^2
dw2[dw2 < 1e-4] <- 0

# create the matrix of covariates
X1.t <- X2.t <- array(1, dim = c(ns, nt, 2))
for (t in 1:nt) {
  time <- (t - nt / 2) / nt
  X1.t[, t, 2] <- X2.t[, t, 2] <- time
}

bw.basis.t <- 0.2
B.t <- makeW(dw2 = dw2, rho = bw.basis.t)
X1.t <- add.basis.X(X1.t, B.t, time.interact = time.interact)
X2.t <- add.basis.X(X2.t, B.t, time.interact = time.interact)

np <- dim(X1.t)[3]

beta1.t <- rnorm(np, 0, 10)
beta2.t <- rnorm(np, 0, 5)

bw.gp.t <- 0.2
Sigma.t <- exp(-d / bw.gp.t)
tau1.t <- tau2.t <- rgamma(nt, 1, 1)
Qb.t    <- chol2inv(chol(Sigma.t))

Xb1.t <- getXBeta(X = X1.t, beta = beta1.t)
Xb2.t <- getXBeta(X = X2.t, beta = beta2.t)

mu.t <- ls.t <- matrix(0, ns, nt)
for (t in 1:nt) {
  mu.t[, t] <- Xb1.t[, t] + t(chol(Sigma.t)) %*% rnorm(ns) / sqrt(tau1.t[t])
  ls.t[, t] <- Xb2.t[, t] + t(chol(Sigma.t)) %*% rnorm(ns) / sqrt(tau2.t[t])
}

xi.t <- 0.01
y.t <- rgev(n = ns * nt, loc = mu.t, scale = exp(ls.t), shape = xi.t)

# initialize values
beta1 <- beta2 <- rep(0, np)
beta1.sd <- beta2.sd <- 100
bw.basis <- 0.4
bw.basis.min <- quantile(dw2, 0.01)
bw.basis.max <- quantile(dw2, 0.99)
B <- makeW(dw2 = dw2, rho = bw.basis)
X1 <- rep.basis.X(X = X1.t, newB = B, time.interact = time.interact)
X2 <- rep.basis.X(X = X2.t, newB = B, time.interact = time.interact)
Xb1 <- getXBeta(X = X1, beta = beta1)
Xb2 <- getXBeta(X = X2, beta = beta1)
SS1 <- getGPSS(Qb = Qb.t, param = mu.t, Xb = Xb1)
SS2 <- getGPSS(Qb = Qb.t, param = ls.t, Xb = Xb2)

niters <- 30000
burn   <- 25000

# storage
bw.basis.keep <- rep(0, niters)
beta1.keep <- beta2.keep <- matrix(0, niters, np)
beta1.sd.keep <- beta2.sd.keep <- rep(0, niters)
acc.bw.basis <- att.bw.basis <- MH.bw.basis <- 0.1

set.seed(3366)  # demo
for (iter in 1:niters) {
  this.update <- updateGPBeta(beta.sd = beta1.sd, Qb = Qb.t,
                              param = mu.t, X = X1, SS = SS1, tau = tau1.t)
  beta1 <- this.update$beta
  Xb1   <- this.update$Xb
  SS1   <- this.update$SS
  beta1.keep[iter, ] <- beta1

  this.update <- updateGPBetaSD(beta = beta1, tau.a = 0.5, tau.b = 0.5)
  beta1.sd <- this.update$beta.sd
  beta1.sd.keep[iter] <- beta1.sd

  this.update <- updateGPBeta(beta.sd = beta2.sd, Qb = Qb.t,
                              param = ls.t, X = X2, SS = SS2, tau = tau2.t)
  beta2 <- this.update$beta
  Xb2   <- this.update$Xb
  SS2   <- this.update$SS
  beta2.keep[iter, ] <- beta2

  this.update <- updateGPBetaSD(beta = beta2, tau.a = 0.5, tau.b = 0.5)
  beta2.sd <- this.update$beta.sd
  beta2.sd.keep[iter] <- beta2.sd

  this.update <- updateXBasisBW(bw = bw.basis, bw.min = bw.basis.min,
                                bw.max = bw.basis.max,
                                X1 = X1, beta1 = beta1, Xb1 = Xb1,
                                mu = mu.t, tau1 = tau1.t, SS1 = SS1,
                                X2 = X2, beta2 = beta2, Xb2 = Xb2,
                                ls = ls.t, tau2 = tau2.t, SS2 = SS2,
                                Qb = Qb.t, dw2 = dw2,
                                time.interact = time.interact,
                                acc = acc.bw.basis, att = att.bw.basis,
                                MH = MH.bw.basis)
  bw.basis <- this.update$bw
  X1  <- this.update$X1
  Xb1 <- this.update$Xb1
  SS1 <- this.update$SS1
  X2  <- this.update$X2
  Xb2 <- this.update$Xb2
  SS2 <- this.update$SS2
  acc.bw.basis <- this.update$acc
  att.bw.basis <- this.update$att

  bw.basis.keep[iter] <- bw.basis

  if (iter < burn / 2) {
    this.update <- mhUpdate(acc = acc.bw.basis, att = att.bw.basis,
                            MH = MH.bw.basis,
                            target.min = 0.3, target.max = 0.6,
                            lower = 0.8, higher = 1.2)
    acc.bw.basis <- this.update$acc
    att.bw.basis <- this.update$att
    MH.bw.basis  <- this.update$MH
  }

  if (iter %% 500 == 0) {
    par(mfrow = c(3, 3))
    plot(bw.basis.keep[100:iter], type = "l",
         main = paste("BW basis: ", bw.basis.t, sep = ""),
         ylab = round(acc.bw.basis / att.bw.basis, 3),
         xlab = MH.bw.basis)
    for(i in 1:4) {
      plot(beta1.keep[100:iter, i], type = "l",
           main = paste("beta 1: ", beta1.t[i], sep = ""))
    }
    for (i in 1:4) {
      plot(beta2.keep[100:iter, i], type = "l",
           main = paste("beta 2: ", beta2.t[i], sep = ""))
    }

  }
}

#### testing basis bandwidth, beta1, beta2, mu, and logsig ####
rm(list=ls())
source("../../../usefulR/usefulfunctions.R", chdir = TRUE)
openblas.set.num.threads(3)
source("auxfunctions.R", chdir = TRUE)
source("updatemodel.R", chdir = TRUE)

set.seed(2000)
# setting np later after X is created
ns <- 400
nt <- 12
nknots <- 5
time.interact <- TRUE

s <- cbind(runif(ns), runif(ns))
knots <- as.matrix(cover.design(R = s, nd = nknots)$design)
d <- rdist(s)
dw2 <- rdist(s, knots)^2
dw2[dw2 < 1e-4] <- 0

# create the matrix of covariates
X1.t <- X2.t <- array(1, dim = c(ns, nt, 2))
for (t in 1:nt) {
  time <- (t - nt / 2) / nt
  X1.t[, t, 2] <- X2.t[, t, 2] <- time
}

bw.basis.t <- 0.2
B.t <- makeW(dw2 = dw2, rho = bw.basis.t)
X1.t <- add.basis.X(X1.t, B.t, time.interact = time.interact)
X2.t <- add.basis.X(X2.t, B.t, time.interact = time.interact)

np <- dim(X1.t)[3]

beta1.t <- rnorm(np, 0, 10)
beta2.t <- rnorm(np, 0, 5)

bw.gp.t <- 0.2
Sigma.t <- exp(-d / bw.gp.t)
tau1.t <- tau2.t <- rgamma(nt, 1, 1)
Qb.t    <- chol2inv(chol(Sigma.t))

Xb1.t <- getXBeta(X = X1.t, beta = beta1.t)
Xb2.t <- getXBeta(X = X2.t, beta = beta2.t)

mu.t <- ls.t <- matrix(0, ns, nt)
for (t in 1:nt) {
  mu.t[, t] <- Xb1.t[, t] + t(chol(Sigma.t)) %*% rnorm(ns) / sqrt(tau1.t[t])
  ls.t[, t] <- Xb2.t[, t] + t(chol(Sigma.t)) %*% rnorm(ns) / sqrt(tau2.t[t])
}

xi.t <- 0.01
y.t <- rgev(n = ns * nt, loc = mu.t, scale = exp(ls.t), shape = xi.t)

# initialize values
mu <- matrix(mu.t + rnorm(ns * nt), ns, nt)
ls <- matrix(ls.t + rnorm(ns * nt), ns, nt)
beta1 <- beta2 <- rep(0, np)
beta1.sd <- beta2.sd <- 100
bw.basis <- 0.4
bw.basis.min <- quantile(dw2, 0.01)
bw.basis.max <- quantile(dw2, 0.99)
B <- makeW(dw2 = dw2, rho = bw.basis)
X1 <- rep.basis.X(X = X1.t, newB = B, time.interact = time.interact)
X2 <- rep.basis.X(X = X2.t, newB = B, time.interact = time.interact)
Xb1 <- getXBeta(X = X1, beta = beta1)
Xb2 <- getXBeta(X = X2, beta = beta1)
SS1 <- getGPSS(Qb = Qb.t, param = mu.t, Xb = Xb1)
SS2 <- getGPSS(Qb = Qb.t, param = ls.t, Xb = Xb2)

curll <- matrix(0, ns, nt)
for (t in 1:nt) {
  curll[, t] <- dgev(x = y.t[, t], loc = mu[, t], exp(ls[, t]), xi.t,
                     log = TRUE)
}

niters <- 30000
burn   <- 25000

# storage
bw.basis.keep <- rep(0, niters)
beta1.keep <- beta2.keep <- matrix(0, niters, np)
beta1.sd.keep <- beta2.sd.keep <- rep(0, niters)
mu.keep <- array(0, dim = c(niters, ns, nt))
ls.keep <- array(0, dim = c(niters, ns, nt))
acc.bw.basis <- att.bw.basis <- MH.bw.basis <- 0.1
acc.mu <- att.mu <- MH.mu <- matrix(0.2, ns, nt)
acc.ls <- att.ls <- MH.ls <- matrix(0.1, ns, nt)

set.seed(3366)  # demo
for (iter in 1:niters) {
  this.update <- updateGPBeta(beta.sd = beta1.sd, Qb = Qb.t,
                              param = mu, X = X1, SS = SS1, tau = tau1.t)
  beta1 <- this.update$beta
  Xb1   <- this.update$Xb
  SS1   <- this.update$SS
  beta1.keep[iter, ] <- beta1

  this.update <- updateGPBetaSD(beta = beta1, tau.a = 0.5, tau.b = 0.5)
  beta1.sd <- this.update$beta.sd
  beta1.sd.keep[iter] <- beta1.sd

  this.update <- updateGPBeta(beta.sd = beta2.sd, Qb = Qb.t,
                              param = ls, X = X2, SS = SS2, tau = tau2.t)
  beta2 <- this.update$beta
  Xb2   <- this.update$Xb
  SS2   <- this.update$SS
  beta2.keep[iter, ] <- beta2

  this.update <- updateGPBetaSD(beta = beta2, tau.a = 0.5, tau.b = 0.5)
  beta2.sd <- this.update$beta.sd
  beta2.sd.keep[iter] <- beta2.sd

  this.update <- updateXBasisBW(bw = bw.basis, bw.min = bw.basis.min,
                                bw.max = bw.basis.max,
                                X1 = X1, beta1 = beta1, Xb1 = Xb1,
                                mu = mu, tau1 = tau1.t, SS1 = SS1,
                                X2 = X2, beta2 = beta2, Xb2 = Xb2,
                                ls = ls, tau2 = tau2.t, SS2 = SS2,
                                Qb = Qb.t, dw2 = dw2,
                                time.interact = time.interact,
                                acc = acc.bw.basis, att = att.bw.basis,
                                MH = MH.bw.basis)
  bw.basis <- this.update$bw
  X1  <- this.update$X1
  Xb1 <- this.update$Xb1
  SS1 <- this.update$SS1
  X2  <- this.update$X2
  Xb2 <- this.update$Xb2
  SS2 <- this.update$SS2
  acc.bw.basis <- this.update$acc
  att.bw.basis <- this.update$att
  bw.basis.keep[iter] <- bw.basis

  this.update <- updateMuTest(mu = mu, tau = tau1.t, Xb = Xb1, SS = SS1,
                              y = y.t, ls = ls, xi = xi.t,
                              Qb = Qb.t, curll = curll, acc = acc.mu,
                              att = att.mu, MH = MH.mu)
  mu    <- this.update$mu
  SS1   <- this.update$SS
  curll <- this.update$curll
  acc.mu <- this.update$acc
  att.mu <- this.update$att
  mu.keep[iter, , ] <- mu

  this.update <- updateLSTest(ls = ls, tau = tau2.t, Xb = Xb2, SS = SS2,
                              y = y.t, mu = mu, xi = xi.t,
                              Qb = Qb.t, curll = curll, acc = acc.ls,
                              att = att.ls, MH = MH.ls)
  ls     <- this.update$ls
  SS2    <- this.update$SS
  curll  <- this.update$curll
  acc.ls <- this.update$acc
  att.ls <- this.update$att
  ls.keep[iter, , ] <- ls

  if (iter < burn / 2) {
    this.update <- mhUpdate(acc = acc.bw.basis, att = att.bw.basis,
                            MH = MH.bw.basis,
                            target.min = 0.3, target.max = 0.6,
                            lower = 0.8, higher = 1.2)
    acc.bw.basis <- this.update$acc
    att.bw.basis <- this.update$att
    MH.bw.basis  <- this.update$MH

    this.update <- mhUpdate(acc = acc.mu, att = att.mu, MH = MH.mu,
                            target.min = 0.4, target.max = 0.7,
                            nattempts = 400)
    acc.mu <- this.update$acc
    att.mu <- this.update$att
    MH.mu  <- this.update$MH

    this.update <- mhUpdate(acc = acc.ls, att = att.ls, MH = MH.ls,
                            target.min = 0.4, target.max = 0.7,
                            nattempts = 400)
    acc.ls <- this.update$acc
    att.ls <- this.update$att
    MH.ls  <- this.update$MH
  }

  if (iter %% 500 == 0) {
    par(mfrow = c(3, 3))
    plot(bw.basis.keep[100:iter], type = "l",
         main = paste("BW basis: ", bw.basis.t, sep = ""),
         ylab = round(acc.bw.basis / att.bw.basis, 3),
         xlab = MH.bw.basis)
    for(i in 1:4) {
      plot(beta1.keep[100:iter, i], type = "l",
           main = paste("beta 1: ", beta1.t[i], sep = ""))
    }
    for (i in 1:4) {
      plot(beta2.keep[100:iter, i], type = "l",
           main = paste("beta 2: ", beta2.t[i], sep = ""))
    }

  }
}


#### Verify gradients - with residual dependence ####
rm(list=ls())
source("../../../usefulR/usefulfunctions.R", chdir = TRUE)
openblas.set.num.threads(3)
source("auxfunctions.R", chdir = TRUE)
source("updatemodel.R", chdir = TRUE)

set.seed(2000)
ns <- 10
nt <- 3
np <- 6
X1 <- rX(ns, nt, np)
X2 <- rX(ns, nt, np)
beta1.t <- rnorm(np, 0, 1)
beta2.t <- rnorm(np, 0, 0.1)

phi.t <- 0.2
s <- cbind(runif(ns), runif(ns))
d <- rdist(s)
Sigma.t <- exp(-d / phi.t)
tau.t   <- rgamma(nt, 1, 1)
Qb.t    <- chol2inv(chol(Sigma.t))

Xb1.t <- getXBeta(X = X1, beta = beta1.t)
Xb2.t <- getXBeta(X = X2, beta = beta2.t)

mu.t <- ls.t <- matrix(0, ns, nt)
for (t in 1:nt) {
  mu.t[, t] <- Xb1.t[, t] + t(chol(Sigma.t)) %*% rnorm(ns) / sqrt(tau.t[t])
  ls.t[, t] <- Xb2.t[, t] + t(chol(Sigma.t)) %*% rnorm(ns) / sqrt(tau.t[t])
}

xi.t <- 0.1
y.t <- rgev(n = ns * nt, loc = mu.t, scale = exp(ls.t), xi.t)

nknots <- 4
theta.t <- matrix(abs(rnorm(ns * nt)), ns, nt)
alpha.t <- 0.4
thresh.t <- matrix(median(y.t), ns, nt)

xi.t <- 0
lp.mu <- logpost.mu(mu = mu.t[, t], Xb = Xb1.t[, t], tau = tau.t[t],
                    Qb = Qb.t, y = y.t[, t], ls = ls.t[, t], xi = xi.t,
                    theta = theta.t[, t], thresh = thresh.t[, t],
                    alpha = alpha.t)

mean(grad(func = logpost.mu, x = mu.t[, t], Xb = Xb1.t[, t], tau = tau.t[t],
          Qb = Qb.t, y = y.t[, t], ls = ls.t[, t], xi = xi.t,
          theta = theta.t[, t], thresh = thresh.t[, t], alpha = alpha.t) /
       logpost.mu.grad(mu = mu.t[, t], Xb = Xb1.t[, t], tau = tau.t[t],
                       Qb = Qb.t, y = y.t[, t], ls = ls.t[, t],
                       xi = xi.t, theta = theta.t[, t], thresh = thresh.t[, t],
                       alpha = alpha.t))

sd(grad(func = logpost.mu, x = mu.t[, t], Xb = Xb1.t[, t], tau = tau.t[t],
        Qb = Qb.t, y = y.t[, t], ls = ls.t[, t], xi = xi.t,
        theta = theta.t[, t], thresh = thresh.t[, t], alpha = alpha.t) /
     logpost.mu.grad(mu = mu.t[, t], Xb = Xb1.t[, t], tau = tau.t[t],
                     Qb = Qb.t, y = y.t[, t], ls = ls.t[, t],
                     xi = xi.t, theta = theta.t[, t], thresh = thresh.t[, t],
                     alpha = alpha.t))

lp.logsig <- logpost.logsig(ls = ls.t[, t], Xb = Xb2.t[, t], tau = tau.t[t],
                            Qb = Qb.t, y = y.t[, t], mu = mu.t[, t],
                            xi = xi.t, theta = theta.t[, t],
                            thresh = thresh.t[, t], alpha = alpha.t)

mean(grad(func = logpost.logsig, x = ls.t[, t], Xb = Xb2.t[, t],
          tau = tau.t[t], Qb = Qb.t, y = y.t[, t], mu = mu.t[, t], xi = xi.t,
          theta = theta.t[, t], thresh = thresh.t[, t], alpha = alpha.t) /
       logpost.logsig.grad(ls = ls.t[, t], Xb = Xb2.t[, t],
                           tau = tau.t[t], Qb = Qb.t, y = y.t[, t],
                           mu = mu.t[, t], xi = xi.t,
                           theta = theta.t[, t], thresh = thresh.t[, t],
                           alpha = alpha.t))

sd(grad(func = logpost.logsig, x = ls.t[, t], Xb = Xb2.t[, t],
        tau = tau.t[t], Qb = Qb.t, y = y.t[, t], mu = mu.t[, t], xi = xi.t,
        theta = theta.t[, t], thresh = thresh.t[, t], alpha = alpha.t) /
     logpost.logsig.grad(ls = ls.t[, t], Xb = Xb2.t[, t],
                         tau = tau.t[t], Qb = Qb.t, y = y.t[, t],
                         mu = mu.t[, t], xi = xi.t,
                         theta = theta.t[, t], thresh = thresh.t[, t],
                         alpha = alpha.t))


#### testing xi ####
rm(list=ls())
source("../../../usefulR/usefulfunctions.R", chdir = TRUE)
openblas.set.num.threads(3)
source("auxfunctions.R", chdir = TRUE)
source("updatemodel.R", chdir = TRUE)

set.seed(2000)
ns <- 400
nt <- 12
np <- 6
X1.t <- rX(ns, nt, np)
X2.t <- rX(ns, nt, np)
beta1.t <- rnorm(np, 0, 10)
beta2.t <- rnorm(np, 0, 5)

phi.t <- 0.2
s <- cbind(runif(ns), runif(ns))
d <- rdist(s)
Sigma.t <- exp(-d / phi.t)
tau.t   <- rgamma(nt, 1, 1)
Qb.t    <- chol2inv(chol(Sigma.t))
Xb1.t   <- getXBeta(X = X1.t, beta = beta1.t)
Xb2.t   <- getXBeta(X = X2.t, beta = beta2.t)

mu.t <- ls.t <- matrix(0, ns, nt)
for (t in 1:nt) {
  mu.t[, t] <- Xb1.t[, t] + t(chol(Sigma.t)) %*% rnorm(ns) / sqrt(tau.t[t])
  ls.t[, t] <- Xb2.t[, t] + t(chol(Sigma.t)) %*% rnorm(ns) / sqrt(tau.t[t])
}

xi.t <- -0.7
y.t <- rgev(n = ns * nt, loc = mu.t, exp(ls.t), xi.t)

# initialize values
xi <- 0.1

thresh.t <- matrix(-Inf, ns, nt)
theta.t <- matrix(1, ns, nt)
alpha.t <- 1
curll <- loglike(y = y.t, mu = mu.t, ls = ls.t, xi = xi,
                 theta = theta.t, thresh = thresh.t, alpha = alpha.t)


niters <- 10000
burn   <- 8000
xi.keep <- rep(0, niters)
acc.xi <- att.xi <- MH.xi <- 0.01


set.seed(3366)  # demo
for (iter in 1:niters) {
  this.update <- updateXi(xi = xi, xi.min = -2, xi.max = 2,
                          xi.mn = 0, xi.sd = 0.5, y = y.t, mu = mu.t,
                          ls = ls.t, curll = curll, theta = theta.t,
                          thresh = thresh.t, alpha = alpha.t, acc = acc.xi,
                          att = att.xi, MH = MH.xi)
  xi <- this.update$xi
  curll <- this.update$curll
  acc.xi <- this.update$acc
  att.xi <- this.update$att
  xi.keep[iter] <- xi

  if (iter < burn / 2) {
    this.update <- mhUpdate(acc = acc.xi, att = att.xi, MH = MH.xi,
                            target.min = 0.3, target.max = 0.6,
                            nattempts = 400)
    acc.xi <- this.update$acc
    att.xi <- this.update$att
    MH.xi  <- this.update$MH
  }

  if (iter %% 500 == 0) {
    start <- max(1, iter - 20000)
    plot(xi.keep[start:iter], type = "l", main = paste("xi: ", xi.t),
         ylab = round(acc.xi / att.xi, 3),
         xlab = MH.xi)
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
acc.beta1  <- att.beta1 <- MH.beta1 <- matrix(0.1, ns, nt)

set.seed(3366)  # demo
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
                            target.min = 0.3, target.max = 0.6,
                            nattempts = 400)
    acc.beta1 <- this.update$acc
    att.beta1 <- this.update$att
    MH.beta1  <- this.update$MH
  }

  if (iter %% 500 == 0) {
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

#### testing beta1 ####
rm(list=ls())
library(compiler)
enableJIT(3)
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
beta1.int.mn <- 0
beta1.time <- beta1.time.t + rnorm(ns * nt)
beta1.time.mn <- 0
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
keep.beta1.int.mn <- keep.beta1.time.mn <- rep(0)
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

rm(list=ls())
library(fields)
library(Rcpp)
library(emulator)
library(microbenchmark)
library(SpatialExtremes)
library(numDeriv)
library(fields)

#### testing beta ####
source("../../../usefulR/usefulfunctions.R", chdir = TRUE)
source("auxfunctions.R", chdir = TRUE)
source("updatemodel.R", chdir = TRUE)

set.seed(2000)
ns <- 400
nt <- 2
np <- 6
X <- rX(ns, nt, np)
beta.t <- rnorm(np, 0, 10)

phi <- 0.2
s <- cbind(runif(ns), runif(ns))
d <- rdist(s)
Sigma <- exp(-d / phi)
tau   <- rgamma(nt, 0.5, 0.5)
Qb <- chol2inv(chol(Sigma))
Xb <- getXBeta(X = X, beta = beta.t)

if (nt == 1) {
  Xb <- matrix(Xb, ns, nt)
}

mu <- matrix(0, ns, nt)
for (t in 1:nt) {
  mu[, t] <- Xb[, t] + t(chol(Sigma)) %*% rnorm(ns) / sqrt(tau[t])
}

# initialize values
SS <- diag(quad.form(Qb, mu - Xb))

niters <- 10000
beta.keep <- matrix(0, niters, np)
beta <- rep(0, np)
set.seed(3366)  # demo
for (iter in 1:niters) {
  this.update <- updateGPBeta(beta.sd = 100, Qb = Qb,
                              param = mu, X = X, SS = SS, tau = tau)
  beta <- this.update$beta
  Xb   <- this.update$Xb
  SS   <- this.update$SS
  beta.keep[iter, ] <- beta

  if (iter %% 500 == 0) {
    start <- max(1, iter - 2000)
    par(mfrow = c(2, np / 2))
    for (i in 1:6) {
      plot(beta.keep[start:iter, i], type = "l",
           main = paste("Beta = ", round(beta.t[i], 3)))
    }
  }
}

#### testing phi ####
rm(list=ls())
source("../../../usefulR/usefulfunctions.R", chdir = TRUE)
source("auxfunctions.R", chdir = TRUE)
source("updatemodel.R", chdir = TRUE)

set.seed(2000)
ns <- 400
nt <- 10
np <- 6
X <- rX(ns, nt, np)
beta.t <- rnorm(np, 0, 10)

phi.t <- 0.2
s <- cbind(runif(ns), runif(ns))
d <- rdist(s)
Sigma.t <- exp(-d / phi.t)
tau   <- rgamma(nt, 0.5, 0.5)
Qb <- chol2inv(chol(Sigma.t))
Xb <- getXBeta(X = X, beta = beta.t)

if (nt == 1) {
  Xb <- matrix(Xb, ns, nt)
}

mu <- matrix(0, ns, nt)
for (t in 1:nt) {
  mu[, t] <- Xb[, t] + t(chol(Sigma.t)) %*% rnorm(ns) / sqrt(tau[t])
}

# initialize values
phi <- 0.05
Qb <- chol2inv(chol(exp(-d / phi)))
SS <- diag(quad.form(Qb, mu - Xb))
phi <- 0.05

niters <- 10000
phi.keep <- rep(0, niters)
acc.phi <- att.phi <- MH.phi <- 0.1
set.seed(3366)  # demo
for (iter in 1:niters) {
  this.update <- updateGPBW(bw = phi, bw.min = 0.01, bw.max = 1.2,
                            bw.mn = 0, bw.sd = 1, Qb = Qb, d = d,
                            mu = mu, Xb1 = Xb, tau1 = tau, SS1 = SS,
                            ls = mu, Xb2 = Xb, tau2 = tau, SS2 = SS,
                            acc = acc.phi, att = att.phi, MH = MH.phi)
  phi <- this.update$bw
  Qb <- this.update$Qb
  SS <- this.update$SS1
  acc.phi <- this.update$acc
  att.phi <- this.update$att

  this.update <- mhUpdate(acc = acc.phi, att = att.phi, MH = MH.phi)
  acc.phi <- this.update$acc
  att.phi <- this.update$att
  MH.phi  <- this.update$MH

  phi.keep[iter] <- phi
  if (iter %% 500 == 0) {
    start <- max(1, iter - 2000)
    plot(phi.keep[start:iter], type = "l")
  }
}

#### testing tau ####
rm(list=ls())
source("../../../usefulR/usefulfunctions.R", chdir = TRUE)
source("auxfunctions.R", chdir = TRUE)
source("updatemodel.R", chdir = TRUE)

set.seed(2000)
ns <- 400
nt <- 12
np <- 6
X <- rX(ns, nt, np)
beta.t <- rnorm(np, 0, 10)

phi.t <- 0.2
s <- cbind(runif(ns), runif(ns))
d <- rdist(s)
Sigma.t <- exp(-d / phi.t)
tau.t   <- rgamma(nt, 0.5, 0.5)
Qb <- chol2inv(chol(Sigma.t))
Xb <- getXBeta(X = X, beta = beta.t)

if (nt == 1) {
  Xb <- matrix(Xb, ns, nt)
}

mu <- matrix(0, ns, nt)
for (t in 1:nt) {
  mu[, t] <- Xb[, t] + t(chol(Sigma.t)) %*% rnorm(ns) / sqrt(tau.t[t])
}

SS <- getGPSS(Qb = Qb, param = mu, Xb = Xb)

# initialize values
niters <- 10000
tau <- rep(1, nt)
tau.keep <- matrix(0, niters, nt)
set.seed(3366)  # demo
for (iter in 1:niters) {
  this.update <- updateGPTau(SS = SS, tau.a = 0.1, tau.b = 0.1,
                             ns = ns)
  tau <- this.update$tau

  tau.keep[iter, ] <- tau
  if (iter %% 500 == 0) {
    start <- max(1, iter - 2000)
    par(mfrow = c(4, 3))
    for (t in 1:nt) {
      plot(tau.keep[start:iter, t], type = "l",
           main = paste("tau = ", round(tau.t[t], 3)))
    }
  }
}

#### testing tau, phi, and beta ####
rm(list=ls())
source("../../../usefulR/usefulfunctions.R", chdir = TRUE)
source("auxfunctions.R", chdir = TRUE)
source("updatemodel.R", chdir = TRUE)

set.seed(2000)
ns <- 400
nt <- 12
np <- 6
X <- rX(ns, nt, np)
beta.t <- rnorm(np, 0, 10)

phi.t <- 0.2
s <- cbind(runif(ns), runif(ns))
d <- rdist(s)
Sigma.t <- exp(-d / phi.t)
tau.t   <- rgamma(nt, 1, 1)
Qb.t    <- chol2inv(chol(Sigma.t))
Xb.t    <- getXBeta(X = X, beta = beta.t)


mu <- matrix(0, ns, nt)
for (t in 1:nt) {
  mu[, t] <- Xb.t[, t] + t(chol(Sigma.t)) %*% rnorm(ns) / sqrt(tau.t[t])
}

# initialize values
beta <- rep(0, np)
Xb   <- getXBeta(X = X, beta = beta)
tau <- rep(1, nt)
phi <- 0.05
Qb <- chol2inv(chol(exp(-d / phi)))
SS <- getGPSS(Qb = Qb, param = mu, Xb = Xb)

niters <- 2000
burn   <- 1500
beta.sd <- 100
beta <- rep(0, np)
beta.keep <- matrix(0, niters, np)
beta.sd.keep <- rep(0, niters)
tau <- rep(1, nt)
tau.keep <- matrix(0, niters, nt)
phi.keep <- rep(0, niters)
acc.phi <- att.phi <- MH.phi <- 0.1
set.seed(3366)  # demo
for (iter in 1:niters) {
  this.update <- updateGPBeta(beta.sd = beta.sd, Qb = Qb,
                              param = mu, X = X, SS = SS, tau = tau)
  beta <- this.update$beta
  Xb   <- this.update$Xb
  SS   <- this.update$SS
  beta.keep[iter, ] <- beta

  this.update <- updateGPBetaSD(beta = beta, tau.a = 0.1, tau.b = 0.1)
  beta.sd <- this.update$beta.sd

  beta.sd.keep[iter] <- beta.sd

  this.update <- updateGPTau(SS = SS, tau.a = 0.1, tau.b = 0.1,
                             ns = ns)
  tau <- this.update$tau

  tau.keep[iter, ] <- tau

  this.update <- updateGPBW(bw = phi, bw.min = 0.01, bw.max = 1.2,
                            bw.mn = 0, bw.sd = 1, Qb = Qb, d = d,
                            mu = mu, Xb1 = Xb, tau1 = tau, SS1 = SS,
                            ls = mu, Xb2 = Xb, tau2 = tau, SS2 = SS,
                            acc = acc.phi, att = att.phi, MH = MH.phi)
  phi     <- this.update$bw
  Qb      <- this.update$Qb
  SS      <- this.update$SS1
  acc.phi <- this.update$acc
  att.phi <- this.update$att
  phi.keep[iter] <- phi

  if (iter < burn / 2) {
    this.update <- mhUpdate(acc = acc.phi, att = att.phi, MH = MH.phi)
    acc.phi <- this.update$acc
    att.phi <- this.update$att
    MH.phi  <- this.update$MH
  }

  if (iter %% 100 == 0) {
    par(mfrow = c(5, 3))
    for (i in 1:np) {
      plot(beta.keep[1:iter, i], type = "l",
           main = paste("beta: ", round(beta.t[i], 3)))
    }
    plot(beta.sd.keep[1:iter], type = "l", main = "beta sd")
    plot(phi.keep[1:iter], type = "l", main = paste("phi: ", phi.t))
    for(i in 1:7) {
      plot(tau.keep[1:iter, i], type = "l",
           main = paste("tau: ", round(tau.t[i], 3)))
    }
  }
}

#### testing mu ####
rm(list=ls())
source("../../../usefulR/usefulfunctions.R", chdir = TRUE)
source("auxfunctions.R", chdir = TRUE)
source("updatemodel.R", chdir = TRUE)

set.seed(2000)
ns <- 400
nt <- 12
np <- 6
X <- rX(ns, nt, np)
beta.t <- rnorm(np, 0, 10)

phi.t <- 0.2
s <- cbind(runif(ns), runif(ns))
d <- rdist(s)
Sigma.t <- exp(-d / phi.t)
tau.t   <- rgamma(nt, 1, 1)
Qb.t    <- chol2inv(chol(Sigma.t))
Xb.t    <- getXBeta(X = X, beta = beta.t)

mu.t <- matrix(0, ns, nt)
for (t in 1:nt) {
  mu.t[, t] <- Xb.t[, t] + t(chol(Sigma.t)) %*% rnorm(ns) / sqrt(tau.t[t])
}

ls.t <- matrix(0, ns, nt)
xi.t <- 0.01
y.t <- rgev(n = ns * nt, loc = mu.t, exp(ls.t), xi.t)

Sigma <- solve(Qb.t * tau.t[t])

# initialize values
mu <- matrix(mu.t + rnorm(ns * nt), ns, nt)
SS <- getGPSS(Qb = Qb.t, param = mu, Xb = Xb.t)
curll <- matrix(0, ns, nt)
for (t in 1:nt) {
  curll[, t] <- dgev(x = y.t[, t], loc = mu[, t], exp(ls.t[, t]), xi.t,
                     log = TRUE)
}

niters <- 10000
burn   <- 8000
mu.keep <- array(0, dim = c(niters, ns, nt))
acc.mu <- att.mu <- MH.mu <- matrix(0.2, ns, nt)

set.seed(3366)  # demo
for (iter in 1:niters) {
  this.update <- updateMuTest(mu = mu, Qb = Qb.t, tau = tau.t, Xb = Xb.t,
                              y = y.t, ls = ls.t, xi = xi.t,
                              SS = SS, curll = curll, acc = acc.mu,
                              att = att.mu, MH = MH.mu)
  mu <- this.update$mu
  SS <- this.update$SS
  curll <- this.update$curll
  acc.mu <- this.update$acc
  att.mu <- this.update$att
  mu.keep[iter, , ] <- mu

  if (iter < burn / 2) {
    this.update <- mhUpdate(acc = acc.mu, att = att.mu, MH = MH.mu,
                            target.min = 0.4, target.max = 0.7,
                            nattempts = 400)
    acc.mu <- this.update$acc
    att.mu <- this.update$att
    MH.mu  <- this.update$MH
  }

  if (iter %% 500 == 0) {
    par(mfrow = c(3, 3))
    start <- max(1, iter - 20000)
    for (i in 1:3) {
      for (j in 1:3) {
        plot(mu.keep[start:iter, i, j], type = "l",
             main = paste("mu: ", round(mu.t[i, j], 3)),
             ylab = round(acc.mu[i, j] / att.mu[i, j], 3),
             xlab = MH.mu[i, j])
      }
    }
  }
}


#### testing mu and tau ####
rm(list=ls())
source("../../../usefulR/usefulfunctions.R", chdir = TRUE)
source("auxfunctions.R", chdir = TRUE)
source("updatemodel.R", chdir = TRUE)

set.seed(2000)
ns <- 400
nt <- 12
np <- 6
X <- rX(ns, nt, np)
beta.t <- rnorm(np, 0, 10)

phi.t <- 0.2
s <- cbind(runif(ns), runif(ns))
d <- rdist(s)
Sigma.t <- exp(-d / phi.t)
tau.t   <- rgamma(nt, 1, 1)
Qb.t    <- chol2inv(chol(Sigma.t))
Xb.t    <- getXBeta(X = X, beta = beta.t)

mu.t <- matrix(0, ns, nt)
for (t in 1:nt) {
  mu.t[, t] <- Xb.t[, t] + t(chol(Sigma.t)) %*% rnorm(ns) / sqrt(tau.t[t])
}

ls.t <- matrix(0, ns, nt)
xi.t <- 0.01
y.t <- rgev(n = ns * nt, loc = mu.t, exp(ls.t), xi.t)

# initialize values
mu <- matrix(mu.t + rnorm(ns * nt), ns, nt)
SS <- getGPSS(Qb = Qb.t, param = mu, Xb = Xb.t)
curll <- matrix(0, ns, nt)
for (t in 1:nt) {
  curll[, t] <- dgev(x = y.t[, t], loc = mu[, t], exp(ls.t[, t]), xi.t,
                     log = TRUE)
}

niters <- 60000
burn   <- 50000
mu.keep <- array(0, dim = c(niters, ns, nt))
tau <- rep(1, nt)
tau.keep <- matrix(0, niters, nt)
acc.mu <- att.mu <- MH.mu <- matrix(0.1, ns, nt)

set.seed(3366)  # demo
for (iter in 1:niters) {
  this.update <- updateMuTest(mu = mu, Qb = Qb.t, tau = tau, Xb = Xb.t,
                              y = y.t, ls = ls.t, xi = xi.t,
                              SS = SS, curll = curll, acc = acc.mu,
                              att = att.mu, MH = MH.mu)
  mu <- this.update$mu
  SS <- this.update$SS
  curll <- this.update$curll
  acc.mu <- this.update$acc
  att.mu <- this.update$att
  mu.keep[iter, , ] <- mu

  this.update <- updateGPTau(SS = SS, tau.a = 0.1, tau.b = 0.1,
                             ns = ns)
  tau <- this.update$tau
  tau.keep[iter, ] <- tau

  if (iter < burn / 2) {
    this.update <- mhUpdate(acc = acc.mu, att = att.mu, MH = MH.mu,
                            target.min = 0.5, target.max = 0.7,
                            nattempts = 200)
    acc.mu <- this.update$acc
    att.mu <- this.update$att
    MH.mu  <- this.update$MH
  }

  if (iter %% 1000 == 0) {
    par(mfrow = c(4, 3))
    start <- max(1, iter - 20000)
    for (i in 1:3) {
      for (j in 1:3) {
        plot(mu.keep[start:iter, i, j], type = "l",
             main = paste("mu: ", round(mu.t[i, j], 3)),
             ylab = round(acc.mu[i, j] / att.mu[i, j], 3),
             xlab = MH.mu[i, j])
      }
    }
    for (i in 1:3) {
      plot(tau.keep[start:iter, i], type = "l",
           main = paste("tau: ", round(tau.t[i], 3)))
    }
  }
}

#### testing mu, tau, and beta ####
rm(list=ls())
source("../../../usefulR/usefulfunctions.R", chdir = TRUE)
source("auxfunctions.R", chdir = TRUE)
source("updatemodel.R", chdir = TRUE)

set.seed(2000)
ns <- 400
nt <- 12
np <- 6
X <- rX(ns, nt, np)
beta.t <- rnorm(np, 0, 10)

phi.t <- 0.2
s <- cbind(runif(ns), runif(ns))
d <- rdist(s)
Sigma.t <- exp(-d / phi.t)
tau.t   <- rgamma(nt, 1, 1)
Qb.t    <- chol2inv(chol(Sigma.t))
Xb.t    <- getXBeta(X = X, beta = beta.t)

mu.t <- matrix(0, ns, nt)
for (t in 1:nt) {
  mu.t[, t] <- Xb.t[, t] + t(chol(Sigma.t)) %*% rnorm(ns) / sqrt(tau.t[t])
}

ls.t <- matrix(0, ns, nt)
xi.t <- 0.01
y.t <- rgev(n = ns * nt, loc = mu.t, exp(ls.t), xi.t)

# initialize values
mu <- matrix(mu.t + rnorm(ns * nt), ns, nt)
SS <- getGPSS(Qb = Qb.t, param = mu, Xb = Xb.t)
curll <- matrix(0, ns, nt)
for (t in 1:nt) {
  curll[, t] <- dgev(x = y.t[, t], loc = mu[, t], exp(ls.t[, t]), xi.t,
                     log = TRUE)
}

niters <- 20000
burn   <- 15000
beta.sd <- 100
beta <- rep(0, np)
beta.keep <- matrix(0, niters, np)
beta.sd.keep <- rep(0, niters)
mu.keep <- array(0, dim = c(niters, ns, nt))
tau <- rep(1, nt)
tau.keep <- matrix(0, niters, nt)
acc.mu <- att.mu <- MH.mu <- matrix(0.1, ns, nt)

set.seed(3366)  # demo
for (iter in 1:niters) {
  this.update <- updateGPBeta(beta.sd = beta.sd, Qb = Qb.t,
                              param = mu, X = X, SS = SS, tau = tau)
  beta <- this.update$beta
  Xb   <- this.update$Xb
  SS   <- this.update$SS
  beta.keep[iter, ] <- beta

  this.update <- updateGPBetaSD(beta = beta, tau.a = 0.1, tau.b = 1)
  beta.sd <- this.update$beta.sd

  beta.sd.keep[iter] <- beta.sd

  this.update <- updateMuTest(mu = mu, Qb = Qb.t, tau = tau.t, Xb = Xb,
                              y = y.t, ls = ls.t, xi = xi.t,
                              SS = SS, curll = curll, acc = acc.mu,
                              att = att.mu, MH = MH.mu)
  mu <- this.update$mu
  SS <- this.update$SS
  curll <- this.update$curll
  acc.mu <- this.update$acc
  att.mu <- this.update$att
  mu.keep[iter, , ] <- mu

  this.update <- updateGPTau(SS = SS, tau.a = 0.1, tau.b = 0.1,
                             ns = ns)
  tau <- this.update$tau
  tau.keep[iter, ] <- tau

  if (iter < burn / 2) {
    this.update <- mhUpdate(acc = acc.mu, att = att.mu, MH = MH.mu,
                            target.min = 0.5, target.max = 0.7,
                            nattempts = 200)
    acc.mu <- this.update$acc
    att.mu <- this.update$att
    MH.mu  <- this.update$MH
  }

  if (iter %% 1000 == 0) {
    par(mfrow = c(4, 3))
    start <- max(1, iter - 20000)
    for (i in 1:2) {
      for (j in 1:3) {
        plot(mu.keep[start:iter, i, j], type = "l",
             main = paste("mu: ", round(mu.t[i, j], 3)),
             ylab = round(acc.mu[i, j] / att.mu[i, j], 3),
             xlab = MH.mu[i, j])
      }
    }

    for (i in 1:3) {
      plot(beta.keep[start:iter, i], type = "l",
           main = paste("beta: ", round(beta.t[i], 3)))
    }

    for (i in 1:3) {
      plot(tau.keep[start:iter, i], type = "l",
           main = paste("tau: ", round(tau.t[i], 3)))
    }
  }
}

#### Verify gradients - no residual dependence ####
rm(list=ls())
source("../../../usefulR/usefulfunctions.R", chdir = TRUE)
source("auxfunctions.R", chdir = TRUE)
source("updatemodel.R", chdir = TRUE)

set.seed(2000)
ns <- 10
nt <- 3
np <- 6
X1 <- rX(ns, nt, np)
X2 <- rX(ns, nt, np)
beta1.t <- rnorm(np, 0, 1)
beta2.t <- rnorm(np, 0, 0.1)

phi.t <- 0.2
s <- cbind(runif(ns), runif(ns))
d <- rdist(s)
Sigma.t <- exp(-d / phi.t)
tau.t   <- rgamma(nt, 1, 1)
Qb.t    <- chol2inv(chol(Sigma.t))

Xb1.t <- getXBeta(X = X1, beta = beta1.t)
Xb2.t <- getXBeta(X = X2, beta = beta2.t)

mu.t <- ls.t <- matrix(0, ns, nt)
for (t in 1:nt) {
  mu.t[, t] <- Xb1.t[, t] + t(chol(Sigma.t)) %*% rnorm(ns) / sqrt(tau.t[t])
  ls.t[, t] <- Xb2.t[, t] + t(chol(Sigma.t)) %*% rnorm(ns) / sqrt(tau.t[t])
}

xi.t <- 0.1
y.t <- rgev(n = ns * nt, loc = mu.t, scale = exp(ls.t), xi.t)

lp.mu <- logpost.mu.test(mu = mu.t[, t], Xb = Xb1.t[, t], tau = tau.t[t],
                         Qb = Qb.t, y = y.t[, t], ls = ls.t[, t], xi = xi.t)

mean(grad(func = logpost.mu.test, x = mu.t[, t], Xb = Xb1.t[, t],
          tau = tau.t[t], Qb = Qb.t, y = y.t[, t], ls = ls.t[, t], xi = xi.t) /
       logpost.mu.grad.test(mu = mu.t[, t], Xb = Xb1.t[, t], tau = tau.t[t],
                            Qb = Qb.t, y = y.t[, t], ls = ls.t[, t],
                            xi = xi.t))

sd(grad(func = logpost.mu.test, x = mu.t[, t], Xb = Xb1.t[, t], tau = tau.t[t],
        Qb = Qb.t, y = y.t[, t], ls = ls.t[, t], xi = xi.t) /
     logpost.mu.grad.test(mu = mu.t[, t], Xb = Xb1.t[, t], tau = tau.t[t],
                          Qb = Qb.t, y = y.t[, t], ls = ls.t[, t],
                          xi = xi.t))

lp.logsig <- logpost.logsig.test(ls = ls.t[, t], Xb = Xb2.t[, t],
                                 tau = tau.t[t], Qb = Qb.t, y = y.t[, t],
                                 mu = mu.t[, t], xi = xi.t)

mean(grad(func = logpost.logsig.test, x = ls.t[, t], Xb = Xb2.t[, t],
          tau = tau.t[t], Qb = Qb.t, y = y.t[, t], mu = mu.t[, t], xi = xi.t) /
       logpost.logsig.grad.test(ls = ls.t[, t], Xb = Xb2.t[, t],
                                tau = tau.t[t], Qb = Qb.t, y = y.t[, t],
                                mu = mu.t[, t], xi = xi.t))

sd(grad(func = logpost.logsig.test, x = ls.t[, t], Xb = Xb2.t[, t],
        tau = tau.t[t], Qb = Qb.t, y = y.t[, t], mu = mu.t[, t], xi = xi.t) /
     logpost.logsig.grad.test(ls = ls.t[, t], Xb = Xb2.t[, t],
                              tau = tau.t[t], Qb = Qb.t, y = y.t[, t],
                              mu = mu.t[, t], xi = xi.t))

#### testing logsig ####
rm(list=ls())
source("../../../usefulR/usefulfunctions.R", chdir = TRUE)
source("auxfunctions.R", chdir = TRUE)
source("updatemodel.R", chdir = TRUE)

set.seed(2000)
ns <- 400
nt <- 12
np <- 6
X1.t <- rX(ns, nt, np)
X2.t <- rX(ns, nt, np)
beta1.t <- rnorm(np, 0, 10)
beta2.t <- rnorm(np, 0, 5)

phi.t <- 0.2
s <- cbind(runif(ns), runif(ns))
d <- rdist(s)
Sigma.t <- exp(-d / phi.t)
tau.t   <- rgamma(nt, 1, 1)
Qb.t    <- chol2inv(chol(Sigma.t))
Xb1.t   <- getXBeta(X = X1.t, beta = beta1.t)
Xb2.t   <- getXBeta(X = X2.t, beta = beta2.t)

mu.t <- ls.t <- matrix(0, ns, nt)
for (t in 1:nt) {
  mu.t[, t] <- Xb1.t[, t] + t(chol(Sigma.t)) %*% rnorm(ns) / sqrt(tau.t[t])
  ls.t[, t] <- Xb2.t[, t] + t(chol(Sigma.t)) %*% rnorm(ns) / sqrt(tau.t[t])
}

xi.t <- 0.01
y.t <- rgev(n = ns * nt, loc = mu.t, exp(ls.t), xi.t)

# initialize values
ls <- matrix(ls.t + rnorm(ns * nt, 0, 0.1), ns, nt)
SS <- getGPSS(Qb = Qb.t, param = ls, Xb = Xb2.t)
curll <- matrix(0, ns, nt)
for (t in 1:nt) {
  curll[, t] <- dgev(x = y.t[, t], loc = mu.t[, t], exp(ls[, t]), xi.t,
                     log = TRUE)
}

niters <- 10000
burn   <- 8000
ls.keep <- array(0, dim = c(niters, ns, nt))
acc.ls <- att.ls <- MH.ls <- matrix(0.1, ns, nt)

set.seed(3366)  # demo
for (iter in 1:niters) {
  this.update <- updateLSTest(ls = ls, tau = tau.t, Xb = Xb2.t, SS = SS,
                              y = y.t, mu = mu.t, xi = xi.t,
                              Qb = Qb.t, curll = curll, acc = acc.ls,
                              att = att.ls, MH = MH.ls)
  ls <- this.update$ls
  SS <- this.update$SS
  curll <- this.update$curll
  acc.ls <- this.update$acc
  att.ls <- this.update$att
  ls.keep[iter, , ] <- ls

  if (iter < burn / 2) {
    this.update <- mhUpdate(acc = acc.ls, att = att.ls, MH = MH.ls,
                            target.min = 0.4, target.max = 0.7,
                            nattempts = 400)
    acc.ls <- this.update$acc
    att.ls <- this.update$att
    MH.ls  <- this.update$MH
  }

  if (iter %% 500 == 0) {
    par(mfrow = c(3, 3))
    start <- max(1, iter - 20000)
    for (i in 1:3) {
      for (j in 1:3) {
        plot(ls.keep[start:iter, i, j], type = "l",
             main = paste("logsig: ", round(ls.t[i, j], 3)),
             ylab = round(acc.ls[i, j] / att.ls[i, j], 3),
             xlab = MH.ls[i, j])
      }
    }
  }
}

#### testing logsig, tau, and beta ####
rm(list=ls())
source("../../../usefulR/usefulfunctions.R", chdir = TRUE)
source("auxfunctions.R", chdir = TRUE)
source("updatemodel.R", chdir = TRUE)

set.seed(2000)
ns <- 400
nt <- 12
np <- 6
X1.t <- rX(ns, nt, np)
X2.t <- rX(ns, nt, np)
beta1.t <- rnorm(np, 0, 10)
beta2.t <- rnorm(np, 0, 5)

phi.t <- 0.2
s <- cbind(runif(ns), runif(ns))
d <- rdist(s)
Sigma.t <- exp(-d / phi.t)
tau.t   <- rgamma(nt, 1, 1)
Qb.t    <- chol2inv(chol(Sigma.t))
Xb1.t   <- getXBeta(X = X1.t, beta = beta1.t)
Xb2.t   <- getXBeta(X = X2.t, beta = beta2.t)

mu.t <- ls.t <- matrix(0, ns, nt)
for (t in 1:nt) {
  mu.t[, t] <- Xb1.t[, t] + t(chol(Sigma.t)) %*% rnorm(ns) / sqrt(tau.t[t])
  ls.t[, t] <- Xb2.t[, t] + t(chol(Sigma.t)) %*% rnorm(ns) / sqrt(tau.t[t])
}

xi.t <- 0.01
y.t <- rgev(n = ns * nt, loc = mu.t, exp(ls.t), xi.t)

# initialize values
ls <- matrix(ls.t + rnorm(ns * nt, 0, 0.1), ns, nt)
SS <- getGPSS(Qb = Qb.t, param = ls, Xb = Xb2.t)
curll <- matrix(0, ns, nt)
for (t in 1:nt) {
  curll[, t] <- dgev(x = y.t[, t], loc = mu.t[, t], exp(ls[, t]), xi.t,
                     log = TRUE)
}

niters <- 30000
burn   <- 25000
beta.sd <- 100
beta <- rep(0, np)
beta.keep <- matrix(0, niters, np)
beta.sd.keep <- rep(0, niters)
ls.keep <- array(0, dim = c(niters, ns, nt))
tau <- rep(1, nt)
tau.keep <- matrix(0, niters, nt)
acc.ls <- att.ls <- MH.ls <- matrix(0.1, ns, nt)

set.seed(3366)  # demo
for (iter in 1:niters) {
  this.update <- updateGPBeta(beta.sd = beta.sd, Qb = Qb.t,
                              param = ls, X = X2.t, SS = SS, tau = tau)
  beta <- this.update$beta
  Xb2  <- this.update$Xb
  SS   <- this.update$SS
  beta.keep[iter, ] <- beta

  this.update <- updateGPBetaSD(beta = beta, tau.a = 0.1, tau.b = 1)
  beta.sd <- this.update$beta.sd

  beta.sd.keep[iter] <- beta.sd

  this.update <- updateLSTest(ls = ls, tau = tau.t, Xb = Xb2, SS = SS,
                              y = y.t, mu = mu.t, xi = xi.t,
                              Qb = Qb.t, curll = curll,
                              acc = acc.ls, att = att.ls, MH = MH.ls)
  ls    <- this.update$ls
  SS    <- this.update$SS
  curll <- this.update$curll
  acc.ls <- this.update$acc
  att.ls <- this.update$att
  ls.keep[iter, , ] <- ls

  this.update <- updateGPTau(SS = SS, tau.a = 0.1, tau.b = 0.1,
                             ns = ns)
  tau <- this.update$tau
  tau.keep[iter, ] <- tau

  if (iter < burn / 2) {
    this.update <- mhUpdate(acc = acc.ls, att = att.ls, MH = MH.ls,
                            target.min = 0.4, target.max = 0.7,
                            nattempts = 400)
    acc.ls <- this.update$acc
    att.ls <- this.update$att
    MH.ls  <- this.update$MH
  }

  if (iter %% 1000 == 0) {
    par(mfrow = c(4, 3))
    start <- max(1, iter - 20000)
    for (i in 1:2) {
      for (j in 1:3) {
        plot(ls.keep[start:iter, i, j], type = "l",
             main = paste("logsig: ", round(ls.t[i, j], 3)),
             ylab = round(acc.ls[i, j] / att.ls[i, j], 3),
             xlab = MH.ls[i, j])
      }
    }

    for (i in 1:3) {
      plot(beta.keep[start:iter, i], type = "l",
           main = paste("beta: ", round(beta2.t[i], 3)))
    }

    for (i in 1:3) {
      plot(tau.keep[start:iter, i], type = "l",
           main = paste("tau: ", round(tau.t[i], 3)))
    }
  }
}

#### testing basis bandwidth update ####
rm(list=ls())
source("../../../usefulR/usefulfunctions.R", chdir = TRUE)
source("auxfunctions.R", chdir = TRUE)
source("updatemodel.R", chdir = TRUE)

set.seed(2000)
# setting np later after X is created
ns <- 400
nt <- 12
nknots <- 5
time.interact <- TRUE

s <- cbind(runif(ns), runif(ns))
knots <- as.matrix(cover.design(R = s, nd = nknots)$design)
d <- rdist(s)
dw2 <- rdist(s, knots)^2
dw2[dw2 < 1e-4] <- 0

# create the matrix of covariates
X1.t <- X2.t <- array(1, dim = c(ns, nt, 2))
for (t in 1:nt) {
  time <- (t - nt / 2) / nt
  X1.t[, t, 2] <- X2.t[, t, 2] <- time
}

bw.basis.t <- 0.2
B.t <- makeW(dw2 = dw2, rho = bw.basis.t)
X1.t <- add.basis.X(X1.t, B.t, time.interact = time.interact)
X2.t <- add.basis.X(X2.t, B.t, time.interact = time.interact)

np <- dim(X1.t)[3]

beta1.t <- rnorm(np, 0, 10)
beta2.t <- rnorm(np, 0, 5)

bw.gp.t <- 0.2
Sigma.t <- exp(-d / bw.gp.t)
tau1.t <- tau2.t <- rgamma(nt, 1, 1)
Qb.t    <- chol2inv(chol(Sigma.t))

Xb1.t <- getXBeta(X = X1.t, beta = beta1.t)
Xb2.t <- getXBeta(X = X2.t, beta = beta2.t)

mu.t <- ls.t <- matrix(0, ns, nt)
for (t in 1:nt) {
  mu.t[, t] <- Xb1.t[, t] + t(chol(Sigma.t)) %*% rnorm(ns) / sqrt(tau1.t[t])
  ls.t[, t] <- Xb2.t[, t] + t(chol(Sigma.t)) %*% rnorm(ns) / sqrt(tau2.t[t])
}

xi.t <- 0.01
y.t <- rgev(n = ns * nt, loc = mu.t, scale = exp(ls.t), shape = xi.t)

# initialize values
bw.basis <- 0.4
bw.basis.min <- quantile(dw2, 0.01)
bw.basis.max <- quantile(dw2, 0.99)
B <- makeW(dw2 = dw2, rho = bw.basis)
X1 <- rep.basis.X(X = X1.t, newB = B, time.interact = time.interact)
X2 <- rep.basis.X(X = X2.t, newB = B, time.interact = time.interact)
Xb1 <- getXBeta(X = X1, beta = beta1.t)
Xb2 <- getXBeta(X = X2, beta = beta1.t)
SS1 <- getGPSS(Qb = Qb.t, param = mu.t, Xb = Xb1)
SS2 <- getGPSS(Qb = Qb.t, param = ls.t, Xb = Xb2)

niters <- 30000
burn   <- 25000

# storage
bw.basis.keep <- rep(0, niters)
Xb1.keep     <- array(0, dim = c(niters, ns, nt))
acc.bw.basis <- att.bw.basis <- MH.bw.basis <- 0.1

set.seed(3366)  # demo
for (iter in 1:niters) {
  this.update <- updateXBasisBW(bw = bw.basis, bw.min = bw.basis.min,
                                bw.max = bw.basis.max,
                                X1 = X1, beta1 = beta1.t, Xb1 = Xb1,
                                mu = mu.t, tau1 = tau1.t, SS1 = SS1,
                                X2 = X2, beta2 = beta2.t, Xb2 = Xb2,
                                ls = ls.t, tau2 = tau2.t, SS2 = SS2,
                                Qb = Qb.t, dw2 = dw2,
                                time.interact = time.interact,
                                acc = acc.bw.basis, att = att.bw.basis,
                                MH = MH.bw.basis)
  bw.basis <- this.update$bw
  X1  <- this.update$X1
  Xb1 <- this.update$Xb1
  SS1 <- this.update$SS1
  X2  <- this.update$X2
  Xb2 <- this.update$Xb2
  SS2 <- this.update$SS2
  acc.bw.basis <- this.update$acc
  att.bw.basis <- this.update$att

  bw.basis.keep[iter] <- bw.basis
  Xb1.keep[iter, , ] <- Xb1

  if (iter < burn / 2) {
    this.update <- mhUpdate(acc = acc.bw.basis, att = att.bw.basis,
                            MH = MH.bw.basis,
                            target.min = 0.3, target.max = 0.6,
                            lower = 0.8, higher = 1.2)
    acc.bw.basis <- this.update$acc
    att.bw.basis <- this.update$att
    MH.bw.basis  <- this.update$MH
  }

  if (iter %% 1000 == 0) {
    par(mfrow = c(2, 5))
    plot(bw.basis.keep[100:iter], type = "l",
         main = paste("BW basis: ", bw.basis.t, sep = ""),
         ylab = round(acc.bw.basis / att.bw.basis, 3),
         xlab = MH.bw.basis)
    for(i in 1:3) { for (j in 1:3) {
      plot(Xb1.keep[100:iter, i, j], type = "l",
           main = paste("Xb1: ", Xb1.t[i, j], sep = ""))
    }}
  }
}

#### testing basis bandwidth, beta1, beta2 ####
rm(list=ls())
source("../../../usefulR/usefulfunctions.R", chdir = TRUE)
source("auxfunctions.R", chdir = TRUE)
source("updatemodel.R", chdir = TRUE)

set.seed(2000)
# setting np later after X is created
ns <- 400
nt <- 12
nknots <- 5
time.interact <- TRUE

s <- cbind(runif(ns), runif(ns))
knots <- as.matrix(cover.design(R = s, nd = nknots)$design)
d <- rdist(s)
dw2 <- rdist(s, knots)^2
dw2[dw2 < 1e-4] <- 0

# create the matrix of covariates
X1.t <- X2.t <- array(1, dim = c(ns, nt, 2))
for (t in 1:nt) {
  time <- (t - nt / 2) / nt
  X1.t[, t, 2] <- X2.t[, t, 2] <- time
}

bw.basis.t <- 0.2
B.t <- makeW(dw2 = dw2, rho = bw.basis.t)
X1.t <- add.basis.X(X1.t, B.t, time.interact = time.interact)
X2.t <- add.basis.X(X2.t, B.t, time.interact = time.interact)

np <- dim(X1.t)[3]

beta1.t <- rnorm(np, 0, 10)
beta2.t <- rnorm(np, 0, 5)

bw.gp.t <- 0.2
Sigma.t <- exp(-d / bw.gp.t)
tau1.t <- tau2.t <- rgamma(nt, 1, 1)
Qb.t    <- chol2inv(chol(Sigma.t))

Xb1.t <- getXBeta(X = X1.t, beta = beta1.t)
Xb2.t <- getXBeta(X = X2.t, beta = beta2.t)

mu.t <- ls.t <- matrix(0, ns, nt)
for (t in 1:nt) {
  mu.t[, t] <- Xb1.t[, t] + t(chol(Sigma.t)) %*% rnorm(ns) / sqrt(tau1.t[t])
  ls.t[, t] <- Xb2.t[, t] + t(chol(Sigma.t)) %*% rnorm(ns) / sqrt(tau2.t[t])
}

xi.t <- 0.01
y.t <- rgev(n = ns * nt, loc = mu.t, scale = exp(ls.t), shape = xi.t)

# initialize values
beta1 <- beta2 <- rep(0, np)
beta1.sd <- beta2.sd <- 100
bw.basis <- 0.4
bw.basis.min <- quantile(dw2, 0.01)
bw.basis.max <- quantile(dw2, 0.99)
B <- makeW(dw2 = dw2, rho = bw.basis)
X1 <- rep.basis.X(X = X1.t, newB = B, time.interact = time.interact)
X2 <- rep.basis.X(X = X2.t, newB = B, time.interact = time.interact)
Xb1 <- getXBeta(X = X1, beta = beta1)
Xb2 <- getXBeta(X = X2, beta = beta1)
SS1 <- getGPSS(Qb = Qb.t, param = mu.t, Xb = Xb1)
SS2 <- getGPSS(Qb = Qb.t, param = ls.t, Xb = Xb2)

niters <- 30000
burn   <- 25000

# storage
bw.basis.keep <- rep(0, niters)
beta1.keep <- beta2.keep <- matrix(0, niters, np)
beta1.sd.keep <- beta2.sd.keep <- rep(0, niters)
acc.bw.basis <- att.bw.basis <- MH.bw.basis <- 0.1

set.seed(3366)  # demo
for (iter in 1:niters) {
  this.update <- updateGPBeta(beta.sd = beta1.sd, Qb = Qb.t,
                              param = mu.t, X = X1, SS = SS1, tau = tau1.t)
  beta1 <- this.update$beta
  Xb1   <- this.update$Xb
  SS1   <- this.update$SS
  beta1.keep[iter, ] <- beta1

  this.update <- updateGPBetaSD(beta = beta1, tau.a = 0.5, tau.b = 0.5)
  beta1.sd <- this.update$beta.sd
  beta1.sd.keep[iter] <- beta1.sd

  this.update <- updateGPBeta(beta.sd = beta2.sd, Qb = Qb.t,
                              param = ls.t, X = X2, SS = SS2, tau = tau2.t)
  beta2 <- this.update$beta
  Xb2   <- this.update$Xb
  SS2   <- this.update$SS
  beta2.keep[iter, ] <- beta2

  this.update <- updateGPBetaSD(beta = beta2, tau.a = 0.5, tau.b = 0.5)
  beta2.sd <- this.update$beta.sd
  beta2.sd.keep[iter] <- beta2.sd

  this.update <- updateXBasisBW(bw = bw.basis, bw.min = bw.basis.min,
                                bw.max = bw.basis.max,
                                X1 = X1, beta1 = beta1, Xb1 = Xb1,
                                mu = mu.t, tau1 = tau1.t, SS1 = SS1,
                                X2 = X2, beta2 = beta2, Xb2 = Xb2,
                                ls = ls.t, tau2 = tau2.t, SS2 = SS2,
                                Qb = Qb.t, dw2 = dw2,
                                time.interact = time.interact,
                                acc = acc.bw.basis, att = att.bw.basis,
                                MH = MH.bw.basis)
  bw.basis <- this.update$bw
  X1  <- this.update$X1
  Xb1 <- this.update$Xb1
  SS1 <- this.update$SS1
  X2  <- this.update$X2
  Xb2 <- this.update$Xb2
  SS2 <- this.update$SS2
  acc.bw.basis <- this.update$acc
  att.bw.basis <- this.update$att

  bw.basis.keep[iter] <- bw.basis

  if (iter < burn / 2) {
    this.update <- mhUpdate(acc = acc.bw.basis, att = att.bw.basis,
                            MH = MH.bw.basis,
                            target.min = 0.3, target.max = 0.6,
                            lower = 0.8, higher = 1.2)
    acc.bw.basis <- this.update$acc
    att.bw.basis <- this.update$att
    MH.bw.basis  <- this.update$MH
  }

  if (iter %% 500 == 0) {
    par(mfrow = c(3, 3))
    plot(bw.basis.keep[100:iter], type = "l",
         main = paste("BW basis: ", bw.basis.t, sep = ""),
         ylab = round(acc.bw.basis / att.bw.basis, 3),
         xlab = MH.bw.basis)
    for(i in 1:4) {
      plot(beta1.keep[100:iter, i], type = "l",
           main = paste("beta 1: ", beta1.t[i], sep = ""))
    }
    for (i in 1:4) {
      plot(beta2.keep[100:iter, i], type = "l",
           main = paste("beta 2: ", beta2.t[i], sep = ""))
    }

  }
}

#### testing basis bandwidth, beta1, beta2, mu, and logsig ####
rm(list=ls())
source("../../../usefulR/usefulfunctions.R", chdir = TRUE)
source("auxfunctions.R", chdir = TRUE)
source("updatemodel.R", chdir = TRUE)

set.seed(2000)
# setting np later after X is created
ns <- 400
nt <- 12
nknots <- 5
time.interact <- TRUE

s <- cbind(runif(ns), runif(ns))
knots <- as.matrix(cover.design(R = s, nd = nknots)$design)
d <- rdist(s)
dw2 <- rdist(s, knots)^2
dw2[dw2 < 1e-4] <- 0

# create the matrix of covariates
X1.t <- X2.t <- array(1, dim = c(ns, nt, 2))
for (t in 1:nt) {
  time <- (t - nt / 2) / nt
  X1.t[, t, 2] <- X2.t[, t, 2] <- time
}

bw.basis.t <- 0.2
B.t <- makeW(dw2 = dw2, rho = bw.basis.t)
X1.t <- add.basis.X(X1.t, B.t, time.interact = time.interact)
X2.t <- add.basis.X(X2.t, B.t, time.interact = time.interact)

np <- dim(X1.t)[3]

beta1.t <- rnorm(np, 0, 10)
beta2.t <- rnorm(np, 0, 5)

bw.gp.t <- 0.2
Sigma.t <- exp(-d / bw.gp.t)
tau1.t <- tau2.t <- rgamma(nt, 1, 1)
Qb.t    <- chol2inv(chol(Sigma.t))

Xb1.t <- getXBeta(X = X1.t, beta = beta1.t)
Xb2.t <- getXBeta(X = X2.t, beta = beta2.t)

mu.t <- ls.t <- matrix(0, ns, nt)
for (t in 1:nt) {
  mu.t[, t] <- Xb1.t[, t] + t(chol(Sigma.t)) %*% rnorm(ns) / sqrt(tau1.t[t])
  ls.t[, t] <- Xb2.t[, t] + t(chol(Sigma.t)) %*% rnorm(ns) / sqrt(tau2.t[t])
}

xi.t <- 0.01
y.t <- rgev(n = ns * nt, loc = mu.t, scale = exp(ls.t), shape = xi.t)

# initialize values
mu <- matrix(mu.t + rnorm(ns * nt), ns, nt)
ls <- matrix(ls.t + rnorm(ns * nt), ns, nt)
beta1 <- beta2 <- rep(0, np)
beta1.sd <- beta2.sd <- 100
bw.basis <- 0.4
bw.basis.min <- quantile(dw2, 0.01)
bw.basis.max <- quantile(dw2, 0.99)
B <- makeW(dw2 = dw2, rho = bw.basis)
X1 <- rep.basis.X(X = X1.t, newB = B, time.interact = time.interact)
X2 <- rep.basis.X(X = X2.t, newB = B, time.interact = time.interact)
Xb1 <- getXBeta(X = X1, beta = beta1)
Xb2 <- getXBeta(X = X2, beta = beta1)
SS1 <- getGPSS(Qb = Qb.t, param = mu.t, Xb = Xb1)
SS2 <- getGPSS(Qb = Qb.t, param = ls.t, Xb = Xb2)

curll <- matrix(0, ns, nt)
for (t in 1:nt) {
  curll[, t] <- dgev(x = y.t[, t], loc = mu[, t], exp(ls[, t]), xi.t,
                     log = TRUE)
}

niters <- 30000
burn   <- 25000

# storage
bw.basis.keep <- rep(0, niters)
beta1.keep <- beta2.keep <- matrix(0, niters, np)
beta1.sd.keep <- beta2.sd.keep <- rep(0, niters)
mu.keep <- array(0, dim = c(niters, ns, nt))
ls.keep <- array(0, dim = c(niters, ns, nt))
acc.bw.basis <- att.bw.basis <- MH.bw.basis <- 0.1
acc.mu <- att.mu <- MH.mu <- matrix(0.2, ns, nt)
acc.ls <- att.ls <- MH.ls <- matrix(0.1, ns, nt)

set.seed(3366)  # demo
for (iter in 1:niters) {
  this.update <- updateGPBeta(beta.sd = beta1.sd, Qb = Qb.t,
                              param = mu, X = X1, SS = SS1, tau = tau1.t)
  beta1 <- this.update$beta
  Xb1   <- this.update$Xb
  SS1   <- this.update$SS
  beta1.keep[iter, ] <- beta1

  this.update <- updateGPBetaSD(beta = beta1, tau.a = 0.5, tau.b = 0.5)
  beta1.sd <- this.update$beta.sd
  beta1.sd.keep[iter] <- beta1.sd

  this.update <- updateGPBeta(beta.sd = beta2.sd, Qb = Qb.t,
                              param = ls, X = X2, SS = SS2, tau = tau2.t)
  beta2 <- this.update$beta
  Xb2   <- this.update$Xb
  SS2   <- this.update$SS
  beta2.keep[iter, ] <- beta2

  this.update <- updateGPBetaSD(beta = beta2, tau.a = 0.5, tau.b = 0.5)
  beta2.sd <- this.update$beta.sd
  beta2.sd.keep[iter] <- beta2.sd

  this.update <- updateXBasisBW(bw = bw.basis, bw.min = bw.basis.min,
                                bw.max = bw.basis.max,
                                X1 = X1, beta1 = beta1, Xb1 = Xb1,
                                mu = mu, tau1 = tau1.t, SS1 = SS1,
                                X2 = X2, beta2 = beta2, Xb2 = Xb2,
                                ls = ls, tau2 = tau2.t, SS2 = SS2,
                                Qb = Qb.t, dw2 = dw2,
                                time.interact = time.interact,
                                acc = acc.bw.basis, att = att.bw.basis,
                                MH = MH.bw.basis)
  bw.basis <- this.update$bw
  X1  <- this.update$X1
  Xb1 <- this.update$Xb1
  SS1 <- this.update$SS1
  X2  <- this.update$X2
  Xb2 <- this.update$Xb2
  SS2 <- this.update$SS2
  acc.bw.basis <- this.update$acc
  att.bw.basis <- this.update$att
  bw.basis.keep[iter] <- bw.basis

  this.update <- updateMuTest(mu = mu, tau = tau1.t, Xb = Xb1, SS = SS1,
                              y = y.t, ls = ls, xi = xi.t,
                              Qb = Qb.t, curll = curll, acc = acc.mu,
                              att = att.mu, MH = MH.mu)
  mu    <- this.update$mu
  SS1   <- this.update$SS
  curll <- this.update$curll
  acc.mu <- this.update$acc
  att.mu <- this.update$att
  mu.keep[iter, , ] <- mu

  this.update <- updateLSTest(ls = ls, tau = tau2.t, Xb = Xb2, SS = SS2,
                              y = y.t, mu = mu, xi = xi.t,
                              Qb = Qb.t, curll = curll, acc = acc.ls,
                              att = att.ls, MH = MH.ls)
  ls     <- this.update$ls
  SS2    <- this.update$SS
  curll  <- this.update$curll
  acc.ls <- this.update$acc
  att.ls <- this.update$att
  ls.keep[iter, , ] <- ls

  if (iter < burn / 2) {
    this.update <- mhUpdate(acc = acc.bw.basis, att = att.bw.basis,
                            MH = MH.bw.basis,
                            target.min = 0.3, target.max = 0.6,
                            lower = 0.8, higher = 1.2)
    acc.bw.basis <- this.update$acc
    att.bw.basis <- this.update$att
    MH.bw.basis  <- this.update$MH

    this.update <- mhUpdate(acc = acc.mu, att = att.mu, MH = MH.mu,
                            target.min = 0.4, target.max = 0.7,
                            nattempts = 400)
    acc.mu <- this.update$acc
    att.mu <- this.update$att
    MH.mu  <- this.update$MH

    this.update <- mhUpdate(acc = acc.ls, att = att.ls, MH = MH.ls,
                            target.min = 0.4, target.max = 0.7,
                            nattempts = 400)
    acc.ls <- this.update$acc
    att.ls <- this.update$att
    MH.ls  <- this.update$MH
  }

  if (iter %% 500 == 0) {
    par(mfrow = c(3, 3))
    plot(bw.basis.keep[100:iter], type = "l",
         main = paste("BW basis: ", bw.basis.t, sep = ""),
         ylab = round(acc.bw.basis / att.bw.basis, 3),
         xlab = MH.bw.basis)
    for(i in 1:4) {
      plot(beta1.keep[100:iter, i], type = "l",
           main = paste("beta 1: ", beta1.t[i], sep = ""))
    }
    for (i in 1:4) {
      plot(beta2.keep[100:iter, i], type = "l",
           main = paste("beta 2: ", beta2.t[i], sep = ""))
    }

  }
}


#### Verify gradients - with residual dependence ####
rm(list=ls())
source("../../../usefulR/usefulfunctions.R", chdir = TRUE)
source("auxfunctions.R", chdir = TRUE)
source("updatemodel.R", chdir = TRUE)

set.seed(2000)
ns <- 10
nt <- 3
np <- 6
X1 <- rX(ns, nt, np)
X2 <- rX(ns, nt, np)
beta1.t <- rnorm(np, 0, 1)
beta2.t <- rnorm(np, 0, 0.1)

phi.t <- 0.2
s <- cbind(runif(ns), runif(ns))
d <- rdist(s)
Sigma.t <- exp(-d / phi.t)
tau.t   <- rgamma(nt, 1, 1)
Qb.t    <- chol2inv(chol(Sigma.t))

Xb1.t <- getXBeta(X = X1, beta = beta1.t)
Xb2.t <- getXBeta(X = X2, beta = beta2.t)

mu.t <- ls.t <- matrix(0, ns, nt)
for (t in 1:nt) {
  mu.t[, t] <- Xb1.t[, t] + t(chol(Sigma.t)) %*% rnorm(ns) / sqrt(tau.t[t])
  ls.t[, t] <- Xb2.t[, t] + t(chol(Sigma.t)) %*% rnorm(ns) / sqrt(tau.t[t])
}

xi.t <- 0.1
y.t <- rgev(n = ns * nt, loc = mu.t, scale = exp(ls.t), xi.t)

nknots <- 4
theta.t <- matrix(abs(rnorm(ns * nt)), ns, nt)
alpha.t <- 0.4
thresh.t <- matrix(median(y.t), ns, nt)

xi.t <- 0
lp.mu <- logpost.mu(mu = mu.t[, t], Xb = Xb1.t[, t], tau = tau.t[t],
                    Qb = Qb.t, y = y.t[, t], ls = ls.t[, t], xi = xi.t,
                    theta = theta.t[, t], thresh = thresh.t[, t],
                    alpha = alpha.t)

mean(grad(func = logpost.mu, x = mu.t[, t], Xb = Xb1.t[, t], tau = tau.t[t],
          Qb = Qb.t, y = y.t[, t], ls = ls.t[, t], xi = xi.t,
          theta = theta.t[, t], thresh = thresh.t[, t], alpha = alpha.t) /
       logpost.mu.grad(mu = mu.t[, t], Xb = Xb1.t[, t], tau = tau.t[t],
                       Qb = Qb.t, y = y.t[, t], ls = ls.t[, t],
                       xi = xi.t, theta = theta.t[, t], thresh = thresh.t[, t],
                       alpha = alpha.t))

sd(grad(func = logpost.mu, x = mu.t[, t], Xb = Xb1.t[, t], tau = tau.t[t],
        Qb = Qb.t, y = y.t[, t], ls = ls.t[, t], xi = xi.t,
        theta = theta.t[, t], thresh = thresh.t[, t], alpha = alpha.t) /
     logpost.mu.grad(mu = mu.t[, t], Xb = Xb1.t[, t], tau = tau.t[t],
                     Qb = Qb.t, y = y.t[, t], ls = ls.t[, t],
                     xi = xi.t, theta = theta.t[, t], thresh = thresh.t[, t],
                     alpha = alpha.t))

lp.logsig <- logpost.logsig(ls = ls.t[, t], Xb = Xb2.t[, t], tau = tau.t[t],
                            Qb = Qb.t, y = y.t[, t], mu = mu.t[, t],
                            xi = xi.t, theta = theta.t[, t],
                            thresh = thresh.t[, t], alpha = alpha.t)

mean(grad(func = logpost.logsig, x = ls.t[, t], Xb = Xb2.t[, t],
          tau = tau.t[t], Qb = Qb.t, y = y.t[, t], mu = mu.t[, t], xi = xi.t,
          theta = theta.t[, t], thresh = thresh.t[, t], alpha = alpha.t) /
       logpost.logsig.grad(ls = ls.t[, t], Xb = Xb2.t[, t],
                           tau = tau.t[t], Qb = Qb.t, y = y.t[, t],
                           mu = mu.t[, t], xi = xi.t,
                           theta = theta.t[, t], thresh = thresh.t[, t],
                           alpha = alpha.t))

sd(grad(func = logpost.logsig, x = ls.t[, t], Xb = Xb2.t[, t],
        tau = tau.t[t], Qb = Qb.t, y = y.t[, t], mu = mu.t[, t], xi = xi.t,
        theta = theta.t[, t], thresh = thresh.t[, t], alpha = alpha.t) /
     logpost.logsig.grad(ls = ls.t[, t], Xb = Xb2.t[, t],
                         tau = tau.t[t], Qb = Qb.t, y = y.t[, t],
                         mu = mu.t[, t], xi = xi.t,
                         theta = theta.t[, t], thresh = thresh.t[, t],
                         alpha = alpha.t))


#### testing xi ####
rm(list=ls())
source("../../../usefulR/usefulfunctions.R", chdir = TRUE)
source("auxfunctions.R", chdir = TRUE)
source("updatemodel.R", chdir = TRUE)

set.seed(2000)
ns <- 400
nt <- 12
np <- 6
X1.t <- rX(ns, nt, np)
X2.t <- rX(ns, nt, np)
beta1.t <- rnorm(np, 0, 10)
beta2.t <- rnorm(np, 0, 5)

phi.t <- 0.2
s <- cbind(runif(ns), runif(ns))
d <- rdist(s)
Sigma.t <- exp(-d / phi.t)
tau.t   <- rgamma(nt, 1, 1)
Qb.t    <- chol2inv(chol(Sigma.t))
Xb1.t   <- getXBeta(X = X1.t, beta = beta1.t)
Xb2.t   <- getXBeta(X = X2.t, beta = beta2.t)

mu.t <- ls.t <- matrix(0, ns, nt)
for (t in 1:nt) {
  mu.t[, t] <- Xb1.t[, t] + t(chol(Sigma.t)) %*% rnorm(ns) / sqrt(tau.t[t])
  ls.t[, t] <- Xb2.t[, t] + t(chol(Sigma.t)) %*% rnorm(ns) / sqrt(tau.t[t])
}

xi.t <- -0.7
y.t <- rgev(n = ns * nt, loc = mu.t, exp(ls.t), xi.t)

# initialize values
xi <- 0.1

thresh.t <- matrix(-Inf, ns, nt)
theta.t <- matrix(1, ns, nt)
alpha.t <- 1
curll <- loglike(y = y.t, mu = mu.t, ls = ls.t, xi = xi,
                 theta = theta.t, thresh = thresh.t, alpha = alpha.t)


niters <- 10000
burn   <- 8000
xi.keep <- rep(0, niters)
acc.xi <- att.xi <- MH.xi <- 0.01


set.seed(3366)  # demo
for (iter in 1:niters) {
  this.update <- updateXi(xi = xi, xi.min = -2, xi.max = 2,
                          xi.mn = 0, xi.sd = 0.5, y = y.t, mu = mu.t,
                          ls = ls.t, curll = curll, theta = theta.t,
                          thresh = thresh.t, alpha = alpha.t, acc = acc.xi,
                          att = att.xi, MH = MH.xi)
  xi <- this.update$xi
  curll <- this.update$curll
  acc.xi <- this.update$acc
  att.xi <- this.update$att
  xi.keep[iter] <- xi

  if (iter < burn / 2) {
    this.update <- mhUpdate(acc = acc.xi, att = att.xi, MH = MH.xi,
                            target.min = 0.3, target.max = 0.6,
                            nattempts = 400)
    acc.xi <- this.update$acc
    att.xi <- this.update$att
    MH.xi  <- this.update$MH
  }

  if (iter %% 500 == 0) {
    start <- max(1, iter - 20000)
    plot(xi.keep[start:iter], type = "l", main = paste("xi: ", xi.t),
         ylab = round(acc.xi / att.xi, 3),
         xlab = MH.xi)
  }
}

#### testing beta1 ####
rm(list=ls())
source("../../../usefulR/usefulfunctions.R", chdir = TRUE)
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
acc.beta1  <- att.beta1 <- MH.beta1 <- matrix(0.1, ns, nt)

set.seed(3366)  # demo
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
                            target.min = 0.3, target.max = 0.6,
                            nattempts = 400)
    acc.beta1 <- this.update$acc
    att.beta1 <- this.update$att
    MH.beta1  <- this.update$MH
  }

  if (iter %% 500 == 0) {
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

#### Verify gradients - on beta ####
rm(list=ls())
source("../../../usefulR/usefulfunctions.R", chdir = TRUE)
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

#### testing beta1 ####
rm(list=ls())
library(compiler)
enableJIT(3)
source("../../../usefulR/usefulfunctions.R", chdir = TRUE)
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
beta1.int.mn <- 0
beta1.time <- beta1.time.t + rnorm(ns * nt)
beta1.time.mn <- 0
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
keep.beta1.int.mn <- keep.beta1.time.mn <- rep(0)
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

rm(list=ls())
library(fields)
library(Rcpp)
library(emulator)
library(microbenchmark)
library(SpatialExtremes)
library(numDeriv)
library(fields)

#### testing beta ####
source("../../../usefulR/usefulfunctions.R", chdir = TRUE)
source("auxfunctions.R", chdir = TRUE)
source("updatemodel.R", chdir = TRUE)

set.seed(2000)
ns <- 400
nt <- 2
np <- 6
X <- rX(ns, nt, np)
beta.t <- rnorm(np, 0, 10)

phi <- 0.2
s <- cbind(runif(ns), runif(ns))
d <- rdist(s)
Sigma <- exp(-d / phi)
tau   <- rgamma(nt, 0.5, 0.5)
Qb <- chol2inv(chol(Sigma))
Xb <- getXBeta(X = X, beta = beta.t)

if (nt == 1) {
  Xb <- matrix(Xb, ns, nt)
}

mu <- matrix(0, ns, nt)
for (t in 1:nt) {
  mu[, t] <- Xb[, t] + t(chol(Sigma)) %*% rnorm(ns) / sqrt(tau[t])
}

# initialize values
SS <- diag(quad.form(Qb, mu - Xb))

niters <- 10000
beta.keep <- matrix(0, niters, np)
beta <- rep(0, np)
set.seed(3366)  # demo
for (iter in 1:niters) {
  this.update <- updateGPBeta(beta.sd = 100, Qb = Qb,
                              param = mu, X = X, SS = SS, tau = tau)
  beta <- this.update$beta
  Xb   <- this.update$Xb
  SS   <- this.update$SS
  beta.keep[iter, ] <- beta

  if (iter %% 500 == 0) {
    start <- max(1, iter - 2000)
    par(mfrow = c(2, np / 2))
    for (i in 1:6) {
      plot(beta.keep[start:iter, i], type = "l",
           main = paste("Beta = ", round(beta.t[i], 3)))
    }
  }
}

#### testing phi ####
rm(list=ls())
source("../../../usefulR/usefulfunctions.R", chdir = TRUE)
source("auxfunctions.R", chdir = TRUE)
source("updatemodel.R", chdir = TRUE)

set.seed(2000)
ns <- 400
nt <- 10
np <- 6
X <- rX(ns, nt, np)
beta.t <- rnorm(np, 0, 10)

phi.t <- 0.2
s <- cbind(runif(ns), runif(ns))
d <- rdist(s)
Sigma.t <- exp(-d / phi.t)
tau   <- rgamma(nt, 0.5, 0.5)
Qb <- chol2inv(chol(Sigma.t))
Xb <- getXBeta(X = X, beta = beta.t)

if (nt == 1) {
  Xb <- matrix(Xb, ns, nt)
}

mu <- matrix(0, ns, nt)
for (t in 1:nt) {
  mu[, t] <- Xb[, t] + t(chol(Sigma.t)) %*% rnorm(ns) / sqrt(tau[t])
}

# initialize values
phi <- 0.05
Qb <- chol2inv(chol(exp(-d / phi)))
SS <- diag(quad.form(Qb, mu - Xb))
phi <- 0.05

niters <- 10000
phi.keep <- rep(0, niters)
acc.phi <- att.phi <- MH.phi <- 0.1
set.seed(3366)  # demo
for (iter in 1:niters) {
  this.update <- updateGPBW(bw = phi, bw.min = 0.01, bw.max = 1.2,
                            bw.mn = 0, bw.sd = 1, Qb = Qb, d = d,
                            mu = mu, Xb1 = Xb, tau1 = tau, SS1 = SS,
                            ls = mu, Xb2 = Xb, tau2 = tau, SS2 = SS,
                            acc = acc.phi, att = att.phi, MH = MH.phi)
  phi <- this.update$bw
  Qb <- this.update$Qb
  SS <- this.update$SS1
  acc.phi <- this.update$acc
  att.phi <- this.update$att

  this.update <- mhUpdate(acc = acc.phi, att = att.phi, MH = MH.phi)
  acc.phi <- this.update$acc
  att.phi <- this.update$att
  MH.phi  <- this.update$MH

  phi.keep[iter] <- phi
  if (iter %% 500 == 0) {
    start <- max(1, iter - 2000)
    plot(phi.keep[start:iter], type = "l")
  }
}

#### testing tau ####
rm(list=ls())
source("../../../usefulR/usefulfunctions.R", chdir = TRUE)
source("auxfunctions.R", chdir = TRUE)
source("updatemodel.R", chdir = TRUE)

set.seed(2000)
ns <- 400
nt <- 12
np <- 6
X <- rX(ns, nt, np)
beta.t <- rnorm(np, 0, 10)

phi.t <- 0.2
s <- cbind(runif(ns), runif(ns))
d <- rdist(s)
Sigma.t <- exp(-d / phi.t)
tau.t   <- rgamma(nt, 0.5, 0.5)
Qb <- chol2inv(chol(Sigma.t))
Xb <- getXBeta(X = X, beta = beta.t)

if (nt == 1) {
  Xb <- matrix(Xb, ns, nt)
}

mu <- matrix(0, ns, nt)
for (t in 1:nt) {
  mu[, t] <- Xb[, t] + t(chol(Sigma.t)) %*% rnorm(ns) / sqrt(tau.t[t])
}

SS <- getGPSS(Qb = Qb, param = mu, Xb = Xb)

# initialize values
niters <- 10000
tau <- rep(1, nt)
tau.keep <- matrix(0, niters, nt)
set.seed(3366)  # demo
for (iter in 1:niters) {
  this.update <- updateGPTau(SS = SS, tau.a = 0.1, tau.b = 0.1,
                             ns = ns)
  tau <- this.update$tau

  tau.keep[iter, ] <- tau
  if (iter %% 500 == 0) {
    start <- max(1, iter - 2000)
    par(mfrow = c(4, 3))
    for (t in 1:nt) {
      plot(tau.keep[start:iter, t], type = "l",
           main = paste("tau = ", round(tau.t[t], 3)))
    }
  }
}

#### testing tau, phi, and beta ####
rm(list=ls())
source("../../../usefulR/usefulfunctions.R", chdir = TRUE)
source("auxfunctions.R", chdir = TRUE)
source("updatemodel.R", chdir = TRUE)

set.seed(2000)
ns <- 400
nt <- 12
np <- 6
X <- rX(ns, nt, np)
beta.t <- rnorm(np, 0, 10)

phi.t <- 0.2
s <- cbind(runif(ns), runif(ns))
d <- rdist(s)
Sigma.t <- exp(-d / phi.t)
tau.t   <- rgamma(nt, 1, 1)
Qb.t    <- chol2inv(chol(Sigma.t))
Xb.t    <- getXBeta(X = X, beta = beta.t)


mu <- matrix(0, ns, nt)
for (t in 1:nt) {
  mu[, t] <- Xb.t[, t] + t(chol(Sigma.t)) %*% rnorm(ns) / sqrt(tau.t[t])
}

# initialize values
beta <- rep(0, np)
Xb   <- getXBeta(X = X, beta = beta)
tau <- rep(1, nt)
phi <- 0.05
Qb <- chol2inv(chol(exp(-d / phi)))
SS <- getGPSS(Qb = Qb, param = mu, Xb = Xb)

niters <- 2000
burn   <- 1500
beta.sd <- 100
beta <- rep(0, np)
beta.keep <- matrix(0, niters, np)
beta.sd.keep <- rep(0, niters)
tau <- rep(1, nt)
tau.keep <- matrix(0, niters, nt)
phi.keep <- rep(0, niters)
acc.phi <- att.phi <- MH.phi <- 0.1
set.seed(3366)  # demo
for (iter in 1:niters) {
  this.update <- updateGPBeta(beta.sd = beta.sd, Qb = Qb,
                              param = mu, X = X, SS = SS, tau = tau)
  beta <- this.update$beta
  Xb   <- this.update$Xb
  SS   <- this.update$SS
  beta.keep[iter, ] <- beta

  this.update <- updateGPBetaSD(beta = beta, tau.a = 0.1, tau.b = 0.1)
  beta.sd <- this.update$beta.sd

  beta.sd.keep[iter] <- beta.sd

  this.update <- updateGPTau(SS = SS, tau.a = 0.1, tau.b = 0.1,
                             ns = ns)
  tau <- this.update$tau

  tau.keep[iter, ] <- tau

  this.update <- updateGPBW(bw = phi, bw.min = 0.01, bw.max = 1.2,
                            bw.mn = 0, bw.sd = 1, Qb = Qb, d = d,
                            mu = mu, Xb1 = Xb, tau1 = tau, SS1 = SS,
                            ls = mu, Xb2 = Xb, tau2 = tau, SS2 = SS,
                            acc = acc.phi, att = att.phi, MH = MH.phi)
  phi     <- this.update$bw
  Qb      <- this.update$Qb
  SS      <- this.update$SS1
  acc.phi <- this.update$acc
  att.phi <- this.update$att
  phi.keep[iter] <- phi

  if (iter < burn / 2) {
    this.update <- mhUpdate(acc = acc.phi, att = att.phi, MH = MH.phi)
    acc.phi <- this.update$acc
    att.phi <- this.update$att
    MH.phi  <- this.update$MH
  }

  if (iter %% 100 == 0) {
    par(mfrow = c(5, 3))
    for (i in 1:np) {
      plot(beta.keep[1:iter, i], type = "l",
           main = paste("beta: ", round(beta.t[i], 3)))
    }
    plot(beta.sd.keep[1:iter], type = "l", main = "beta sd")
    plot(phi.keep[1:iter], type = "l", main = paste("phi: ", phi.t))
    for(i in 1:7) {
      plot(tau.keep[1:iter, i], type = "l",
           main = paste("tau: ", round(tau.t[i], 3)))
    }
  }
}

#### testing mu ####
rm(list=ls())
source("../../../usefulR/usefulfunctions.R", chdir = TRUE)
source("auxfunctions.R", chdir = TRUE)
source("updatemodel.R", chdir = TRUE)

set.seed(2000)
ns <- 400
nt <- 12
np <- 6
X <- rX(ns, nt, np)
beta.t <- rnorm(np, 0, 10)

phi.t <- 0.2
s <- cbind(runif(ns), runif(ns))
d <- rdist(s)
Sigma.t <- exp(-d / phi.t)
tau.t   <- rgamma(nt, 1, 1)
Qb.t    <- chol2inv(chol(Sigma.t))
Xb.t    <- getXBeta(X = X, beta = beta.t)

mu.t <- matrix(0, ns, nt)
for (t in 1:nt) {
  mu.t[, t] <- Xb.t[, t] + t(chol(Sigma.t)) %*% rnorm(ns) / sqrt(tau.t[t])
}

ls.t <- matrix(0, ns, nt)
xi.t <- 0.01
y.t <- rgev(n = ns * nt, loc = mu.t, exp(ls.t), xi.t)

Sigma <- solve(Qb.t * tau.t[t])

# initialize values
mu <- matrix(mu.t + rnorm(ns * nt), ns, nt)
SS <- getGPSS(Qb = Qb.t, param = mu, Xb = Xb.t)
curll <- matrix(0, ns, nt)
for (t in 1:nt) {
  curll[, t] <- dgev(x = y.t[, t], loc = mu[, t], exp(ls.t[, t]), xi.t,
                     log = TRUE)
}

niters <- 10000
burn   <- 8000
mu.keep <- array(0, dim = c(niters, ns, nt))
acc.mu <- att.mu <- MH.mu <- matrix(0.2, ns, nt)

set.seed(3366)  # demo
for (iter in 1:niters) {
  this.update <- updateMuTest(mu = mu, Qb = Qb.t, tau = tau.t, Xb = Xb.t,
                              y = y.t, ls = ls.t, xi = xi.t,
                              SS = SS, curll = curll, acc = acc.mu,
                              att = att.mu, MH = MH.mu)
  mu <- this.update$mu
  SS <- this.update$SS
  curll <- this.update$curll
  acc.mu <- this.update$acc
  att.mu <- this.update$att
  mu.keep[iter, , ] <- mu

  if (iter < burn / 2) {
    this.update <- mhUpdate(acc = acc.mu, att = att.mu, MH = MH.mu,
                            target.min = 0.4, target.max = 0.7,
                            nattempts = 400)
    acc.mu <- this.update$acc
    att.mu <- this.update$att
    MH.mu  <- this.update$MH
  }

  if (iter %% 500 == 0) {
    par(mfrow = c(3, 3))
    start <- max(1, iter - 20000)
    for (i in 1:3) {
      for (j in 1:3) {
        plot(mu.keep[start:iter, i, j], type = "l",
             main = paste("mu: ", round(mu.t[i, j], 3)),
             ylab = round(acc.mu[i, j] / att.mu[i, j], 3),
             xlab = MH.mu[i, j])
      }
    }
  }
}


#### testing mu and tau ####
rm(list=ls())
source("../../../usefulR/usefulfunctions.R", chdir = TRUE)
source("auxfunctions.R", chdir = TRUE)
source("updatemodel.R", chdir = TRUE)

set.seed(2000)
ns <- 400
nt <- 12
np <- 6
X <- rX(ns, nt, np)
beta.t <- rnorm(np, 0, 10)

phi.t <- 0.2
s <- cbind(runif(ns), runif(ns))
d <- rdist(s)
Sigma.t <- exp(-d / phi.t)
tau.t   <- rgamma(nt, 1, 1)
Qb.t    <- chol2inv(chol(Sigma.t))
Xb.t    <- getXBeta(X = X, beta = beta.t)

mu.t <- matrix(0, ns, nt)
for (t in 1:nt) {
  mu.t[, t] <- Xb.t[, t] + t(chol(Sigma.t)) %*% rnorm(ns) / sqrt(tau.t[t])
}

ls.t <- matrix(0, ns, nt)
xi.t <- 0.01
y.t <- rgev(n = ns * nt, loc = mu.t, exp(ls.t), xi.t)

# initialize values
mu <- matrix(mu.t + rnorm(ns * nt), ns, nt)
SS <- getGPSS(Qb = Qb.t, param = mu, Xb = Xb.t)
curll <- matrix(0, ns, nt)
for (t in 1:nt) {
  curll[, t] <- dgev(x = y.t[, t], loc = mu[, t], exp(ls.t[, t]), xi.t,
                     log = TRUE)
}

niters <- 60000
burn   <- 50000
mu.keep <- array(0, dim = c(niters, ns, nt))
tau <- rep(1, nt)
tau.keep <- matrix(0, niters, nt)
acc.mu <- att.mu <- MH.mu <- matrix(0.1, ns, nt)

set.seed(3366)  # demo
for (iter in 1:niters) {
  this.update <- updateMuTest(mu = mu, Qb = Qb.t, tau = tau, Xb = Xb.t,
                              y = y.t, ls = ls.t, xi = xi.t,
                              SS = SS, curll = curll, acc = acc.mu,
                              att = att.mu, MH = MH.mu)
  mu <- this.update$mu
  SS <- this.update$SS
  curll <- this.update$curll
  acc.mu <- this.update$acc
  att.mu <- this.update$att
  mu.keep[iter, , ] <- mu

  this.update <- updateGPTau(SS = SS, tau.a = 0.1, tau.b = 0.1,
                             ns = ns)
  tau <- this.update$tau
  tau.keep[iter, ] <- tau

  if (iter < burn / 2) {
    this.update <- mhUpdate(acc = acc.mu, att = att.mu, MH = MH.mu,
                            target.min = 0.5, target.max = 0.7,
                            nattempts = 200)
    acc.mu <- this.update$acc
    att.mu <- this.update$att
    MH.mu  <- this.update$MH
  }

  if (iter %% 1000 == 0) {
    par(mfrow = c(4, 3))
    start <- max(1, iter - 20000)
    for (i in 1:3) {
      for (j in 1:3) {
        plot(mu.keep[start:iter, i, j], type = "l",
             main = paste("mu: ", round(mu.t[i, j], 3)),
             ylab = round(acc.mu[i, j] / att.mu[i, j], 3),
             xlab = MH.mu[i, j])
      }
    }
    for (i in 1:3) {
      plot(tau.keep[start:iter, i], type = "l",
           main = paste("tau: ", round(tau.t[i], 3)))
    }
  }
}

#### testing mu, tau, and beta ####
rm(list=ls())
source("../../../usefulR/usefulfunctions.R", chdir = TRUE)
source("auxfunctions.R", chdir = TRUE)
source("updatemodel.R", chdir = TRUE)

set.seed(2000)
ns <- 400
nt <- 12
np <- 6
X <- rX(ns, nt, np)
beta.t <- rnorm(np, 0, 10)

phi.t <- 0.2
s <- cbind(runif(ns), runif(ns))
d <- rdist(s)
Sigma.t <- exp(-d / phi.t)
tau.t   <- rgamma(nt, 1, 1)
Qb.t    <- chol2inv(chol(Sigma.t))
Xb.t    <- getXBeta(X = X, beta = beta.t)

mu.t <- matrix(0, ns, nt)
for (t in 1:nt) {
  mu.t[, t] <- Xb.t[, t] + t(chol(Sigma.t)) %*% rnorm(ns) / sqrt(tau.t[t])
}

ls.t <- matrix(0, ns, nt)
xi.t <- 0.01
y.t <- rgev(n = ns * nt, loc = mu.t, exp(ls.t), xi.t)

# initialize values
mu <- matrix(mu.t + rnorm(ns * nt), ns, nt)
SS <- getGPSS(Qb = Qb.t, param = mu, Xb = Xb.t)
curll <- matrix(0, ns, nt)
for (t in 1:nt) {
  curll[, t] <- dgev(x = y.t[, t], loc = mu[, t], exp(ls.t[, t]), xi.t,
                     log = TRUE)
}

niters <- 20000
burn   <- 15000
beta.sd <- 100
beta <- rep(0, np)
beta.keep <- matrix(0, niters, np)
beta.sd.keep <- rep(0, niters)
mu.keep <- array(0, dim = c(niters, ns, nt))
tau <- rep(1, nt)
tau.keep <- matrix(0, niters, nt)
acc.mu <- att.mu <- MH.mu <- matrix(0.1, ns, nt)

set.seed(3366)  # demo
for (iter in 1:niters) {
  this.update <- updateGPBeta(beta.sd = beta.sd, Qb = Qb.t,
                              param = mu, X = X, SS = SS, tau = tau)
  beta <- this.update$beta
  Xb   <- this.update$Xb
  SS   <- this.update$SS
  beta.keep[iter, ] <- beta

  this.update <- updateGPBetaSD(beta = beta, tau.a = 0.1, tau.b = 1)
  beta.sd <- this.update$beta.sd

  beta.sd.keep[iter] <- beta.sd

  this.update <- updateMuTest(mu = mu, Qb = Qb.t, tau = tau.t, Xb = Xb,
                              y = y.t, ls = ls.t, xi = xi.t,
                              SS = SS, curll = curll, acc = acc.mu,
                              att = att.mu, MH = MH.mu)
  mu <- this.update$mu
  SS <- this.update$SS
  curll <- this.update$curll
  acc.mu <- this.update$acc
  att.mu <- this.update$att
  mu.keep[iter, , ] <- mu

  this.update <- updateGPTau(SS = SS, tau.a = 0.1, tau.b = 0.1,
                             ns = ns)
  tau <- this.update$tau
  tau.keep[iter, ] <- tau

  if (iter < burn / 2) {
    this.update <- mhUpdate(acc = acc.mu, att = att.mu, MH = MH.mu,
                            target.min = 0.5, target.max = 0.7,
                            nattempts = 200)
    acc.mu <- this.update$acc
    att.mu <- this.update$att
    MH.mu  <- this.update$MH
  }

  if (iter %% 1000 == 0) {
    par(mfrow = c(4, 3))
    start <- max(1, iter - 20000)
    for (i in 1:2) {
      for (j in 1:3) {
        plot(mu.keep[start:iter, i, j], type = "l",
             main = paste("mu: ", round(mu.t[i, j], 3)),
             ylab = round(acc.mu[i, j] / att.mu[i, j], 3),
             xlab = MH.mu[i, j])
      }
    }

    for (i in 1:3) {
      plot(beta.keep[start:iter, i], type = "l",
           main = paste("beta: ", round(beta.t[i], 3)))
    }

    for (i in 1:3) {
      plot(tau.keep[start:iter, i], type = "l",
           main = paste("tau: ", round(tau.t[i], 3)))
    }
  }
}

#### Verify gradients - no residual dependence ####
rm(list=ls())
source("../../../usefulR/usefulfunctions.R", chdir = TRUE)
source("auxfunctions.R", chdir = TRUE)
source("updatemodel.R", chdir = TRUE)

set.seed(2000)
ns <- 10
nt <- 3
np <- 6
X1 <- rX(ns, nt, np)
X2 <- rX(ns, nt, np)
beta1.t <- rnorm(np, 0, 1)
beta2.t <- rnorm(np, 0, 0.1)

phi.t <- 0.2
s <- cbind(runif(ns), runif(ns))
d <- rdist(s)
Sigma.t <- exp(-d / phi.t)
tau.t   <- rgamma(nt, 1, 1)
Qb.t    <- chol2inv(chol(Sigma.t))

Xb1.t <- getXBeta(X = X1, beta = beta1.t)
Xb2.t <- getXBeta(X = X2, beta = beta2.t)

mu.t <- ls.t <- matrix(0, ns, nt)
for (t in 1:nt) {
  mu.t[, t] <- Xb1.t[, t] + t(chol(Sigma.t)) %*% rnorm(ns) / sqrt(tau.t[t])
  ls.t[, t] <- Xb2.t[, t] + t(chol(Sigma.t)) %*% rnorm(ns) / sqrt(tau.t[t])
}

xi.t <- 0.1
y.t <- rgev(n = ns * nt, loc = mu.t, scale = exp(ls.t), xi.t)

lp.mu <- logpost.mu.test(mu = mu.t[, t], Xb = Xb1.t[, t], tau = tau.t[t],
                         Qb = Qb.t, y = y.t[, t], ls = ls.t[, t], xi = xi.t)

mean(grad(func = logpost.mu.test, x = mu.t[, t], Xb = Xb1.t[, t],
          tau = tau.t[t], Qb = Qb.t, y = y.t[, t], ls = ls.t[, t], xi = xi.t) /
       logpost.mu.grad.test(mu = mu.t[, t], Xb = Xb1.t[, t], tau = tau.t[t],
                            Qb = Qb.t, y = y.t[, t], ls = ls.t[, t],
                            xi = xi.t))

sd(grad(func = logpost.mu.test, x = mu.t[, t], Xb = Xb1.t[, t], tau = tau.t[t],
        Qb = Qb.t, y = y.t[, t], ls = ls.t[, t], xi = xi.t) /
     logpost.mu.grad.test(mu = mu.t[, t], Xb = Xb1.t[, t], tau = tau.t[t],
                          Qb = Qb.t, y = y.t[, t], ls = ls.t[, t],
                          xi = xi.t))

lp.logsig <- logpost.logsig.test(ls = ls.t[, t], Xb = Xb2.t[, t],
                                 tau = tau.t[t], Qb = Qb.t, y = y.t[, t],
                                 mu = mu.t[, t], xi = xi.t)

mean(grad(func = logpost.logsig.test, x = ls.t[, t], Xb = Xb2.t[, t],
          tau = tau.t[t], Qb = Qb.t, y = y.t[, t], mu = mu.t[, t], xi = xi.t) /
       logpost.logsig.grad.test(ls = ls.t[, t], Xb = Xb2.t[, t],
                                tau = tau.t[t], Qb = Qb.t, y = y.t[, t],
                                mu = mu.t[, t], xi = xi.t))

sd(grad(func = logpost.logsig.test, x = ls.t[, t], Xb = Xb2.t[, t],
        tau = tau.t[t], Qb = Qb.t, y = y.t[, t], mu = mu.t[, t], xi = xi.t) /
     logpost.logsig.grad.test(ls = ls.t[, t], Xb = Xb2.t[, t],
                              tau = tau.t[t], Qb = Qb.t, y = y.t[, t],
                              mu = mu.t[, t], xi = xi.t))

#### testing logsig ####
rm(list=ls())
source("../../../usefulR/usefulfunctions.R", chdir = TRUE)
source("auxfunctions.R", chdir = TRUE)
source("updatemodel.R", chdir = TRUE)

set.seed(2000)
ns <- 400
nt <- 12
np <- 6
X1.t <- rX(ns, nt, np)
X2.t <- rX(ns, nt, np)
beta1.t <- rnorm(np, 0, 10)
beta2.t <- rnorm(np, 0, 5)

phi.t <- 0.2
s <- cbind(runif(ns), runif(ns))
d <- rdist(s)
Sigma.t <- exp(-d / phi.t)
tau.t   <- rgamma(nt, 1, 1)
Qb.t    <- chol2inv(chol(Sigma.t))
Xb1.t   <- getXBeta(X = X1.t, beta = beta1.t)
Xb2.t   <- getXBeta(X = X2.t, beta = beta2.t)

mu.t <- ls.t <- matrix(0, ns, nt)
for (t in 1:nt) {
  mu.t[, t] <- Xb1.t[, t] + t(chol(Sigma.t)) %*% rnorm(ns) / sqrt(tau.t[t])
  ls.t[, t] <- Xb2.t[, t] + t(chol(Sigma.t)) %*% rnorm(ns) / sqrt(tau.t[t])
}

xi.t <- 0.01
y.t <- rgev(n = ns * nt, loc = mu.t, exp(ls.t), xi.t)

# initialize values
ls <- matrix(ls.t + rnorm(ns * nt, 0, 0.1), ns, nt)
SS <- getGPSS(Qb = Qb.t, param = ls, Xb = Xb2.t)
curll <- matrix(0, ns, nt)
for (t in 1:nt) {
  curll[, t] <- dgev(x = y.t[, t], loc = mu.t[, t], exp(ls[, t]), xi.t,
                     log = TRUE)
}

niters <- 10000
burn   <- 8000
ls.keep <- array(0, dim = c(niters, ns, nt))
acc.ls <- att.ls <- MH.ls <- matrix(0.1, ns, nt)

set.seed(3366)  # demo
for (iter in 1:niters) {
  this.update <- updateLSTest(ls = ls, tau = tau.t, Xb = Xb2.t, SS = SS,
                              y = y.t, mu = mu.t, xi = xi.t,
                              Qb = Qb.t, curll = curll, acc = acc.ls,
                              att = att.ls, MH = MH.ls)
  ls <- this.update$ls
  SS <- this.update$SS
  curll <- this.update$curll
  acc.ls <- this.update$acc
  att.ls <- this.update$att
  ls.keep[iter, , ] <- ls

  if (iter < burn / 2) {
    this.update <- mhUpdate(acc = acc.ls, att = att.ls, MH = MH.ls,
                            target.min = 0.4, target.max = 0.7,
                            nattempts = 400)
    acc.ls <- this.update$acc
    att.ls <- this.update$att
    MH.ls  <- this.update$MH
  }

  if (iter %% 500 == 0) {
    par(mfrow = c(3, 3))
    start <- max(1, iter - 20000)
    for (i in 1:3) {
      for (j in 1:3) {
        plot(ls.keep[start:iter, i, j], type = "l",
             main = paste("logsig: ", round(ls.t[i, j], 3)),
             ylab = round(acc.ls[i, j] / att.ls[i, j], 3),
             xlab = MH.ls[i, j])
      }
    }
  }
}

#### testing logsig, tau, and beta ####
rm(list=ls())
source("../../../usefulR/usefulfunctions.R", chdir = TRUE)
source("auxfunctions.R", chdir = TRUE)
source("updatemodel.R", chdir = TRUE)

set.seed(2000)
ns <- 400
nt <- 12
np <- 6
X1.t <- rX(ns, nt, np)
X2.t <- rX(ns, nt, np)
beta1.t <- rnorm(np, 0, 10)
beta2.t <- rnorm(np, 0, 5)

phi.t <- 0.2
s <- cbind(runif(ns), runif(ns))
d <- rdist(s)
Sigma.t <- exp(-d / phi.t)
tau.t   <- rgamma(nt, 1, 1)
Qb.t    <- chol2inv(chol(Sigma.t))
Xb1.t   <- getXBeta(X = X1.t, beta = beta1.t)
Xb2.t   <- getXBeta(X = X2.t, beta = beta2.t)

mu.t <- ls.t <- matrix(0, ns, nt)
for (t in 1:nt) {
  mu.t[, t] <- Xb1.t[, t] + t(chol(Sigma.t)) %*% rnorm(ns) / sqrt(tau.t[t])
  ls.t[, t] <- Xb2.t[, t] + t(chol(Sigma.t)) %*% rnorm(ns) / sqrt(tau.t[t])
}

xi.t <- 0.01
y.t <- rgev(n = ns * nt, loc = mu.t, exp(ls.t), xi.t)

# initialize values
ls <- matrix(ls.t + rnorm(ns * nt, 0, 0.1), ns, nt)
SS <- getGPSS(Qb = Qb.t, param = ls, Xb = Xb2.t)
curll <- matrix(0, ns, nt)
for (t in 1:nt) {
  curll[, t] <- dgev(x = y.t[, t], loc = mu.t[, t], exp(ls[, t]), xi.t,
                     log = TRUE)
}

niters <- 30000
burn   <- 25000
beta.sd <- 100
beta <- rep(0, np)
beta.keep <- matrix(0, niters, np)
beta.sd.keep <- rep(0, niters)
ls.keep <- array(0, dim = c(niters, ns, nt))
tau <- rep(1, nt)
tau.keep <- matrix(0, niters, nt)
acc.ls <- att.ls <- MH.ls <- matrix(0.1, ns, nt)

set.seed(3366)  # demo
for (iter in 1:niters) {
  this.update <- updateGPBeta(beta.sd = beta.sd, Qb = Qb.t,
                              param = ls, X = X2.t, SS = SS, tau = tau)
  beta <- this.update$beta
  Xb2  <- this.update$Xb
  SS   <- this.update$SS
  beta.keep[iter, ] <- beta

  this.update <- updateGPBetaSD(beta = beta, tau.a = 0.1, tau.b = 1)
  beta.sd <- this.update$beta.sd

  beta.sd.keep[iter] <- beta.sd

  this.update <- updateLSTest(ls = ls, tau = tau.t, Xb = Xb2, SS = SS,
                              y = y.t, mu = mu.t, xi = xi.t,
                              Qb = Qb.t, curll = curll,
                              acc = acc.ls, att = att.ls, MH = MH.ls)
  ls    <- this.update$ls
  SS    <- this.update$SS
  curll <- this.update$curll
  acc.ls <- this.update$acc
  att.ls <- this.update$att
  ls.keep[iter, , ] <- ls

  this.update <- updateGPTau(SS = SS, tau.a = 0.1, tau.b = 0.1,
                             ns = ns)
  tau <- this.update$tau
  tau.keep[iter, ] <- tau

  if (iter < burn / 2) {
    this.update <- mhUpdate(acc = acc.ls, att = att.ls, MH = MH.ls,
                            target.min = 0.4, target.max = 0.7,
                            nattempts = 400)
    acc.ls <- this.update$acc
    att.ls <- this.update$att
    MH.ls  <- this.update$MH
  }

  if (iter %% 1000 == 0) {
    par(mfrow = c(4, 3))
    start <- max(1, iter - 20000)
    for (i in 1:2) {
      for (j in 1:3) {
        plot(ls.keep[start:iter, i, j], type = "l",
             main = paste("logsig: ", round(ls.t[i, j], 3)),
             ylab = round(acc.ls[i, j] / att.ls[i, j], 3),
             xlab = MH.ls[i, j])
      }
    }

    for (i in 1:3) {
      plot(beta.keep[start:iter, i], type = "l",
           main = paste("beta: ", round(beta2.t[i], 3)))
    }

    for (i in 1:3) {
      plot(tau.keep[start:iter, i], type = "l",
           main = paste("tau: ", round(tau.t[i], 3)))
    }
  }
}

#### testing basis bandwidth update ####
rm(list=ls())
source("../../../usefulR/usefulfunctions.R", chdir = TRUE)
source("auxfunctions.R", chdir = TRUE)
source("updatemodel.R", chdir = TRUE)

set.seed(2000)
# setting np later after X is created
ns <- 400
nt <- 12
nknots <- 5
time.interact <- TRUE

s <- cbind(runif(ns), runif(ns))
knots <- as.matrix(cover.design(R = s, nd = nknots)$design)
d <- rdist(s)
dw2 <- rdist(s, knots)^2
dw2[dw2 < 1e-4] <- 0

# create the matrix of covariates
X1.t <- X2.t <- array(1, dim = c(ns, nt, 2))
for (t in 1:nt) {
  time <- (t - nt / 2) / nt
  X1.t[, t, 2] <- X2.t[, t, 2] <- time
}

bw.basis.t <- 0.2
B.t <- makeW(dw2 = dw2, rho = bw.basis.t)
X1.t <- add.basis.X(X1.t, B.t, time.interact = time.interact)
X2.t <- add.basis.X(X2.t, B.t, time.interact = time.interact)

np <- dim(X1.t)[3]

beta1.t <- rnorm(np, 0, 10)
beta2.t <- rnorm(np, 0, 5)

bw.gp.t <- 0.2
Sigma.t <- exp(-d / bw.gp.t)
tau1.t <- tau2.t <- rgamma(nt, 1, 1)
Qb.t    <- chol2inv(chol(Sigma.t))

Xb1.t <- getXBeta(X = X1.t, beta = beta1.t)
Xb2.t <- getXBeta(X = X2.t, beta = beta2.t)

mu.t <- ls.t <- matrix(0, ns, nt)
for (t in 1:nt) {
  mu.t[, t] <- Xb1.t[, t] + t(chol(Sigma.t)) %*% rnorm(ns) / sqrt(tau1.t[t])
  ls.t[, t] <- Xb2.t[, t] + t(chol(Sigma.t)) %*% rnorm(ns) / sqrt(tau2.t[t])
}

xi.t <- 0.01
y.t <- rgev(n = ns * nt, loc = mu.t, scale = exp(ls.t), shape = xi.t)

# initialize values
bw.basis <- 0.4
bw.basis.min <- quantile(dw2, 0.01)
bw.basis.max <- quantile(dw2, 0.99)
B <- makeW(dw2 = dw2, rho = bw.basis)
X1 <- rep.basis.X(X = X1.t, newB = B, time.interact = time.interact)
X2 <- rep.basis.X(X = X2.t, newB = B, time.interact = time.interact)
Xb1 <- getXBeta(X = X1, beta = beta1.t)
Xb2 <- getXBeta(X = X2, beta = beta1.t)
SS1 <- getGPSS(Qb = Qb.t, param = mu.t, Xb = Xb1)
SS2 <- getGPSS(Qb = Qb.t, param = ls.t, Xb = Xb2)

niters <- 30000
burn   <- 25000

# storage
bw.basis.keep <- rep(0, niters)
Xb1.keep     <- array(0, dim = c(niters, ns, nt))
acc.bw.basis <- att.bw.basis <- MH.bw.basis <- 0.1

set.seed(3366)  # demo
for (iter in 1:niters) {
  this.update <- updateXBasisBW(bw = bw.basis, bw.min = bw.basis.min,
                                bw.max = bw.basis.max,
                                X1 = X1, beta1 = beta1.t, Xb1 = Xb1,
                                mu = mu.t, tau1 = tau1.t, SS1 = SS1,
                                X2 = X2, beta2 = beta2.t, Xb2 = Xb2,
                                ls = ls.t, tau2 = tau2.t, SS2 = SS2,
                                Qb = Qb.t, dw2 = dw2,
                                time.interact = time.interact,
                                acc = acc.bw.basis, att = att.bw.basis,
                                MH = MH.bw.basis)
  bw.basis <- this.update$bw
  X1  <- this.update$X1
  Xb1 <- this.update$Xb1
  SS1 <- this.update$SS1
  X2  <- this.update$X2
  Xb2 <- this.update$Xb2
  SS2 <- this.update$SS2
  acc.bw.basis <- this.update$acc
  att.bw.basis <- this.update$att

  bw.basis.keep[iter] <- bw.basis
  Xb1.keep[iter, , ] <- Xb1

  if (iter < burn / 2) {
    this.update <- mhUpdate(acc = acc.bw.basis, att = att.bw.basis,
                            MH = MH.bw.basis,
                            target.min = 0.3, target.max = 0.6,
                            lower = 0.8, higher = 1.2)
    acc.bw.basis <- this.update$acc
    att.bw.basis <- this.update$att
    MH.bw.basis  <- this.update$MH
  }

  if (iter %% 1000 == 0) {
    par(mfrow = c(2, 5))
    plot(bw.basis.keep[100:iter], type = "l",
         main = paste("BW basis: ", bw.basis.t, sep = ""),
         ylab = round(acc.bw.basis / att.bw.basis, 3),
         xlab = MH.bw.basis)
    for(i in 1:3) { for (j in 1:3) {
      plot(Xb1.keep[100:iter, i, j], type = "l",
           main = paste("Xb1: ", Xb1.t[i, j], sep = ""))
    }}
  }
}

#### testing basis bandwidth, beta1, beta2 ####
rm(list=ls())
source("../../../usefulR/usefulfunctions.R", chdir = TRUE)
source("auxfunctions.R", chdir = TRUE)
source("updatemodel.R", chdir = TRUE)

set.seed(2000)
# setting np later after X is created
ns <- 400
nt <- 12
nknots <- 5
time.interact <- TRUE

s <- cbind(runif(ns), runif(ns))
knots <- as.matrix(cover.design(R = s, nd = nknots)$design)
d <- rdist(s)
dw2 <- rdist(s, knots)^2
dw2[dw2 < 1e-4] <- 0

# create the matrix of covariates
X1.t <- X2.t <- array(1, dim = c(ns, nt, 2))
for (t in 1:nt) {
  time <- (t - nt / 2) / nt
  X1.t[, t, 2] <- X2.t[, t, 2] <- time
}

bw.basis.t <- 0.2
B.t <- makeW(dw2 = dw2, rho = bw.basis.t)
X1.t <- add.basis.X(X1.t, B.t, time.interact = time.interact)
X2.t <- add.basis.X(X2.t, B.t, time.interact = time.interact)

np <- dim(X1.t)[3]

beta1.t <- rnorm(np, 0, 10)
beta2.t <- rnorm(np, 0, 5)

bw.gp.t <- 0.2
Sigma.t <- exp(-d / bw.gp.t)
tau1.t <- tau2.t <- rgamma(nt, 1, 1)
Qb.t    <- chol2inv(chol(Sigma.t))

Xb1.t <- getXBeta(X = X1.t, beta = beta1.t)
Xb2.t <- getXBeta(X = X2.t, beta = beta2.t)

mu.t <- ls.t <- matrix(0, ns, nt)
for (t in 1:nt) {
  mu.t[, t] <- Xb1.t[, t] + t(chol(Sigma.t)) %*% rnorm(ns) / sqrt(tau1.t[t])
  ls.t[, t] <- Xb2.t[, t] + t(chol(Sigma.t)) %*% rnorm(ns) / sqrt(tau2.t[t])
}

xi.t <- 0.01
y.t <- rgev(n = ns * nt, loc = mu.t, scale = exp(ls.t), shape = xi.t)

# initialize values
beta1 <- beta2 <- rep(0, np)
beta1.sd <- beta2.sd <- 100
bw.basis <- 0.4
bw.basis.min <- quantile(dw2, 0.01)
bw.basis.max <- quantile(dw2, 0.99)
B <- makeW(dw2 = dw2, rho = bw.basis)
X1 <- rep.basis.X(X = X1.t, newB = B, time.interact = time.interact)
X2 <- rep.basis.X(X = X2.t, newB = B, time.interact = time.interact)
Xb1 <- getXBeta(X = X1, beta = beta1)
Xb2 <- getXBeta(X = X2, beta = beta1)
SS1 <- getGPSS(Qb = Qb.t, param = mu.t, Xb = Xb1)
SS2 <- getGPSS(Qb = Qb.t, param = ls.t, Xb = Xb2)

niters <- 30000
burn   <- 25000

# storage
bw.basis.keep <- rep(0, niters)
beta1.keep <- beta2.keep <- matrix(0, niters, np)
beta1.sd.keep <- beta2.sd.keep <- rep(0, niters)
acc.bw.basis <- att.bw.basis <- MH.bw.basis <- 0.1

set.seed(3366)  # demo
for (iter in 1:niters) {
  this.update <- updateGPBeta(beta.sd = beta1.sd, Qb = Qb.t,
                              param = mu.t, X = X1, SS = SS1, tau = tau1.t)
  beta1 <- this.update$beta
  Xb1   <- this.update$Xb
  SS1   <- this.update$SS
  beta1.keep[iter, ] <- beta1

  this.update <- updateGPBetaSD(beta = beta1, tau.a = 0.5, tau.b = 0.5)
  beta1.sd <- this.update$beta.sd
  beta1.sd.keep[iter] <- beta1.sd

  this.update <- updateGPBeta(beta.sd = beta2.sd, Qb = Qb.t,
                              param = ls.t, X = X2, SS = SS2, tau = tau2.t)
  beta2 <- this.update$beta
  Xb2   <- this.update$Xb
  SS2   <- this.update$SS
  beta2.keep[iter, ] <- beta2

  this.update <- updateGPBetaSD(beta = beta2, tau.a = 0.5, tau.b = 0.5)
  beta2.sd <- this.update$beta.sd
  beta2.sd.keep[iter] <- beta2.sd

  this.update <- updateXBasisBW(bw = bw.basis, bw.min = bw.basis.min,
                                bw.max = bw.basis.max,
                                X1 = X1, beta1 = beta1, Xb1 = Xb1,
                                mu = mu.t, tau1 = tau1.t, SS1 = SS1,
                                X2 = X2, beta2 = beta2, Xb2 = Xb2,
                                ls = ls.t, tau2 = tau2.t, SS2 = SS2,
                                Qb = Qb.t, dw2 = dw2,
                                time.interact = time.interact,
                                acc = acc.bw.basis, att = att.bw.basis,
                                MH = MH.bw.basis)
  bw.basis <- this.update$bw
  X1  <- this.update$X1
  Xb1 <- this.update$Xb1
  SS1 <- this.update$SS1
  X2  <- this.update$X2
  Xb2 <- this.update$Xb2
  SS2 <- this.update$SS2
  acc.bw.basis <- this.update$acc
  att.bw.basis <- this.update$att

  bw.basis.keep[iter] <- bw.basis

  if (iter < burn / 2) {
    this.update <- mhUpdate(acc = acc.bw.basis, att = att.bw.basis,
                            MH = MH.bw.basis,
                            target.min = 0.3, target.max = 0.6,
                            lower = 0.8, higher = 1.2)
    acc.bw.basis <- this.update$acc
    att.bw.basis <- this.update$att
    MH.bw.basis  <- this.update$MH
  }

  if (iter %% 500 == 0) {
    par(mfrow = c(3, 3))
    plot(bw.basis.keep[100:iter], type = "l",
         main = paste("BW basis: ", bw.basis.t, sep = ""),
         ylab = round(acc.bw.basis / att.bw.basis, 3),
         xlab = MH.bw.basis)
    for(i in 1:4) {
      plot(beta1.keep[100:iter, i], type = "l",
           main = paste("beta 1: ", beta1.t[i], sep = ""))
    }
    for (i in 1:4) {
      plot(beta2.keep[100:iter, i], type = "l",
           main = paste("beta 2: ", beta2.t[i], sep = ""))
    }

  }
}

#### testing basis bandwidth, beta1, beta2, mu, and logsig ####
rm(list=ls())
source("../../../usefulR/usefulfunctions.R", chdir = TRUE)
source("auxfunctions.R", chdir = TRUE)
source("updatemodel.R", chdir = TRUE)

set.seed(2000)
# setting np later after X is created
ns <- 400
nt <- 12
nknots <- 5
time.interact <- TRUE

s <- cbind(runif(ns), runif(ns))
knots <- as.matrix(cover.design(R = s, nd = nknots)$design)
d <- rdist(s)
dw2 <- rdist(s, knots)^2
dw2[dw2 < 1e-4] <- 0

# create the matrix of covariates
X1.t <- X2.t <- array(1, dim = c(ns, nt, 2))
for (t in 1:nt) {
  time <- (t - nt / 2) / nt
  X1.t[, t, 2] <- X2.t[, t, 2] <- time
}

bw.basis.t <- 0.2
B.t <- makeW(dw2 = dw2, rho = bw.basis.t)
X1.t <- add.basis.X(X1.t, B.t, time.interact = time.interact)
X2.t <- add.basis.X(X2.t, B.t, time.interact = time.interact)

np <- dim(X1.t)[3]

beta1.t <- rnorm(np, 0, 10)
beta2.t <- rnorm(np, 0, 5)

bw.gp.t <- 0.2
Sigma.t <- exp(-d / bw.gp.t)
tau1.t <- tau2.t <- rgamma(nt, 1, 1)
Qb.t    <- chol2inv(chol(Sigma.t))

Xb1.t <- getXBeta(X = X1.t, beta = beta1.t)
Xb2.t <- getXBeta(X = X2.t, beta = beta2.t)

mu.t <- ls.t <- matrix(0, ns, nt)
for (t in 1:nt) {
  mu.t[, t] <- Xb1.t[, t] + t(chol(Sigma.t)) %*% rnorm(ns) / sqrt(tau1.t[t])
  ls.t[, t] <- Xb2.t[, t] + t(chol(Sigma.t)) %*% rnorm(ns) / sqrt(tau2.t[t])
}

xi.t <- 0.01
y.t <- rgev(n = ns * nt, loc = mu.t, scale = exp(ls.t), shape = xi.t)

# initialize values
mu <- matrix(mu.t + rnorm(ns * nt), ns, nt)
ls <- matrix(ls.t + rnorm(ns * nt), ns, nt)
beta1 <- beta2 <- rep(0, np)
beta1.sd <- beta2.sd <- 100
bw.basis <- 0.4
bw.basis.min <- quantile(dw2, 0.01)
bw.basis.max <- quantile(dw2, 0.99)
B <- makeW(dw2 = dw2, rho = bw.basis)
X1 <- rep.basis.X(X = X1.t, newB = B, time.interact = time.interact)
X2 <- rep.basis.X(X = X2.t, newB = B, time.interact = time.interact)
Xb1 <- getXBeta(X = X1, beta = beta1)
Xb2 <- getXBeta(X = X2, beta = beta1)
SS1 <- getGPSS(Qb = Qb.t, param = mu.t, Xb = Xb1)
SS2 <- getGPSS(Qb = Qb.t, param = ls.t, Xb = Xb2)

curll <- matrix(0, ns, nt)
for (t in 1:nt) {
  curll[, t] <- dgev(x = y.t[, t], loc = mu[, t], exp(ls[, t]), xi.t,
                     log = TRUE)
}

niters <- 30000
burn   <- 25000

# storage
bw.basis.keep <- rep(0, niters)
beta1.keep <- beta2.keep <- matrix(0, niters, np)
beta1.sd.keep <- beta2.sd.keep <- rep(0, niters)
mu.keep <- array(0, dim = c(niters, ns, nt))
ls.keep <- array(0, dim = c(niters, ns, nt))
acc.bw.basis <- att.bw.basis <- MH.bw.basis <- 0.1
acc.mu <- att.mu <- MH.mu <- matrix(0.2, ns, nt)
acc.ls <- att.ls <- MH.ls <- matrix(0.1, ns, nt)

set.seed(3366)  # demo
for (iter in 1:niters) {
  this.update <- updateGPBeta(beta.sd = beta1.sd, Qb = Qb.t,
                              param = mu, X = X1, SS = SS1, tau = tau1.t)
  beta1 <- this.update$beta
  Xb1   <- this.update$Xb
  SS1   <- this.update$SS
  beta1.keep[iter, ] <- beta1

  this.update <- updateGPBetaSD(beta = beta1, tau.a = 0.5, tau.b = 0.5)
  beta1.sd <- this.update$beta.sd
  beta1.sd.keep[iter] <- beta1.sd

  this.update <- updateGPBeta(beta.sd = beta2.sd, Qb = Qb.t,
                              param = ls, X = X2, SS = SS2, tau = tau2.t)
  beta2 <- this.update$beta
  Xb2   <- this.update$Xb
  SS2   <- this.update$SS
  beta2.keep[iter, ] <- beta2

  this.update <- updateGPBetaSD(beta = beta2, tau.a = 0.5, tau.b = 0.5)
  beta2.sd <- this.update$beta.sd
  beta2.sd.keep[iter] <- beta2.sd

  this.update <- updateXBasisBW(bw = bw.basis, bw.min = bw.basis.min,
                                bw.max = bw.basis.max,
                                X1 = X1, beta1 = beta1, Xb1 = Xb1,
                                mu = mu, tau1 = tau1.t, SS1 = SS1,
                                X2 = X2, beta2 = beta2, Xb2 = Xb2,
                                ls = ls, tau2 = tau2.t, SS2 = SS2,
                                Qb = Qb.t, dw2 = dw2,
                                time.interact = time.interact,
                                acc = acc.bw.basis, att = att.bw.basis,
                                MH = MH.bw.basis)
  bw.basis <- this.update$bw
  X1  <- this.update$X1
  Xb1 <- this.update$Xb1
  SS1 <- this.update$SS1
  X2  <- this.update$X2
  Xb2 <- this.update$Xb2
  SS2 <- this.update$SS2
  acc.bw.basis <- this.update$acc
  att.bw.basis <- this.update$att
  bw.basis.keep[iter] <- bw.basis

  this.update <- updateMuTest(mu = mu, tau = tau1.t, Xb = Xb1, SS = SS1,
                              y = y.t, ls = ls, xi = xi.t,
                              Qb = Qb.t, curll = curll, acc = acc.mu,
                              att = att.mu, MH = MH.mu)
  mu    <- this.update$mu
  SS1   <- this.update$SS
  curll <- this.update$curll
  acc.mu <- this.update$acc
  att.mu <- this.update$att
  mu.keep[iter, , ] <- mu

  this.update <- updateLSTest(ls = ls, tau = tau2.t, Xb = Xb2, SS = SS2,
                              y = y.t, mu = mu, xi = xi.t,
                              Qb = Qb.t, curll = curll, acc = acc.ls,
                              att = att.ls, MH = MH.ls)
  ls     <- this.update$ls
  SS2    <- this.update$SS
  curll  <- this.update$curll
  acc.ls <- this.update$acc
  att.ls <- this.update$att
  ls.keep[iter, , ] <- ls

  if (iter < burn / 2) {
    this.update <- mhUpdate(acc = acc.bw.basis, att = att.bw.basis,
                            MH = MH.bw.basis,
                            target.min = 0.3, target.max = 0.6,
                            lower = 0.8, higher = 1.2)
    acc.bw.basis <- this.update$acc
    att.bw.basis <- this.update$att
    MH.bw.basis  <- this.update$MH

    this.update <- mhUpdate(acc = acc.mu, att = att.mu, MH = MH.mu,
                            target.min = 0.4, target.max = 0.7,
                            nattempts = 400)
    acc.mu <- this.update$acc
    att.mu <- this.update$att
    MH.mu  <- this.update$MH

    this.update <- mhUpdate(acc = acc.ls, att = att.ls, MH = MH.ls,
                            target.min = 0.4, target.max = 0.7,
                            nattempts = 400)
    acc.ls <- this.update$acc
    att.ls <- this.update$att
    MH.ls  <- this.update$MH
  }

  if (iter %% 500 == 0) {
    par(mfrow = c(3, 3))
    plot(bw.basis.keep[100:iter], type = "l",
         main = paste("BW basis: ", bw.basis.t, sep = ""),
         ylab = round(acc.bw.basis / att.bw.basis, 3),
         xlab = MH.bw.basis)
    for(i in 1:4) {
      plot(beta1.keep[100:iter, i], type = "l",
           main = paste("beta 1: ", beta1.t[i], sep = ""))
    }
    for (i in 1:4) {
      plot(beta2.keep[100:iter, i], type = "l",
           main = paste("beta 2: ", beta2.t[i], sep = ""))
    }

  }
}


#### Verify gradients - with residual dependence ####
rm(list=ls())
source("../../../usefulR/usefulfunctions.R", chdir = TRUE)
source("auxfunctions.R", chdir = TRUE)
source("updatemodel.R", chdir = TRUE)

set.seed(2000)
ns <- 10
nt <- 3
np <- 6
X1 <- rX(ns, nt, np)
X2 <- rX(ns, nt, np)
beta1.t <- rnorm(np, 0, 1)
beta2.t <- rnorm(np, 0, 0.1)

phi.t <- 0.2
s <- cbind(runif(ns), runif(ns))
d <- rdist(s)
Sigma.t <- exp(-d / phi.t)
tau.t   <- rgamma(nt, 1, 1)
Qb.t    <- chol2inv(chol(Sigma.t))

Xb1.t <- getXBeta(X = X1, beta = beta1.t)
Xb2.t <- getXBeta(X = X2, beta = beta2.t)

mu.t <- ls.t <- matrix(0, ns, nt)
for (t in 1:nt) {
  mu.t[, t] <- Xb1.t[, t] + t(chol(Sigma.t)) %*% rnorm(ns) / sqrt(tau.t[t])
  ls.t[, t] <- Xb2.t[, t] + t(chol(Sigma.t)) %*% rnorm(ns) / sqrt(tau.t[t])
}

xi.t <- 0.1
y.t <- rgev(n = ns * nt, loc = mu.t, scale = exp(ls.t), xi.t)

nknots <- 4
theta.t <- matrix(abs(rnorm(ns * nt)), ns, nt)
alpha.t <- 0.4
thresh.t <- matrix(median(y.t), ns, nt)

xi.t <- 0
lp.mu <- logpost.mu(mu = mu.t[, t], Xb = Xb1.t[, t], tau = tau.t[t],
                    Qb = Qb.t, y = y.t[, t], ls = ls.t[, t], xi = xi.t,
                    theta = theta.t[, t], thresh = thresh.t[, t],
                    alpha = alpha.t)

mean(grad(func = logpost.mu, x = mu.t[, t], Xb = Xb1.t[, t], tau = tau.t[t],
          Qb = Qb.t, y = y.t[, t], ls = ls.t[, t], xi = xi.t,
          theta = theta.t[, t], thresh = thresh.t[, t], alpha = alpha.t) /
       logpost.mu.grad(mu = mu.t[, t], Xb = Xb1.t[, t], tau = tau.t[t],
                       Qb = Qb.t, y = y.t[, t], ls = ls.t[, t],
                       xi = xi.t, theta = theta.t[, t], thresh = thresh.t[, t],
                       alpha = alpha.t))

sd(grad(func = logpost.mu, x = mu.t[, t], Xb = Xb1.t[, t], tau = tau.t[t],
        Qb = Qb.t, y = y.t[, t], ls = ls.t[, t], xi = xi.t,
        theta = theta.t[, t], thresh = thresh.t[, t], alpha = alpha.t) /
     logpost.mu.grad(mu = mu.t[, t], Xb = Xb1.t[, t], tau = tau.t[t],
                     Qb = Qb.t, y = y.t[, t], ls = ls.t[, t],
                     xi = xi.t, theta = theta.t[, t], thresh = thresh.t[, t],
                     alpha = alpha.t))

lp.logsig <- logpost.logsig(ls = ls.t[, t], Xb = Xb2.t[, t], tau = tau.t[t],
                            Qb = Qb.t, y = y.t[, t], mu = mu.t[, t],
                            xi = xi.t, theta = theta.t[, t],
                            thresh = thresh.t[, t], alpha = alpha.t)

mean(grad(func = logpost.logsig, x = ls.t[, t], Xb = Xb2.t[, t],
          tau = tau.t[t], Qb = Qb.t, y = y.t[, t], mu = mu.t[, t], xi = xi.t,
          theta = theta.t[, t], thresh = thresh.t[, t], alpha = alpha.t) /
       logpost.logsig.grad(ls = ls.t[, t], Xb = Xb2.t[, t],
                           tau = tau.t[t], Qb = Qb.t, y = y.t[, t],
                           mu = mu.t[, t], xi = xi.t,
                           theta = theta.t[, t], thresh = thresh.t[, t],
                           alpha = alpha.t))

sd(grad(func = logpost.logsig, x = ls.t[, t], Xb = Xb2.t[, t],
        tau = tau.t[t], Qb = Qb.t, y = y.t[, t], mu = mu.t[, t], xi = xi.t,
        theta = theta.t[, t], thresh = thresh.t[, t], alpha = alpha.t) /
     logpost.logsig.grad(ls = ls.t[, t], Xb = Xb2.t[, t],
                         tau = tau.t[t], Qb = Qb.t, y = y.t[, t],
                         mu = mu.t[, t], xi = xi.t,
                         theta = theta.t[, t], thresh = thresh.t[, t],
                         alpha = alpha.t))


#### testing xi ####
rm(list=ls())
source("../../../usefulR/usefulfunctions.R", chdir = TRUE)
source("auxfunctions.R", chdir = TRUE)
source("updatemodel.R", chdir = TRUE)

set.seed(2000)
ns <- 400
nt <- 12
np <- 6
X1.t <- rX(ns, nt, np)
X2.t <- rX(ns, nt, np)
beta1.t <- rnorm(np, 0, 10)
beta2.t <- rnorm(np, 0, 5)

phi.t <- 0.2
s <- cbind(runif(ns), runif(ns))
d <- rdist(s)
Sigma.t <- exp(-d / phi.t)
tau.t   <- rgamma(nt, 1, 1)
Qb.t    <- chol2inv(chol(Sigma.t))
Xb1.t   <- getXBeta(X = X1.t, beta = beta1.t)
Xb2.t   <- getXBeta(X = X2.t, beta = beta2.t)

mu.t <- ls.t <- matrix(0, ns, nt)
for (t in 1:nt) {
  mu.t[, t] <- Xb1.t[, t] + t(chol(Sigma.t)) %*% rnorm(ns) / sqrt(tau.t[t])
  ls.t[, t] <- Xb2.t[, t] + t(chol(Sigma.t)) %*% rnorm(ns) / sqrt(tau.t[t])
}

xi.t <- -0.7
y.t <- rgev(n = ns * nt, loc = mu.t, exp(ls.t), xi.t)

# initialize values
xi <- 0.1

thresh.t <- matrix(-Inf, ns, nt)
theta.t <- matrix(1, ns, nt)
alpha.t <- 1
curll <- loglike(y = y.t, mu = mu.t, ls = ls.t, xi = xi,
                 theta = theta.t, thresh = thresh.t, alpha = alpha.t)


niters <- 10000
burn   <- 8000
xi.keep <- rep(0, niters)
acc.xi <- att.xi <- MH.xi <- 0.01


set.seed(3366)  # demo
for (iter in 1:niters) {
  this.update <- updateXi(xi = xi, xi.min = -2, xi.max = 2,
                          xi.mn = 0, xi.sd = 0.5, y = y.t, mu = mu.t,
                          ls = ls.t, curll = curll, theta = theta.t,
                          thresh = thresh.t, alpha = alpha.t, acc = acc.xi,
                          att = att.xi, MH = MH.xi)
  xi <- this.update$xi
  curll <- this.update$curll
  acc.xi <- this.update$acc
  att.xi <- this.update$att
  xi.keep[iter] <- xi

  if (iter < burn / 2) {
    this.update <- mhUpdate(acc = acc.xi, att = att.xi, MH = MH.xi,
                            target.min = 0.3, target.max = 0.6,
                            nattempts = 400)
    acc.xi <- this.update$acc
    att.xi <- this.update$att
    MH.xi  <- this.update$MH
  }

  if (iter %% 500 == 0) {
    start <- max(1, iter - 20000)
    plot(xi.keep[start:iter], type = "l", main = paste("xi: ", xi.t),
         ylab = round(acc.xi / att.xi, 3),
         xlab = MH.xi)
  }
}

#### testing beta1 ####
rm(list=ls())
source("../../../usefulR/usefulfunctions.R", chdir = TRUE)
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
acc.beta1  <- att.beta1 <- MH.beta1 <- matrix(0.1, ns, nt)

set.seed(3366)  # demo
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
                            target.min = 0.3, target.max = 0.6,
                            nattempts = 400)
    acc.beta1 <- this.update$acc
    att.beta1 <- this.update$att
    MH.beta1  <- this.update$MH
  }

  if (iter %% 500 == 0) {
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