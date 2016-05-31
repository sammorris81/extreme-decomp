rm(list=ls())
library(fields)
library(Rcpp)
library(emulator)
library(microbenchmark)
library(SpatialExtremes)

source("../../../usefulR/usefulfunctions.R", chdir = TRUE)
source("auxfunctions.R", chdir = TRUE)
source("updatemodel.R", chdir = TRUE)

# testing to see that parameter updates are working properly
set.seed(2000)
ns <- 400
nt <- 2
np <- 6
X <- array(rnorm(ns * nt * np), dim = c(ns, nt, np))
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
  this.update <- updateGPBeta(beta = beta, beta.sd = 100, Qb = Qb,
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

rm(list=ls())
source("../../../usefulR/usefulfunctions.R", chdir = TRUE)
source("auxfunctions.R", chdir = TRUE)
source("updatemodel.R", chdir = TRUE)

# testing to see that parameter updates are working properly
set.seed(2000)
ns <- 400
nt <- 10
np <- 6
X <- array(rnorm(ns * nt * np), dim = c(ns, nt, np))
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
                            logsig = mu, Xb2 = Xb, tau2 = tau, SS2 = SS,
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


rm(list=ls())
source("../../../usefulR/usefulfunctions.R", chdir = TRUE)
source("auxfunctions.R", chdir = TRUE)
source("updatemodel.R", chdir = TRUE)

# testing to see that parameter updates are working properly
set.seed(2000)
ns <- 400
nt <- 12
np <- 6
X <- array(rnorm(ns * nt * np), dim = c(ns, nt, np))
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


# test tau, phi, and beta
rm(list=ls())
source("../../../usefulR/usefulfunctions.R", chdir = TRUE)
source("auxfunctions.R", chdir = TRUE)
source("updatemodel.R", chdir = TRUE)

# testing to see that parameter updates are working properly
set.seed(2000)
ns <- 400
nt <- 12
np <- 6
X <- array(rnorm(ns * nt * np), dim = c(ns, nt, np))
beta.t <- rnorm(np, 0, 10)

phi.t <- 0.2
s <- cbind(runif(ns), runif(ns))
d <- rdist(s)
Sigma.t <- exp(-d / phi.t)
tau.t   <- rgamma(nt, 0.5, 0.5)
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


niters <- 10000
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
  this.update <- updateGPBeta(beta = beta, beta.sd = beta.sd, Qb = Qb,
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
                            logsig = mu, Xb2 = Xb, tau2 = tau, SS2 = SS,
                            acc = acc.phi, att = att.phi, MH = MH.phi)
  phi     <- this.update$bw
  Qb      <- this.update$Qb
  SS      <- this.update$SS1
  acc.phi <- this.update$acc
  att.phi <- this.update$att

  this.update <- mhUpdate(acc = acc.phi, att = att.phi, MH = MH.phi)
  acc.phi <- this.update$acc
  att.phi <- this.update$att
  MH.phi  <- this.update$MH

  phi.keep[iter] <- phi

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

# test mu
rm(list=ls())
source("../../../usefulR/usefulfunctions.R", chdir = TRUE)
source("auxfunctions.R", chdir = TRUE)
source("updatemodel.R", chdir = TRUE)

# testing to see that parameter updates are working properly
set.seed(2000)
ns <- 400
nt <- 12
np <- 6
X <- array(rnorm(ns * nt * np), dim = c(ns, nt, np))
beta.t <- rnorm(np, 0, 10)

phi.t <- 0.2
s <- cbind(runif(ns), runif(ns))
d <- rdist(s)
Sigma.t <- exp(-d / phi.t)
tau.t   <- rgamma(nt, 0.5, 0.5)
Qb.t    <- chol2inv(chol(Sigma.t))
Xb.t    <- getXBeta(X = X, beta = beta.t)

mu.t <- matrix(0, ns, nt)
for (t in 1:nt) {
  mu.t[, t] <- Xb.t[, t] + t(chol(Sigma.t)) %*% rnorm(ns) / sqrt(tau.t[t])
}

y.t <- rgev(n = ns * nt, loc = mu.t, 1, 0.1)

# initialize values
mu <- matrix(mu.t + rnorm(ns * nt), ns, nt)
SS <- getGPSS(Qb = Qb.t, param = mu, Xb = Xb.t)
curll <- matrix(0, ns, nt)
for (t in 1:nt) {
  curll[, t] <- dgev(x = y.t[, t], loc = mu[, t], 1, 0.1, log = TRUE)
}

niters <- 10000
mu.keep <- array(0, dim = c(niters, ns, nt))
acc.mu <- att.mu <- MH.mu <- matrix(0.5, ns, nt)

set.seed(3366)  # demo
for (iter in 1:niters) {
  this.update <- updateMuTest(mu = mu, Qb = Qb.t, tau = tau.t, Xb = Xb.t,
                              y = y.t, SS = SS, curll = curll, acc = acc.mu,
                              att = att.mu, MH = MH.mu)
  mu <- this.update$mu
  SS <- this.update$SS
  curll <- this.update$curll
  acc.mu <- this.update$acc
  att.mu <- this.update$att
  mu.keep[iter, , ] <- mu

  this.update <- mhUpdate(acc = acc.mu, att = att.mu, MH = MH.mu)
  acc.mu <- this.update$acc
  att.mu <- this.update$att
  MH.mu  <- this.update$MH

  if (iter %% 100 == 0) {
    par(mfrow = c(5, 3))
    for (i in 1:5) {
      for (j in 1:3) {
        plot(mu.keep[1:iter, i, j], type = "l",
             main = paste("mu: ", round(mu.t[i, j], 3)),
             ylab = round(acc.mu[i, j] / att.mu[i, j], 3),
             xlab = MH.mu[i, j])
      }
    }
  }
}




























# test tau, phi, and beta
rm(list=ls())
source("../../../usefulR/usefulfunctions.R", chdir = TRUE)
source("auxfunctions.R", chdir = TRUE)
source("updatemodel.R", chdir = TRUE)

# testing to see that parameter updates are working properly
set.seed(2000)
ns <- 400
nt <- 12
np <- 6
X <- array(rnorm(ns * nt * np), dim = c(ns, nt, np))
beta.t <- rnorm(np, 0, 10)

phi.t <- 0.2
s <- cbind(runif(ns), runif(ns))
d <- rdist(s)
Sigma.t <- exp(-d / phi.t)
tau.t   <- rgamma(nt, 0.5, 0.5)
Qb.t    <- chol2inv(chol(Sigma.t))
Xb      <- getXBeta(X = X, beta = beta.t)

if (nt == 1) {
  Xb <- matrix(Xb, ns, nt)
}

mu.t <- matrix(0, ns, nt)
for (t in 1:nt) {
  mu.t[, t] <- Xb[, t] + t(chol(Sigma.t)) %*% rnorm(ns) / sqrt(tau.t[t])
}

y <- rgev(n = ns * nt, loc = mu.t, 1, 0.1)

# initialize values
beta <- rep(0, np)
Xb   <- getXBeta(X = X, beta = beta)
mu <- matrix(mu.t + rnorm(ns * nt), ns, nt)
tau <- rep(1, nt)
phi <- 0.05
Qb <- chol2inv(chol(exp(-d / phi)))
SS <- getGPSS(Qb = Qb, param = mu, Xb = Xb)
curll <- matrix(0, ns, nt)
for (t in 1:nt) {
  curll[, t] <- dgev(x = y[, t], loc = mu[, t], 1, 0.1, log = TRUE)
}

niters <- 10000
mu.keep <- array(0, dim = c(niters, ns, nt))
acc.mu <- att.mu <- MH.mu <- matrix(0.5, ns, nt)
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
  this.update <- updateMuTest(mu = mu, Qb = Qb, tau = tau, Xb = Xb,
                              y = y, SS = SS, curll = curll, acc = acc.mu,
                              att = att.mu, MH = MH.mu)
  mu <- this.update$mu
  SS <- this.update$SS
  curll <- this.update$curll
  acc.mu <- this.update$acc
  att.mu <- this.update$att
  mu.keep[iter, , ] <- mu

  this.update <- mhUpdate(acc = acc.mu, att = att.mu, MH = MH.mu)
  acc.mu <- this.update$acc
  att.mu <- this.update$att
  MH.mu  <- this.update$MH

  this.update <- updateGPBeta(beta = beta, beta.sd = beta.sd, Qb = Qb,
                              param = mu, X = X, SS = SS, tau = tau)
  beta <- this.update$beta
  Xb   <- this.update$Xb
  SS   <- this.update$SS
  beta.keep[iter, ] <- beta

  this.update <- updateGPBetaSD(beta = beta, tau.a = 0.1, tau.b = 0.1)
  beta.sd <- this.update$beta.sd

  beta.sd.keep[iter] <- beta.sd

  this.update <- updateGPTau(tau = tau, SS = SS, tau.a = 0.1, tau.b = 0.1,
                             ns = ns)
  tau <- this.update$tau

  tau.keep[iter, ] <- tau

  this.update <- updateGPBW(bw = phi, bw.min = 0.01, bw.max = 1.2,
                            bw.mn = 0, bw.sd = 1, Qb = Qb, d = d,
                            mu = mu, Xb1 = Xb, tau1 = tau, SS1 = SS,
                            logsig = mu, Xb2 = Xb, tau2 = tau, SS2 = SS,
                            acc = acc.phi, att = att.phi, MH = MH.phi)
  phi     <- this.update$bw
  Qb      <- this.update$Qb
  SS      <- this.update$SS1
  acc.phi <- this.update$acc
  att.phi <- this.update$att

  this.update <- mhUpdate(acc = acc.phi, att = att.phi, MH = MH.phi)
  acc.phi <- this.update$acc
  att.phi <- this.update$att
  MH.phi  <- this.update$MH

  phi.keep[iter] <- phi

  if (iter %% 100 == 0) {
    par(mfrow = c(5, 3))
    for (i in 1:3) {
      plot(beta.keep[1:iter, i], type = "l",
           main = paste("beta: ", round(beta.t[i], 3)))
    }
    for (i in 1:3) {
      plot(mu.keep[1:iter, i, i], type = "l",
           main = paste("mu: ", round(mu.t[i, i], 3)),
           ylab = round(acc.mu[i, i] / att.mu[i, i], 3),
           xlab = MH.mu[i, i])
    }
    plot(beta.sd.keep[1:iter], type = "l", main = "beta sd")
    plot(phi.keep[1:iter], type = "l", main = paste("phi: ", phi.t))
    for(i in 1:7) {
      plot(tau.keep[1:iter, i], type = "l",
           main = paste("tau: ", round(tau.t[i], 3)))
    }
  }
}