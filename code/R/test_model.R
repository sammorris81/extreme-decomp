rm(list=ls())
library(fields)
library(Rcpp)
library(emulator)
library(microbenchmark)
library(SpatialExtremes)
library(numDeriv)

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
                            logsig = mu, Xb2 = Xb, tau2 = tau, SS2 = SS,
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
tau.t   <- rgamma(nt, 0.5, 0.5)
Qb.t    <- chol2inv(chol(Sigma.t))
Xb.t    <- getXBeta(X = X, beta = beta.t)

mu.t <- matrix(0, ns, nt)
for (t in 1:nt) {
  mu.t[, t] <- Xb.t[, t] + t(chol(Sigma.t)) %*% rnorm(ns) / sqrt(tau.t[t])
}

logsig.t <- matrix(0, ns, nt)
xi.t <- 0.01
y.t <- rgev(n = ns * nt, loc = mu.t, exp(logsig.t), xi.t)

Sigma <- solve(Qb.t * tau.t[t])

# initialize values
mu <- matrix(mu.t + rnorm(ns * nt), ns, nt)
SS <- getGPSS(Qb = Qb.t, param = mu, Xb = Xb.t)
curll <- matrix(0, ns, nt)
for (t in 1:nt) {
  curll[, t] <- dgev(x = y.t[, t], loc = mu[, t], exp(logsig.t[, t]), xi.t,
                     log = TRUE)
}

niters <- 10000
burn   <- 8000
mu.keep <- array(0, dim = c(niters, ns, nt))
acc.mu <- att.mu <- MH.mu <- matrix(0.2, ns, nt)

set.seed(3366)  # demo
for (iter in 1:niters) {
  this.update <- updateMuTest(mu = mu, Qb = Qb.t, tau = tau.t, Xb = Xb.t,
                              y = y.t, logsig = logsig.t, xi = xi.t,
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
tau.t   <- rgamma(nt, 0.5, 0.5)
Qb.t    <- chol2inv(chol(Sigma.t))
Xb.t    <- getXBeta(X = X, beta = beta.t)

mu.t <- matrix(0, ns, nt)
for (t in 1:nt) {
  mu.t[, t] <- Xb.t[, t] + t(chol(Sigma.t)) %*% rnorm(ns) / sqrt(tau.t[t])
}

logsig.t <- matrix(0, ns, nt)
xi.t <- 0.01
y.t <- rgev(n = ns * nt, loc = mu.t, exp(logsig.t), xi.t)

# initialize values
mu <- matrix(mu.t + rnorm(ns * nt), ns, nt)
SS <- getGPSS(Qb = Qb.t, param = mu, Xb = Xb.t)
curll <- matrix(0, ns, nt)
for (t in 1:nt) {
  curll[, t] <- dgev(x = y.t[, t], loc = mu[, t], exp(logsig.t[, t]), xi.t,
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
                              y = y.t, logsig = logsig.t, xi = xi.t,
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
tau.t   <- rgamma(nt, 0.5, 0.5)
Qb.t    <- chol2inv(chol(Sigma.t))
Xb.t    <- getXBeta(X = X, beta = beta.t)

mu.t <- matrix(0, ns, nt)
for (t in 1:nt) {
  mu.t[, t] <- Xb.t[, t] + t(chol(Sigma.t)) %*% rnorm(ns) / sqrt(tau.t[t])
}

logsig.t <- matrix(0, ns, nt)
xi.t <- 0.01
y.t <- rgev(n = ns * nt, loc = mu.t, exp(logsig.t), xi.t)

# initialize values
mu <- matrix(mu.t + rnorm(ns * nt), ns, nt)
SS <- getGPSS(Qb = Qb.t, param = mu, Xb = Xb.t)
curll <- matrix(0, ns, nt)
for (t in 1:nt) {
  curll[, t] <- dgev(x = y.t[, t], loc = mu[, t], exp(logsig.t[, t]), xi.t,
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
                              y = y.t, logsig = logsig.t, xi = xi.t,
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

#### Verify gradients ####
rm(list=ls())
source("../../../usefulR/usefulfunctions.R", chdir = TRUE)
source("auxfunctions.R", chdir = TRUE)
source("updatemodel.R", chdir = TRUE)

set.seed(2000)
ns <- 10
nt <- 3
np <- 6
X.mu <- rX(ns, nt, np)
X.logsig <- rX(ns, nt, np)
beta.mu.t <- rnorm(np, 0, 1)
beta.logsig.t <- rnorm(np, 0, 0.1)

phi.t <- 0.2
s <- cbind(runif(ns), runif(ns))
d <- rdist(s)
Sigma.t <- exp(-d / phi.t)
tau.t   <- rgamma(nt, 1, 1)
Qb.t    <- chol2inv(chol(Sigma.t))

Xb.mu.t <- getXBeta(X = X.mu, beta = beta.mu.t)
Xb.logsig.t <- getXBeta(X = X.logsig, beta = beta.logsig.t)

mu.t <- logsig.t <- matrix(0, ns, nt)
for (t in 1:nt) {
  mu.t[, t] <- Xb.mu.t[, t] + t(chol(Sigma.t)) %*% rnorm(ns) / sqrt(tau.t[t])
  logsig.t[, t] <- Xb.logsig.t[, t] + t(chol(Sigma.t)) %*% rnorm(ns) / sqrt(tau.t[t])
}

xi.t <- 0.1
y.t <- rgev(n = ns * nt, loc = mu.t, scale = exp(logsig.t), xi.t)

lp.mu <- logpost.mu(mu = mu.t[, t], Xb = Xb.mu.t[, t], tau = tau.t[t],
                    Qb = Qb.t, y = y.t[, t], logsig = logsig.t[, t], xi = xi.t)

mean(grad(func = logpost.mu, x = mu.t[, t], Xb = Xb.mu.t[, t], tau = tau.t[t],
          Qb = Qb.t, y = y.t[, t], logsig = logsig.t[, t], xi = xi.t) /
       logpost.mu.grad(mu = mu.t[, t], Xb = Xb.mu.t[, t], tau = tau.t[t],
                       Qb = Qb.t, y = y.t[, t], logsig = logsig.t[, t],
                       xi = xi.t))

sd(grad(func = logpost.mu, x = mu.t[, t], Xb = Xb.mu.t[, t], tau = tau.t[t],
        Qb = Qb.t, y = y.t[, t], logsig = logsig.t[, t], xi = xi.t) /
     logpost.mu.grad(mu = mu.t[, t], Xb = Xb.mu.t[, t], tau = tau.t[t],
                     Qb = Qb.t, y = y.t[, t], logsig = logsig.t[, t],
                     xi = xi.t))

lp.logsig <- logpost.logsig(mu = mu.t[, t], Xb = Xb.mu.t[, t], tau = tau.t[t],
                            Qb = Qb.t, y = y.t[, t], logsig = logsig.t[, t],
                            xi = xi.t)

mean(grad(func = logpost.logsig, x = logsig.t[, t], Xb = Xb.logsig.t[, t],
          mu = mu.t[, t], tau = tau.t[t], Qb = Qb.t, y = y.t[, t], xi = xi.t) /
       logpost.logsig.grad(mu = mu.t[, t], Xb = Xb.logsig.t[, t],
                           tau = tau.t[t], Qb = Qb.t, y = y.t[, t],
                           logsig = logsig.t[, t], xi = xi.t))

sd(grad(func = logpost.logsig, x = logsig.t[, t], Xb = Xb.logsig.t[, t],
        mu = mu.t[, t], tau = tau.t[t], Qb = Qb.t, y = y.t[, t], xi = xi.t) /
     logpost.logsig.grad(mu = mu.t[, t], Xb = Xb.logsig.t[, t],
                         tau = tau.t[t], Qb = Qb.t, y = y.t[, t],
                         logsig = logsig.t[, t], xi = xi.t))

#### testing logsig ####
rm(list=ls())
source("../../../usefulR/usefulfunctions.R", chdir = TRUE)
source("auxfunctions.R", chdir = TRUE)
source("updatemodel.R", chdir = TRUE)

set.seed(2000)
ns <- 400
nt <- 12
np <- 6
X.mu.t <- rX(ns, nt, np)
X.logsig.t <- rX(ns, nt, np)
beta.mu.t <- rnorm(np, 0, 10)
beta.logsig.t <- rnorm(np, 0, 5)

phi.t <- 0.2
s <- cbind(runif(ns), runif(ns))
d <- rdist(s)
Sigma.t <- exp(-d / phi.t)
tau.t   <- rgamma(nt, 1, 1)
Qb.t    <- chol2inv(chol(Sigma.t))
Xb.mu.t    <- getXBeta(X = X.mu.t, beta = beta.mu.t)
Xb.logsig.t    <- getXBeta(X = X.logsig.t, beta = beta.logsig.t)

mu.t <- logsig.t <- matrix(0, ns, nt)
for (t in 1:nt) {
  mu.t[, t] <- Xb.mu.t[, t] + t(chol(Sigma.t)) %*% rnorm(ns) / sqrt(tau.t[t])
  logsig.t[, t] <- Xb.logsig.t[, t] + t(chol(Sigma.t)) %*% rnorm(ns) / sqrt(tau.t[t])
}


xi.t <- 0.01
y.t <- rgev(n = ns * nt, loc = mu.t, exp(logsig.t), xi.t)

Sigma <- solve(Qb.t * tau.t[t])

# initialize values
logsig <- matrix(logsig.t + rnorm(ns * nt, 0, 0.1), ns, nt)
SS <- getGPSS(Qb = Qb.t, param = logsig, Xb = Xb.logsig.t)
curll <- matrix(0, ns, nt)
for (t in 1:nt) {
  curll[, t] <- dgev(x = y.t[, t], loc = mu.t[, t], exp(logsig[, t]), xi.t,
                     log = TRUE)
}

niters <- 10000
burn   <- 8000
logsig.keep <- array(0, dim = c(niters, ns, nt))
acc.logsig <- att.logsig <- MH.logsig <- matrix(0.1, ns, nt)

set.seed(3366)  # demo
for (iter in 1:niters) {
  this.update <- updateLSTest(mu = mu.t, Qb = Qb.t, tau = tau.t, Xb = Xb.logsig.t,
                              y = y.t, logsig = logsig, xi = xi.t,
                              SS = SS, curll = curll, acc = acc.logsig,
                              att = att.logsig, MH = MH.logsig)
  logsig <- this.update$logsig
  SS <- this.update$SS
  curll <- this.update$curll
  acc.logsig <- this.update$acc
  att.logsig <- this.update$att
  logsig.keep[iter, , ] <- logsig

  if (iter < burn / 2) {
    this.update <- mhUpdate(acc = acc.logsig, att = att.logsig, MH = MH.logsig,
                            target.min = 0.4, target.max = 0.7,
                            nattempts = 400)
    acc.logsig <- this.update$acc
    att.logsig <- this.update$att
    MH.logsig  <- this.update$MH
  }

  if (iter %% 500 == 0) {
    par(mfrow = c(3, 3))
    start <- max(1, iter - 20000)
    for (i in 1:3) {
      for (j in 1:3) {
        plot(logsig.keep[start:iter, i, j], type = "l",
             main = paste("logsig: ", round(logsig.t[i, j], 3)),
             ylab = round(acc.logsig[i, j] / att.logsig[i, j], 3),
             xlab = MH.logsig[i, j])
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
X.mu.t <- rX(ns, nt, np)
X.logsig.t <- rX(ns, nt, np)
beta.mu.t <- rnorm(np, 0, 10)
beta.logsig.t <- rnorm(np, 0, 5)

phi.t <- 0.2
s <- cbind(runif(ns), runif(ns))
d <- rdist(s)
Sigma.t <- exp(-d / phi.t)
tau.t   <- rgamma(nt, 1, 1)
Qb.t    <- chol2inv(chol(Sigma.t))
Xb.mu.t     <- getXBeta(X = X.mu.t, beta = beta.mu.t)
Xb.logsig.t <- getXBeta(X = X.logsig.t, beta = beta.logsig.t)

mu.t <- logsig.t <- matrix(0, ns, nt)
for (t in 1:nt) {
  mu.t[, t] <- Xb.mu.t[, t] + t(chol(Sigma.t)) %*% rnorm(ns) / sqrt(tau.t[t])
  logsig.t[, t] <- Xb.logsig.t[, t] + t(chol(Sigma.t)) %*% rnorm(ns) / sqrt(tau.t[t])
}


xi.t <- 0.01
y.t <- rgev(n = ns * nt, loc = mu.t, exp(logsig.t), xi.t)

Sigma <- solve(Qb.t * tau.t[t])

# initialize values
logsig <- matrix(logsig.t + rnorm(ns * nt, 0, 0.1), ns, nt)
SS <- getGPSS(Qb = Qb.t, param = logsig, Xb = Xb.logsig.t)
curll <- matrix(0, ns, nt)
for (t in 1:nt) {
  curll[, t] <- dgev(x = y.t[, t], loc = mu.t[, t], exp(logsig[, t]), xi.t,
                     log = TRUE)
}

niters <- 30000
burn   <- 25000
beta.sd <- 100
beta <- rep(0, np)
beta.keep <- matrix(0, niters, np)
beta.sd.keep <- rep(0, niters)
logsig.keep <- array(0, dim = c(niters, ns, nt))
tau <- rep(1, nt)
tau.keep <- matrix(0, niters, nt)
acc.logsig <- att.logsig <- MH.logsig <- matrix(0.1, ns, nt)

set.seed(3366)  # demo
for (iter in 1:niters) {
  this.update <- updateGPBeta(beta.sd = beta.sd, Qb = Qb.t,
                              param = logsig, X = X.logsig.t, SS = SS, tau = tau)
  beta      <- this.update$beta
  Xb.logsig <- this.update$Xb
  SS        <- this.update$SS
  beta.keep[iter, ] <- beta

  this.update <- updateGPBetaSD(beta = beta, tau.a = 0.1, tau.b = 1)
  beta.sd <- this.update$beta.sd

  beta.sd.keep[iter] <- beta.sd

  this.update <- updateLSTest(mu = mu.t, Qb = Qb.t, tau = tau.t,
                              Xb = Xb.logsig, y = y.t, logsig = logsig,
                              xi = xi.t, SS = SS, curll = curll,
                              acc = acc.logsig, att = att.logsig,
                              MH = MH.logsig)
  logsig <- this.update$logsig
  SS <- this.update$SS
  curll <- this.update$curll
  acc.logsig <- this.update$acc
  att.logsig <- this.update$att
  logsig.keep[iter, , ] <- logsig

  this.update <- updateGPTau(SS = SS, tau.a = 0.1, tau.b = 0.1,
                             ns = ns)
  tau <- this.update$tau
  tau.keep[iter, ] <- tau

  if (iter < burn / 2) {
    this.update <- mhUpdate(acc = acc.logsig, att = att.logsig, MH = MH.logsig,
                            target.min = 0.4, target.max = 0.7,
                            nattempts = 400)
    acc.logsig <- this.update$acc
    att.logsig <- this.update$att
    MH.logsig  <- this.update$MH
  }

  if (iter %% 1000 == 0) {
    par(mfrow = c(4, 3))
    start <- max(1, iter - 20000)
    for (i in 1:2) {
      for (j in 1:3) {
        plot(logsig.keep[start:iter, i, j], type = "l",
             main = paste("logsig: ", round(logsig.t[i, j], 3)),
             ylab = round(acc.logsig[i, j] / att.logsig[i, j], 3),
             xlab = MH.logsig[i, j])
      }
    }

    for (i in 1:3) {
      plot(beta.keep[start:iter, i], type = "l",
           main = paste("beta: ", round(beta.logsig.t[i], 3)))
    }

    for (i in 1:3) {
      plot(tau.keep[start:iter, i], type = "l",
           main = paste("tau: ", round(tau.t[i], 3)))
    }
  }
}