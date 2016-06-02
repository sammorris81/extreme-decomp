updateA <- function(A, cuts, bins, Ba, theta, y, mu, logsig, xi, thresh, alpha,
                    curll, MH) {
  # update positive stable random effects
  nt <- ncol(A)
  L  <- nrow(A)

  l1    <- getLevelCPP(a = A, cuts = cuts)
  CANA  <- A * exp(MH[l1] * rnorm(nt * L))
  l2    <- getLevelCPP(CANA, cuts)
  q     <- dPS(CANA, alpha, bins) - dPS(A, alpha, bins) +
    dlognormal(A, CANA, matrix(MH[l2], L, nt)) -
    dlognormal(CANA, A, matrix(MH[l1], L, nt))

  for (l in 1:L) {
    canA      <- A
    canA[l, ] <- CANA[l, ]
    cantheta  <- (Ba %*% canA)^alpha
    # cantheta.xi <- cantheta^xi  # theta.xi
    canll     <- loglike(y, cantheta, mu, logsig, xi, thresh, alpha)
    # canll     <- loglike(y, cantheta.xi, mu, logsig, xi, thresh, alpha)

    R    <- colSums(canll - curll) + q[l, ]
    if (all(!is.na(R))) {
      keep <- log(runif(nt)) < R

      A[l, keep]     <- canA[l, keep]
      theta[, keep]  <- cantheta[, keep]
      curll[, keep]  <- canll[, keep]
    }
  }

  results <- list(A = A, l1 = l1, theta = theta, curll = curll)
  return(results)
}

updateXBasisBW <- function(bw, bw.min, bw.max, bw.mn, bw.sd,
                           X.mu, beta1, Xb1, mu, tau1, SS1,
                           X.sig, beta2, Xb2, logsig, tau2, SS2,
                           Qb, dw2, time.interact, acc, att, MH) {
  # update bandwidth to get the basis functions for X
  # TODO: adjust for uniform prior
  # does not impact likelihood
  att <- att + 1
  bw.star <- transform$logit(bw, bw.min, bw.max)
  canbw.star <- rnorm(1, bw.star, MH)
  canbw <- transform$inv.logit(canbw.star, bw.min, bw.max)
  # bw.star <- log(bw)
  # canbw.star <- rnorm(1, bw.star, MH)
  # canbw <- exp(canbw.star)

  canB <- makeW(dw2 = dw2, rho = canbw)
  # canB <- getW(rho  = canbw, dw2 = dw2)
  # print(canbw.star)
  canX.mu  <- rep.basis.X(X = X.mu, newB = canB, time.interact = time.interact)
  canX.sig <- rep.basis.X(X = X.sig, newB = canB, time.interact = time.interact)
  canXb1 <- getXBeta(X = canX.mu, beta = beta1)
  canXb2 <- getXBeta(X = canX.sig, beta = beta2)

  # need curll - Both mu and sigma, but not y
  canll <- curll <- 0
  canSS1 <- getGPSS(Qb = Qb, param = mu, Xb = canXb1)
  canSS2 <- getGPSS(Qb = Qb, param = logsig, Xb = canXb2)
  # canSS1 <- rep(0, nt)
  # canSS2 <- rep(0, nt)
  # for (t in 1:nt) {
  #   canSS1[t] <- quad.form(Qb, mu[, t] - canXb1[, t])
  #   canSS2[t] <- quad.form(Qb, logsig[, t] - canXb2[, t])
  #   # canll <- canll - 0.5 * tau1[t] * canSS1[t]
  #   # canll <- canll - 0.5 * tau2[t] * canSS2[t]
  #   #
  #   # curll <- curll - 0.5 * tau1[t] * quad.form(Qb, mu[, t] - Xb1[, t])
  #   # curll <- curll - 0.5 * tau2[t] * quad.form(Qb, logsig[, t] - Xb2[, t])
  # }
  # need canll - Both mu and sigma
  # print(canbw)
  # print(bw)
  # print(canSS1)
  # print(paste("SS1: ", sum(SS1)))
  # print(paste("SS2: ", sum(SS2)))

  # debug
  # X.mu <<- X.mu
  # canX.mu <<- canX.mu
  # canB  <<- canB
  # tau1 <<- tau1
  # beta1 <<- beta1
  # Xb1   <<- Xb1
  # canXb1 <<- canXb1
  # SS1    <<- SS1
  # canSS1 <<- canSS1
  # stop()

  R <- -0.5 * sum(tau1 * (canSS1 - SS1)) - 0.5 * sum(tau2 * (canSS2 - SS2)) +
    log(canbw - bw.min) + log(bw.max - canbw) -  # Jacobian of the prior
    log(bw - bw.min) - log(bw.max - bw)
    # dnorm(canbw, bw.mn, bw.sd, log = TRUE) + log(canbw) -
    # dnorm(bw, bw.mn, bw.sd, log = TRUE) - log(bw)
  # print(sum(canSS1))
  # print(sum(SS1))
  # print(R)
  if (!is.na(R)) { if (log(runif(1)) < R) {
    acc <- acc + 1
    bw <- canbw
    X.mu <- canX.mu
    Xb1 <- canXb1
    SS1 <- canSS1
    X.sig <- canX.sig
    Xb2 <- canXb2
    SS2 <- canSS2
  }}

  results <- list(bw = bw, X.mu = X.mu, Xb1 = Xb1, SS1 = SS1,
                  X.sig = X.sig, Xb2 = Xb2, SS2 = SS2,
                  acc = acc, att = att)
  return(results)
}

updateMuTest <- function(mu, Qb, tau, Xb, y, logsig, xi, SS, curll, acc, att,
                         MH) {
  # update mu(s, t)
  ns <- nrow(mu)
  nt <- ncol(mu)

  for (t in 1:nt) {
    att[, t] <- att[, t] + 1
    # print(ns)
    # print(mu[, t])
    # print(MH[, t])
    canmu.mn <- mu[, t] + MH[, t]^2 / 2 *
      logpost.mu.grad(mu = mu[, t], Xb = Xb[, t], tau = tau[t], Qb = Qb,
                      y = y[, t], logsig = logsig[, t], xi = xi)
    canmu <- rnorm(ns, canmu.mn, MH[, t])
    if (xi < 0 & any(y[, t] - canmu > -exp(logsig[, t]) / xi, na.rm = TRUE)) {
      R <- -Inf
    } else if (xi > 0 & any(y[, t] - canmu < -exp(logsig[, t]) / xi,
                            na.rm = TRUE)) {
      R <- -Inf
    } else {
      canll <- dgev(x = y[, t], loc = canmu, exp(logsig[, t]), xi, log = TRUE)
      canSS <- getGPSS(Qb = Qb, param = canmu, Xb = Xb[, t])

      curmu.mn <- canmu + MH[, t]^2 / 2 *
        logpost.mu.grad(mu = canmu, Xb = Xb[, t], tau = tau[t], Qb = Qb,
                        y = y[, t], logsig = logsig[, t], xi = xi)

      R <- canll - curll[, t] -
        0.5 * tau[t] * canSS +
        0.5 * tau[t] * SS[t] +
        dnorm(mu[, t], curmu.mn, MH[, t], log = TRUE) -
        dnorm(canmu, canmu.mn, MH[, t], log = TRUE)
      # 0.5 * tau[t] * quad.form(Qb, mu[, t] - Xb[, t])

      # if (!is.na(exp(R))) { if (log(runif(1)) < R) {
      #   mu[, t]    <- canmu
      #   SS[t]      <- canSS
      #   curll[, t] <- canll
      #   acc[, t]   <- acc[, t] + 1
      # }}

      if (!any(is.na(exp(R)))) {
        keep <- log(runif(ns)) < R
        mu[keep, t] <- canmu[keep]
        curll[keep, t] <- canll[keep]
        acc[keep, t]   <- acc[keep, t] + 1
      }
    }
  }

  SS <- getGPSS(Qb = Qb, param = mu, Xb = Xb)


  results <- list(mu = mu, SS = SS, curll = curll, acc = acc, att = att)
  return(results)
}

updateMu <- function(mu, Qb, tau, Xb, y, theta, logsig, xi, thresh, alpha,
                     SS, curll, acc, att, MH) {
  # update mu(s, t)
  ns <- nrow(mu)
  nt <- ncol(mu)

  # for (i in 1:ns) { for (t in 1:nt) {
  #   VVV <- tau[t] * Qb[i, i]
  #   MMM <- tau[t] * Qb[i, i] * Xb[i, t] -
  #     tau[t] * sum(Qb[-i, i] * (mu[-i, t] - Xb[-i, t]))
  #
  #   att[i, t] <- att[i, t] + 1
  #   # canmu <- mu[, t]  # only for the tth day
  #   canmu <- rnorm(1, mu[i, t], MH[i, t])
  #   canll <- loglike(y = y[i, t], theta = theta[i, t], mu = canmu,
  #                    logsig = logsig[i, t], xi = xi, thresh = thresh[i, t],
  #                    alpha = alpha)
  #
  #   R <- sum(canll - curll[i, t]) +
  #     dnorm(canmu, MMM / VVV, 1 / sqrt(VVV), log = TRUE) -
  #     dnorm(mu[i, t], MMM / VVV, 1 / sqrt(VVV), log = TRUE)
  #
  #   if (!is.na(exp(R))) { if (log(runif(1)) < R) {
  #     mu[i, t] <- canmu
  #     curll[i, t] <- canll
  #     acc[i, t] <- acc[i, t] + 1
  #   }}
  # }}

  # VVV <- chol2inv(chol(Qb))          # correlation matrix
  # t.chol.VVV <- t(chol(VVV))
  # MMM <- VVV %*% Xb  # mean

  for (t in 1:nt) {
    att[, t] <- att[, t] + 1
    # VVV <- chol2inv(chol(Qb))
    # MMM <- VVV %*% Xb[, t]

    # canmu <- mu[, t]  # only for the tth day
    canmu <- rnorm(ns, mu[, t], MH[, t])
    # canmu <- Xb[, t] + t.chol.VVV %*% rnorm(ns) / sqrt(tau[t])
    canll <- loglike(y = y[, t], theta = theta[, t], mu = canmu,
                     logsig = logsig[, t], xi = xi, thresh = thresh[, t],
                     alpha = alpha)
    canSS <- getGPSS(Qb = Qb, param = canmu, Xb = Xb[, t])
    # canSS <- quad.form(Qb, canmu - Xb[, t])
    # print(canSS)
    # print(sum(curll[, t]))
    # print(sum(canll))

    R <- sum(canll - curll[, t]) -
      0.5 * tau[t] * canSS +
      0.5 * tau[t] * SS[t]
    # print(R)

    if (!is.na(exp(R))) { if (log(runif(1)) < R) {
      mu[, t]    <- canmu
      SS[t]      <- canSS
      curll[, t] <- canll
      acc[, t]   <- acc[, t] + 1
    }}
  }

  results <- list(mu = mu, SS = SS, curll = curll, acc = acc, att = att)
  return(results)
}

updateLSTest <- function(mu, Qb, tau, Xb, y, logsig, xi, SS, curll, acc, att,
                         MH) {
  # update mu(s, t)
  ns <- nrow(mu)
  nt <- ncol(mu)

  for (t in 1:nt) {
    att[, t] <- att[, t] + 1
    # print(ns)
    # print(mu[, t])
    # print(MH[, t])
    canlogs.mn <- logsig[, t] + MH[, t]^2 / 2 *
      logpost.logsig.grad(mu = mu[, t], Xb = Xb[, t], tau = tau[t], Qb = Qb,
                          y = y[, t], logsig = logsig[, t], xi = xi)
    canlogs <- rnorm(ns, canlogs.mn, MH[, t])
    if (xi < 0 & any(y[, t] - mu[, t] > -exp(canlogs) / xi, na.rm = TRUE)) {
      R <- -Inf
    } else if (xi > 0 & any(y[, t] - mu[, t] < -exp(canlogs) / xi,
                            na.rm = TRUE)) {
      R <- -Inf
    } else {
      canll <- dgev(x = y[, t], loc = mu[, t], exp(canlogs), xi, log = TRUE)
      canSS <- getGPSS(Qb = Qb, param = canlogs, Xb = Xb[, t])

      curlogs.mn <- canlogs + MH[, t]^2 / 2 *
        logpost.logsig.grad(mu = mu[, t], Xb = Xb[, t], tau = tau[t], Qb = Qb,
                            y = y[, t], logsig = canlogs, xi = xi)

      R <- canll - curll[, t] -
        0.5 * tau[t] * canSS +
        0.5 * tau[t] * SS[t] +
        dnorm(logsig[, t], curlogs.mn, MH[, t], log = TRUE) -
        dnorm(canlogs, canlogs.mn, MH[, t], log = TRUE)
      # 0.5 * tau[t] * quad.form(Qb, mu[, t] - Xb[, t])

      # if (!is.na(exp(R))) { if (log(runif(1)) < R) {
      #   mu[, t]    <- canmu
      #   SS[t]      <- canSS
      #   curll[, t] <- canll
      #   acc[, t]   <- acc[, t] + 1
      # }}

      if (!any(is.na(exp(R)))) {
        keep <- log(runif(ns)) < R
        logsig[keep, t] <- canlogs[keep]
        curll[keep, t]  <- canll[keep]
        acc[keep, t]    <- acc[keep, t] + 1
      }
    }
  }

  SS <- getGPSS(Qb = Qb, param = logsig, Xb = Xb)


  results <- list(logsig = logsig, SS = SS, curll = curll, acc = acc, att = att)
  return(results)
}

updateLS <- function(logsig, Qb, tau, Xb, y, theta, mu, xi, thresh, alpha,
                     SS, curll, acc, att, MH) {
  # update logsig(s, t)
  ns <- nrow(logsig)
  nt <- ncol(logsig)

  # update mu(s, t)
  # for (i in 1:ns) { for (t in 1:nt) {
  #   VVV <- tau[t] * Qb[i, i]
  #   MMM <- tau[t] * Qb[i, i] * Xb[i, t] -
  #     tau[t] * sum(Qb[-i, i] * (logsig[-i, t] - Xb[-i, t]))
  #
  #   att[i, t] <- att[i, t] + 1
  #   # canmu <- mu[, t]  # only for the tth day
  #   canlog <- rnorm(1, logsig[i, t], MH[i, t])
  #   canll  <- loglike(y = y[i, t], theta = theta[i, t], mu = mu[i, t],
  #                     logsig = canlog, xi = xi, thresh = thresh[i, t],
  #                     alpha = alpha)
  #
  #   R <- sum(canll - curll[i, t]) +
  #     dnorm(canlog, MMM / VVV, 1 / sqrt(VVV), log = TRUE) -
  #     dnorm(logsig[i, t], MMM / VVV, 1 / sqrt(VVV), log = TRUE)
  #
  #   if (!is.na(exp(R))) { if (log(runif(1)) < R) {
  #     logsig[i, t] <- canlog
  #     curll[i, t] <- canll
  #     acc[i, t] <- acc[i, t] + 1
  #   }}
  # }}

  for (t in 1:nt) {
    att[, t] <- att[, t] + 1
    # VVV <- chol2inv(chol(Qb))
    # MMM <- VVV %*% Xb[, t]

    canlogs <- rnorm(ns, logsig[, t], MH[, t])
    canll <- loglike(y = y[, t], theta = theta[, t], mu = mu[, t],
                     logsig = canlogs, xi = xi, thresh = thresh[, t],
                     alpha = alpha)
    canSS <- getGPSS(Qb = Qb, param = canlogs, Xb = Xb[, t])
    # canSS <- quad.form(Qb, canlogs - Xb[, t])

    R <- sum(canll - curll[, t]) -
      0.5 * tau[t] * canSS +
      0.5 * tau[t] * SS[t]

    if (!is.na(exp(R))) { if (log(runif(1)) < R) {
      logsig[, t] <- canlogs
      SS[t]   <- canSS
      curll[, t] <- canll
      acc[, t] <- acc[, t] + 1
    }}
  }

  results <- list(logsig = logsig, SS = SS, curll = curll, acc = acc, att = att)
  return(results)
}

updateGPBeta <- function(beta.sd, Qb, param, X, SS, tau) {
  # update the beta parameters for the mean of the GPs
  nt <- dim(X)[2]
  np <- dim(X)[3]

  VVV <- diag(1 / (beta.sd^2), np)
  MMM <- rep(0, np)
  for (t in 1:nt) {
    X.t  <- X[, t, ]
    tXQ  <- t(X.t) %*% Qb
    tXQX <- tXQ %*% X.t
    VVV  <- VVV + tau[t] * tXQX
    MMM  <- MMM + tau[t] * tXQ %*% param[, t]
  }
  VVV <- chol2inv(chol(VVV))
  beta <- VVV %*% MMM + t(chol(VVV)) %*% rnorm(np)

  Xb <- getXBeta(X = X, beta = beta)
  SS <- getGPSS(Qb = Qb, param = param, Xb = Xb)

  results <- list(beta = beta, Xb = Xb, SS = SS)
  return(results)
}

updateGPBetaSD <- function(beta, tau.a, tau.b) {
  # update the prior standard deviation on beta
  np <- length(beta)
  beta.var <- 1 / rgamma(1, tau.a + np / 2,
                         tau.b + sum(beta)^2 / 2)

  results <- list(beta.sd = sqrt(beta.var))
  return(results)
}

updateGPTau <- function(SS, tau.a, tau.b, ns) {
  # update the variance parameters for the GPs
  nt <- length(SS)

  for (t in 1:nt) {
    tau[t] <- rgamma(1, ns / 2 + tau.a, SS[t] / 2 + tau.b)
  }

  results <- list(tau = tau)
  return(results)
}

updateGPBW <- function(bw, bw.min, bw.max, bw.mn, bw.sd, Qb, d,
                       mu, Xb1, tau1, SS1, logsig, Xb2, tau2, SS2,
                       acc, att, MH) {
  # update the bandwidth term for the gaussian process
  # TODO: adjust for uniform prior
  ns <- nrow(Xb1)
  nt <- ncol(Xb1)
  att <- att + 1
  bw.star <- transform$logit(bw, bw.min, bw.max)
  canbw.star <- rnorm(1, bw.star, MH)
  canbw <- transform$inv.logit(canbw.star, bw.min, bw.max)

  canQb <- chol2inv(chol(exp(-d / canbw)))

  canSS1 <- getGPSS(Qb = canQb, param = mu, Xb = Xb1)
  canSS2 <- getGPSS(Qb = canQb, param = logsig, Xb = Xb2)
  # res.mu <- mu - Xb1
  # res.sig <- logsig - Xb2
  # canSS1 <- canSS2 <- rep(0, nt)
  # R <- 0
  # for (t in 1:nt) {
  #   canSS1[t] <- quad.form(canQb, res.mu[, t])
  #   canSS2[t] <- quad.form(canQb, res.sig[, t])
  #   # R <- R - 0.5 * tau1[t] * quad.form(canQb, res.mu[, t]) -
  #   #   0.5 * tau2[t] * quad.form(canQb, res.sig[, t]) +
  #   #   0.5 * tau1[t] * quad.form(Qb, res.mu[, t]) +
  #   #   0.5 * tau2[t] * quad.form(Qb, res.sig[, t])
  # }

  # For R, multiply nt * 2 for 2 spatially varying terms
  R <- - 0.5 * sum(tau1 * (canSS1 - SS1)) - 0.5 * sum(tau2 * (canSS2 - SS2)) +
    0.5 * (nt * 2) * (logdet(canQb) - logdet(Qb)) +
    log(canbw - bw.min) + log(bw.max - canbw) -  # Jacobian of the prior
    log(bw - bw.min) - log(bw.max - bw)

  # print(R)

  if (!is.na(exp(R))) { if (log(runif(1)) < R) {
    acc <- acc + 1
    SS1 <- canSS1
    SS2 <- canSS2
    bw <- canbw
    Qb <- canQb
  }}

  results <- list(bw = bw, Qb = Qb, SS1 = SS1, SS2 = SS2, acc = acc, att = att)
  return(results)
}

updateXi <- function(xi, y, mu, logsig, curll, theta, thresh, alpha,
                     acc, att, MH) {
  # update xi term
  att <- att + 1
  canxi  <- rnorm(1, xi, MH)
  if (canxi < 0 & any(y - mu > -exp(logsig) / canxi, na.rm = TRUE)) {
    R <- -Inf
  } else if (canxi > 0 & any(y - mu < -exp(logsig) / canxi, na.rm = TRUE)) {
    R <- -Inf
  } else {
    # cantheta.xi <- theta^xi  # theta.xi
    canll  <- loglike(y, theta, mu, logsig, canxi, thresh, alpha)
    # canll  <- loglike(y, cantheta.xi, mu, logsig, canxi, thresh, alpha)
    R      <- sum(canll - curll) +
      dnorm(canxi, 0, 0.5, log = TRUE) -
      dnorm(xi, 0, 0.5, log = TRUE)
  }

  if (!is.na(R)) { if (log(runif(1)) < R) {
    acc   <- acc + 1
    xi    <- canxi
    # theta.xi <- cantheta.xi  # theta.xi
    curll <- canll
  }}

  results <- list(xi = xi, curll = curll, acc = acc, att = att)
  return(results)
}