updateA <- function(A, cuts, bins, Ba, theta, theta.xi, y, mu, ls, xi, thresh,
                    alpha, curll, MH) {
  # update positive stable random effects
  nt <- ncol(A)
  L  <- nrow(A)

  l1    <- get.level(logs = A, cuts = cuts)
  CANA  <- A * exp(MH[l1] * rnorm(nt * L))
  l2    <- get.level(logs = CANA, cuts = cuts)
  q     <- dPS.Rcpp(CANA, alpha, bins) - dPS.Rcpp(A, alpha, bins) +
    dlognormal(A, CANA, matrix(MH[l2], L, nt)) -
    dlognormal(CANA, A, matrix(MH[l1], L, nt))

  for (l in 1:L) {
    canA        <- A
    canA[l, ]   <- CANA[l, ]
    cantheta    <- (Ba %*% canA)^alpha
    cantheta.xi <- cantheta^xi
    # canll     <- loglike(y = y, mu = mu, ls = ls, xi = xi,
    #                      theta = cantheta, thresh = thresh, alpha = alpha)
    canll     <- loglike(y = y, mu = mu, ls = ls, xi = xi, theta = cantheta,
                         theta.xi = cantheta.xi, thresh = thresh, alpha = alpha)

    R    <- colSums(canll - curll) + q[l, ]
    if (all(!is.na(R))) {
      keep <- log(runif(nt)) < R
      A[l, keep]     <- canA[l, keep]
      # theta[, keep]  <- cantheta[, keep]
      theta.xi[, keep] <- cantheta.xi[, keep]
      theta[, keep]  <- theta.xi[, keep]^(1 / xi)
      curll[, keep]  <- canll[, keep]
    }
  }

  results <- list(A = A, l1 = l1, theta = theta, theta.xi = theta.xi,
                  curll = curll)
  return(results)
}

updateXBasisBW <- function(bw, bw.min, bw.max,
                           beta1, Xb1, mu, tau1, SS1,
                           beta2, Xb2, ls, tau2, SS2,
                           Qb, X, dw2, time.interact, acc, att, MH) {
  # update bandwidth to get the basis functions for X
  # TODO: adjust for uniform prior
  # does not impact likelihood
  att <- att + 1
  bw.star <- transform$logit(bw, bw.min, bw.max)
  canbw.star <- rnorm(1, bw.star, MH)
  canbw <- transform$inv.logit(canbw.star, bw.min, bw.max)

  canB <- makeW(dw2 = dw2, rho = canbw)
  # canX1  <- rep.basis.X(X = X1, newB = canB, time.interact = time.interact)
  # canX2  <- rep.basis.X(X = X2, newB = canB, time.interact = time.interact)
  canX   <- rep.basis.X(X = X, newB = canB, time.interact = time.interact)
  canXb1 <- getXBeta(X = canX, beta = beta1)
  canXb2 <- getXBeta(X = canX, beta = beta2)

  # bw.basis only impacts GP prior NOT the likelihood
  canSS1 <- getGPSS(Qb = Qb, param = mu, Xb = canXb1)
  canSS2 <- getGPSS(Qb = Qb, param = ls, Xb = canXb2)

  R <- -0.5 * sum(tau1 * (canSS1 - SS1)) - 0.5 * sum(tau2 * (canSS2 - SS2)) +
    log(canbw - bw.min) + log(bw.max - canbw) -  # Jacobian of the prior
    log(bw - bw.min) - log(bw.max - bw)

  if (!is.na(R)) { if (log(runif(1)) < R) {
    acc <- acc + 1
    bw  <- canbw
    Xb1 <- canXb1
    SS1 <- canSS1
    Xb2 <- canXb2
    SS2 <- canSS2
    X   <- canX
  }}

  results <- list(bw = bw, X = X,
                  Xb1 = Xb1, SS1 = SS1,
                  Xb2 = Xb2, SS2 = SS2,
                  acc = acc, att = att)
  return(results)
}

updateMu <- function(mu, tau, Xb, SS, y, theta, theta.xi, ls, xi, thresh, alpha,
                     Qb, curll, acc, att, MH) {
  # update mu(s, t)
  ns <- nrow(mu)
  nt <- ncol(mu)

  for (t in 1:nt) {
    att[, t] <- att[, t] + 1
    canmu.mn <- mu[, t] + MH[, t]^2 / 2 *
      logpost.mu.grad(mu = mu[, t], Xb = Xb[, t], tau = tau[t], Qb = Qb,
                      y = y[, t], ls = ls[, t], xi = xi, theta = theta[, t],
                      theta.xi = theta.xi[, t], thresh = thresh[, t],
                      alpha = alpha)

    canmu <- rnorm(ns, canmu.mn, MH[, t])
    if (xi < 0 & any(y[, t] - canmu > -exp(ls[, t]) / xi, na.rm = TRUE)) {
      R <- -Inf
    } else if (xi > 0 & any(y[, t] - canmu < -exp(ls[, t]) / xi,
                            na.rm = TRUE)) {
      R <- -Inf
    } else {
      canll <- loglike(y = y[, t], theta = theta[, t], theta.xi = theta.xi[, t],
                       mu = canmu, ls = ls[, t], xi = xi, thresh = thresh[, t],
                       alpha = alpha)
      canSS <- getGPSS(Qb = Qb, param = canmu, Xb = Xb[, t])

      curmu.mn <- canmu + MH[, t]^2 / 2 *
        logpost.mu.grad(mu = canmu, Xb = Xb[, t], tau = tau[t], Qb = Qb,
                        y = y[, t], ls = ls[, t], xi = xi, theta = theta[, t],
                        theta.xi = theta.xi[, t], thresh = thresh[, t],
                        alpha = alpha)

      R <- canll - curll[, t] -
        0.5 * tau[t] * canSS +
        0.5 * tau[t] * SS[t] +
        dnorm(mu[, t], curmu.mn, MH[, t], log = TRUE) -
        dnorm(canmu, canmu.mn, MH[, t], log = TRUE)

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

updateLS <- function(ls, tau, Xb, SS, y, theta, theta.xi, mu, xi, thresh, alpha,
                     Qb, curll, acc, att, MH) {
  # update logsig(s, t)
  ns <- nrow(ls)
  nt <- ncol(ls)

  for (t in 1:nt) {
    att[, t] <- att[, t] + 1
    canls.mn <- ls[, t] + MH[, t]^2 / 2 *
      logpost.logsig.grad(ls = ls[, t], Xb = Xb[, t], tau = tau[t], Qb = Qb,
                          y = y[, t], mu = mu[, t], xi = xi, theta = theta[, t],
                          theta.xi = theta.xi[, t], thresh = thresh[, t],
                          alpha = alpha)
    canls <- rnorm(ns, canls.mn, MH[, t])
    if (xi < 0 & any(y[, t] - mu[, t] > -exp(canls) / xi, na.rm = TRUE)) {
      R <- -Inf
    } else if (xi > 0 & any(y[, t] - mu[, t] < -exp(canls) / xi,
                            na.rm = TRUE)) {
      R <- -Inf
    } else {
      canll <- loglike(y = y[, t], theta = theta[, t], theta.xi = theta.xi[, t],
                       mu = mu[, t], ls = canls, xi = xi, thresh = thresh[, t],
                       alpha = alpha)
      canSS <- getGPSS(Qb = Qb, param = canls, Xb = Xb[, t])

      curls.mn <- canls + MH[, t]^2 / 2 *
        logpost.logsig.grad(ls = canls, Xb = Xb[, t], tau = tau[t], Qb = Qb,
                            y = y[, t], mu = mu[, t], xi = xi,
                            theta = theta[, t], theta.xi = theta.xi[, t],
                            thresh = thresh[, t],
                            alpha = alpha)

      R <- canll - curll[, t] -
        0.5 * tau[t] * canSS +
        0.5 * tau[t] * SS[t] +
        dnorm(ls[, t], curls.mn, MH[, t], log = TRUE) -
        dnorm(canls, canls.mn, MH[, t], log = TRUE)

      if (!any(is.na(exp(R)))) {
        keep           <- log(runif(ns)) < R
        ls[keep, t]    <- canls[keep]
        curll[keep, t] <- canll[keep]
        acc[keep, t]   <- acc[keep, t] + 1
      }
    }
  }

  SS <- getGPSS(Qb = Qb, param = ls, Xb = Xb)

  results <- list(ls = ls, SS = SS, curll = curll, acc = acc, att = att)
  return(results)
}

updateMuTest <- function(mu, tau, Xb, SS, y, ls, xi, Qb, curll, acc, att,
                         MH) {
  # update mu(s, t)
  ns <- nrow(mu)
  nt <- ncol(mu)

  for (t in 1:nt) {
    att[, t] <- att[, t] + 1
    canmu.mn <- mu[, t] + MH[, t]^2 / 2 *
      logpost.mu.grad.test(mu = mu[, t], Xb = Xb[, t], tau = tau[t], Qb = Qb,
                           y = y[, t], ls = ls[, t], xi = xi)
    canmu <- rnorm(ns, canmu.mn, MH[, t])
    if (xi < 0 & any(y[, t] - canmu > -exp(ls[, t]) / xi, na.rm = TRUE)) {
      R <- -Inf
    } else if (xi > 0 & any(y[, t] - canmu < -exp(ls[, t]) / xi,
                            na.rm = TRUE)) {
      R <- -Inf
    } else {
      canll <- dgev(x = y[, t], loc = canmu, exp(ls[, t]), xi, log = TRUE)
      canSS <- getGPSS(Qb = Qb, param = canmu, Xb = Xb[, t])

      curmu.mn <- canmu + MH[, t]^2 / 2 *
        logpost.mu.grad.test(mu = canmu, Xb = Xb[, t], tau = tau[t], Qb = Qb,
                             y = y[, t], ls = ls[, t], xi = xi)

      R <- canll - curll[, t] -
        0.5 * tau[t] * canSS +
        0.5 * tau[t] * SS[t] +
        dnorm(mu[, t], curmu.mn, MH[, t], log = TRUE) -
        dnorm(canmu, canmu.mn, MH[, t], log = TRUE)

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

updateLSTest <- function(ls, tau, Xb, SS, y, mu, xi, Qb, curll, acc, att,
                         MH) {
  # update mu(s, t)
  ns <- nrow(mu)
  nt <- ncol(mu)

  for (t in 1:nt) {
    att[, t] <- att[, t] + 1
    canls.mn <- ls[, t] + MH[, t]^2 / 2 *
      logpost.logsig.grad.test(mu = mu[, t], Xb = Xb[, t], tau = tau[t],
                               Qb = Qb, y = y[, t], ls = ls[, t], xi = xi)
    canls <- rnorm(ns, canls.mn, MH[, t])
    if (xi < 0 & any(y[, t] - mu[, t] > -exp(canls) / xi, na.rm = TRUE)) {
      R <- -Inf
    } else if (xi > 0 & any(y[, t] - mu[, t] < -exp(canls) / xi,
                            na.rm = TRUE)) {
      R <- -Inf
    } else {
      canll <- dgev(x = y[, t], loc = mu[, t], exp(canls), xi, log = TRUE)
      canSS <- getGPSS(Qb = Qb, param = canls, Xb = Xb[, t])

      curls.mn <- canls + MH[, t]^2 / 2 *
        logpost.logsig.grad.test(mu = mu[, t], Xb = Xb[, t], tau = tau[t],
                                 Qb = Qb, y = y[, t], ls = canls, xi = xi)

      R <- canll - curll[, t] -
        0.5 * tau[t] * canSS +
        0.5 * tau[t] * SS[t] +
        dnorm(ls[, t], curls.mn, MH[, t], log = TRUE) -
        dnorm(canls, canls.mn, MH[, t], log = TRUE)

      if (!any(is.na(exp(R)))) {
        keep           <- log(runif(ns)) < R
        ls[keep, t]    <- canls[keep]
        curll[keep, t] <- canll[keep]
        acc[keep, t]   <- acc[keep, t] + 1
      }
    }
  }

  SS <- getGPSS(Qb = Qb, param = ls, Xb = Xb)

  results <- list(ls = ls, SS = SS, curll = curll, acc = acc, att = att)
  return(results)
}

updateGPBeta <- function(mu, beta1.sd, SS1, tau1,
                         ls, beta2.sd, SS2, tau2, X, Qb) {
  # update the beta parameters for the mean of the GPs
  nt <- dim(X)[2]
  np <- dim(X)[3]

  VVV1 <- diag(1 / (beta1.sd^2), np)
  VVV2 <- diag(1 / (beta1.sd^2), np)
  MMM1 <- MMM2 <- rep(0, np)
  for (t in 1:nt) {
    X.t  <- X[, t, ]
    tXQ  <- crossprod(X.t, Qb)  # t(X.t) %*% Qb
    tXQX <- tXQ %*% X.t
    VVV1 <- VVV1 + tau1[t] * tXQX
    VVV2 <- VVV2 + tau2[t] * tXQX
    MMM1 <- MMM1 + tau1[t] * tXQ %*% mu[, t]
    MMM2 <- MMM2 + tau2[t] * tXQ %*% ls[, t]
  }
  VVV1 <- chol2inv(chol(VVV1))
  VVV2 <- chol2inv(chol(VVV2))
  # beta1 and beta2 are not faster with crossprod
  beta1 <- VVV1 %*% MMM1 + t(chol(VVV1)) %*% rnorm(np)
  beta2 <- VVV2 %*% MMM2 + t(chol(VVV2)) %*% rnorm(np)

  Xb1 <- getXBeta(X = X, beta = beta1)
  SS1 <- getGPSS(Qb = Qb, param = mu, Xb = Xb1)
  Xb2 <- getXBeta(X = X, beta = beta2)
  SS2 <- getGPSS(Qb = Qb, param = ls, Xb = Xb2)

  results <- list(beta1 = beta1, Xb1 = Xb1, SS1 = SS1,
                  beta2 = beta2, Xb2 = Xb2, SS2 = SS2)
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

  tau <- rep(0, nt)
  for (t in 1:nt) {
    tau[t] <- rgamma(1, ns / 2 + tau.a, SS[t] / 2 + tau.b)
  }

  results <- list(tau = tau)
  return(results)
}

updateGPBW <- function(bw, bw.min, bw.max, Qb, logdetQb, d,
                       mu, Xb1, tau1, SS1, ls, Xb2, tau2, SS2,
                       acc, att, MH) {
  # update the bandwidth term for the gaussian process
  ns <- nrow(Xb1)
  nt <- ncol(Xb1)
  att <- att + 1
  bw.star <- transform$logit(bw, bw.min, bw.max)
  canbw.star <- rnorm(1, bw.star, MH)
  canbw <- transform$inv.logit(canbw.star, bw.min, bw.max)

  canQb <- chol2inv(chol(exp(-d / canbw)))
  canlogdetQb <- logdet(canQb)

  canSS1 <- getGPSS(Qb = canQb, param = mu, Xb = Xb1)
  canSS2 <- getGPSS(Qb = canQb, param = ls, Xb = Xb2)

  # For R, multiply nt * 2 for 2 spatially varying terms
  R <- - 0.5 * sum(tau1 * (canSS1 - SS1)) - 0.5 * sum(tau2 * (canSS2 - SS2)) +
    0.5 * (nt * 2) * (canlogdetQb - logdetQb) +
    log(canbw - bw.min) + log(bw.max - canbw) -  # Jacobian of the prior
    log(bw - bw.min) - log(bw.max - bw)

  # print(R)

  if (!is.na(exp(R))) { if (log(runif(1)) < R) {
    acc <- acc + 1
    SS1 <- canSS1
    SS2 <- canSS2
    bw <- canbw
    Qb <- canQb
    logdetQb <- canlogdetQb
  }}

  results <- list(bw = bw, Qb = Qb, logdetQb = logdetQb, SS1 = SS1, SS2 = SS2,
                  acc = acc, att = att)
  return(results)
}

updateXi <- function(xi, xi.min, xi.max, xi.mn, xi.sd, y, mu, ls, curll, theta,
                     theta.xi, thresh, alpha, acc, att, MH) {
  # update xi term
  # using a slightly more complicated update with truncated normals
  # prior := TN(xi.mn, xi.sd, lower = xi.min, upper = xi.max)
  # cand  := TN(xi, MH, lower = xi.min, upper = xi.max)
  att <- att + 1
  # xi.star <- transform$logit(xi, xi.min, xi.max)

  # random walk from truncated normal distribution within bounds
  res <- y - mu
  this.min <- xi.min
  this.max <- xi.max
  if (any(res > 0, na.rm = TRUE)) {
    this.min <- max(xi.min, max(-exp(ls[res > 0]) / res[res > 0], na.rm = TRUE))
  }
  if (any(res < 0, na.rm = TRUE)) {
    this.max <- min(xi.max, min(-exp(ls[res < 0]) / res[res < 0], na.rm = TRUE))
  }

  curlower.U <- pnorm(q = this.min, mean = xi, sd = MH)
  curupper.U <- pnorm(q = this.max, mean = xi, sd = MH)
  if (curlower.U < 1e-6) { curlower.U <- 0 }
  if (curupper.U > 0.999999) { curupper.U < 1 }
  canxi.star <- runif(1, curlower.U + 1e-6, curupper.U - 1e-6)
  canxi   <- qnorm(canxi.star, xi, MH)

  canlog.cand <- dnorm(canxi, xi, MH, log = TRUE) -
    log(curupper.U - curlower.U)  # adjustment for truncation

  canlower.U <- pnorm(q = this.min, mean = canxi, sd = MH)
  canupper.U <- pnorm(q = this.max, mean = canxi, sd = MH)

  curlog.cand <- dnorm(xi, canxi, MH, log = TRUE) -
    log(canupper.U - canlower.U)
  cantheta.xi <- theta^canxi
  # xi.star <- transform$logit(xi, xi.min, xi.max)
  # canxi  <- rnorm(1, xi, MH)
  # canxi.star <- rnorm(1, xi.star, MH)
  # canxi   <- transform$inv.logit(canxi.star, xi.min, xi.max)
  # if (canxi < 0 & any(y - mu > -exp(ls) / canxi, na.rm = TRUE)) {
  #   R <- -Inf
  # } else if (canxi > 0 & any(y - mu < -exp(ls) / canxi, na.rm = TRUE)) {
  #   R <- -Inf
  # } else {
    canll  <- loglike(y = y, mu = mu, ls = ls, xi = canxi, theta = theta,
                      theta.xi = cantheta.xi, thresh = thresh, alpha = alpha)
    # We do not need to account for truncation in prior because the
    # scaling cancels out in R, but we do need to account for asymmetrical
    # candidate distribution
    R      <- sum(canll - curll) +
      dnorm(canxi, xi.mn, xi.sd, log = TRUE) -
      dnorm(xi, xi.mn, xi.sd, log = TRUE) +
      curlog.cand - canlog.cand  # adjust for asymmetrical candidate
      # log(canxi - xi.min) + log(xi.max - canxi) -  # Jacobian of the prior
      # log(xi - xi.min) - log(xi.max - xi)
      # dnorm(canxi, xi.mn, xi.sd, log = TRUE) -
      # dnorm(xi, xi.mn, xi.sd, log = TRUE)
  # }

  if (!is.na(R)) { if (log(runif(1)) < R) {
    acc      <- acc + 1
    xi       <- canxi
    theta.xi <- cantheta.xi
    curll    <- canll
  }}

  results <- list(xi = xi, theta.xi = theta.xi, curll = curll,
                  acc = acc, att = att)
  return(results)
}