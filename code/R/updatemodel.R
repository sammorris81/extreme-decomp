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

updateBeta1.int <- function(beta.int, beta.mn, SS, tau,
                            beta.time, time,
                            y, theta, theta.xi, mu, ls, xi, thresh,
                            alpha, Qb, curll,
                            acc, att, MH) {
  att <- att + 1
  ns  <- length(beta.int)
  nt  <- length(time)

  for (i in 1:ns) {
    canbeta.int <- rnorm(1, beta.int[i], MH[i])
    canmu <- canbeta.int + beta.time[i] * time  # length nt
    if (any(xi * (y[i, ] - canmu) / exp(ls[i, ]) < -1, na.rm = TRUE)) {
      R <- -Inf
    } else {
      canll <- loglike(y = y[i, ], theta = theta[i, ],
                       theta.xi = theta.xi[i, ],
                       mu = canmu, ls = ls[i, ], xi = xi,
                       thresh = thresh[i, ],
                       alpha = alpha)

      cond.sd <- sqrt(Qb[i, i])
      cond.mn <- beta.mn - Qb[i, -i] %*% (beta.int[-i] - beta.mn) / Qb[i, i]

      R <- sum(canll - curll[i, ]) +
        dnorm(canbeta.int, cond.mn, cond.sd, log = TRUE) -
        dnorm(beta.int[i], cond.mn, cond.sd, log = TRUE)
      # print(sum(canll - curll[i, ]))
      # print(dnorm(canbeta.int, cond.mn, cond.sd, log = TRUE))
      # print(dnorm(beta.int[i], cond.mn, cond.sd, log = TRUE))
      # if (i == ns) {
      #   stop()
      # }

      if (!is.na(exp(R))) { if (log(runif(1)) < R) {
        beta.int[i] <- canbeta.int
        mu[i, ]     <- canmu
        curll[i, ]  <- canll
        acc[i]      <- acc[i] + 1
      }}
    }
  }

  SS <- quad.form(Qb, beta.int - beta.mn)
  results <- list(beta.int = beta.int, SS = SS,
                  mu = mu, curll = curll, acc = acc, att = att)

  return(results)
}

updateBeta1.time <- function(beta.time, beta.mn, SS, tau,
                             beta.int, time,
                             y, theta, theta.xi, mu, ls, xi, thresh,
                             alpha, Qb, curll,
                             acc, att, MH) {
  att <- att + 1
  ns  <- length(beta.time)
  nt  <- length(time)

  for (i in 1:ns) {
    canbeta.time <- rnorm(1, beta.time[i], MH[i])
    canmu <- beta.int[i] + canbeta.time * time  # length nt
    if (any(xi * (y[i, ] - canmu) / exp(ls[i, ]) < -1, na.rm = TRUE)) {
      R <- -Inf
    } else {
      canll <- loglike(y = y[i, ], theta = theta[i, ],
                       theta.xi = theta.xi[i, ],
                       mu = canmu, ls = ls[i, ], xi = xi,
                       thresh = thresh[i, ],
                       alpha = alpha)

      cond.sd <- sqrt(Qb[i, i])
      cond.mn <- beta.mn - Qb[i, -i] %*% (beta.time[-i] - beta.mn) / Qb[i, i]

      R <- sum(canll - curll[i, ]) +
        dnorm(canbeta.time, cond.mn, cond.sd, log = TRUE) -
        dnorm(beta.time[i], cond.mn, cond.sd, log = TRUE)

      if (!is.na(exp(R))) { if (log(runif(1)) < R) {
        beta.time[i] <- canbeta.time
        mu[i, ]      <- canmu
        curll[i, ]   <- canll
        acc[i]       <- acc[i] + 1
      }}
    }
  }

  SS <- quad.form(Qb, beta.time - beta.mn)
  results <- list(beta.time = beta.time, SS = SS,
                  mu = mu, curll = curll, acc = acc, att = att)
  return(results)
}

updateBeta1 <- function(beta.int, beta.int.mn, SS.int, tau.int,
                        beta.time, beta.time.mn, SS.time, tau.time,
                        mu, time,
                        y, theta, theta.xi, ls, xi, thresh,
                        alpha, Qb, curll,
                        acc, att, MH) {
  att <- att + 1
  ns  <- nrow(beta.int)
  nt  <- ncol(beta.int)

  for (t in 1:nt) {
    lp.beta1 <- logpost.beta1.grad(beta.int = beta.int[, t],
                                   beta.int.mn = beta.int.mn, tau.int = tau.int,
                                   beta.time = beta.time[, t],
                                   beta.time.mn = beta.time.mn, time = time[t],
                                   tau.time = tau.time, Qb = Qb, y = y[, t],
                                   ls = ls[, t], xi = xi, theta = theta[, t],
                                   theta.xi = theta.xi[, t],
                                   thresh = thresh[, t], alpha = alpha)

    canbeta.int.mn <- beta.int[, t] + MH[, t]^2 / 2 * lp.beta1$grad.beta.int
    canbeta.time.mn <- beta.time[, t] + MH[, t]^2 / 2 * lp.beta1$grad.beta.time

    canbeta.int  <- rnorm(ns, canbeta.int.mn, MH[, t])
    canbeta.time <- rnorm(ns, canbeta.time.mn, MH[, t])

    canmu <- canbeta.int + canbeta.time * time[t]

    if (any(xi * (y[, t] - canmu) / exp(ls[, t]) < -1)) {
      R <- -Inf
    } else {
      canll <- loglike(y = y[, t], theta = theta[, t],
                       theta.xi = theta.xi[, t], mu = canmu, ls = ls[, t],
                       xi = xi, thresh = thresh[, t], alpha = alpha)

      canSS.int  <- quad.form(Qb, canbeta.int - beta.int.mn)
      canSS.time <- quad.form(Qb, canbeta.time - beta.time.mn)

      lp.beta1 <- logpost.beta1.grad(beta.int = canbeta.int,
                                     beta.int.mn = beta.int.mn,
                                     tau.int = tau.int,
                                     beta.time = canbeta.time,
                                     beta.time.mn = beta.time.mn, time = time[t],
                                     tau.time = tau.time, Qb = Qb, y = y[, t],
                                     ls = ls[, t], xi = xi, theta = theta[, t],
                                     theta.xi = theta.xi[, t],
                                     thresh = thresh[, t], alpha = alpha)

      curbeta.int.mn <- canbeta.int + MH[, t]^2 / 2 * lp.beta1$grad.beta.int
      curbeta.time.mn <- canbeta.time + MH[, t]^2 / 2 * lp.beta1$grad.beta.time

      R <- canll - curll[, t] -
        0.5 * tau.int * (canSS.int  - SS.int[t]) -
        0.5 * tau.time * (canSS.time - SS.time[t]) +
        dnorm(beta.int[, t], curbeta.int.mn, MH[, t], log = TRUE) -
        dnorm(canbeta.int, canbeta.int.mn, MH[, t], log = TRUE) +
        dnorm(beta.time[, t], curbeta.time.mn, MH[, t], log = TRUE) -
        dnorm(canbeta.time, canbeta.time.mn, MH[, t], log = TRUE)

      if (all(!is.na(exp(R)))) {
        keep               <- log(runif(ns)) < R
        beta.int[keep, t]  <- canbeta.int[keep]
        beta.time[keep, t] <- canbeta.time[keep]
        mu[keep, t]        <- canmu[keep]
        curll[keep, t]     <- canll[keep]
        acc[keep, t]       <- acc[keep, t] + 1
      }
    }
  }

  for (t in 1:nt) {
    SS.int[t]  <- quad.form(Qb, beta.int[, t] - beta.int.mn)
    SS.time[t] <- quad.form(Qb, beta.time[, t] - beta.time.mn)
  }

  results <- list(beta.int = beta.int, SS.int = SS.int,
                  beta.time = beta.time, SS.time = SS.time,
                  mu = mu, curll = curll, acc = acc, att = att)
  return(results)
}

updateBeta2.int <- function(beta.int, beta.mn, SS, tau,
                            beta.time, time,
                            y, theta, theta.xi, mu, ls, xi, thresh,
                            alpha, Qb, curll,
                            acc, att, MH) {
  att <- att + 1
  ns  <- length(beta.int)
  nt  <- length(time)

  for (i in 1:ns) {
    canbeta.int <- rnorm(1, beta.int[i], MH[i])
    canls <- canbeta.int + beta.time[i] * time  # length nt
    if (any(xi * (y[i, ] - mu[i, ]) / exp(canls) < -1, na.rm = TRUE)) {
      R <- -Inf
    } else {
      canll <- loglike(y = y[i, ], theta = theta[i, ],
                       theta.xi = theta.xi[i, ],
                       mu = mu[i, ], ls = canls, xi = xi,
                       thresh = thresh[i, ],
                       alpha = alpha)

      cond.sd <- sqrt(Qb[i, i])
      cond.mn <- beta.mn - Qb[i, -i] %*% (beta.int[-i] - beta.mn) / Qb[i, i]

      R <- sum(canll - curll[i, ]) +
        dnorm(canbeta.int, cond.mn, cond.sd, log = TRUE) -
        dnorm(beta.int[i], cond.mn, cond.sd, log = TRUE)

      if (!is.na(exp(R))) { if (log(runif(1)) < R) {
        beta.int[i] <- canbeta.int
        ls[i, ]     <- canls
        curll[i, ]  <- canll
        acc[i]      <- acc[i] + 1
      }}
    }
  }

  SS <- quad.form(Qb, beta.int - beta.mn)
  results <- list(beta.int = beta.int, SS = SS,
                  ls = ls, curll = curll, acc = acc, att = att)

  return(results)
}

updateBeta2.time <- function(beta.time, beta.mn, SS, tau,
                             beta.int, time,
                             y, theta, theta.xi, mu, ls, xi, thresh,
                             alpha, Qb, curll,
                             acc, att, MH) {
  att <- att + 1
  ns  <- length(beta.time)
  nt  <- length(time)

  for (i in 1:ns) {
    canbeta.time <- rnorm(1, beta.time[i], MH[i])
    canls <- beta.int[i] + canbeta.time * time  # length nt
    if (any(xi * (y[i, ] - mu[i, ]) / exp(canls) < -1, na.rm = TRUE)) {
      R <- -Inf
    } else {
      canll <- loglike(y = y[i, ], theta = theta[i, ],
                       theta.xi = theta.xi[i, ],
                       mu = mu[i, ], ls = canls, xi = xi,
                       thresh = thresh[i, ],
                       alpha = alpha)

      cond.sd <- sqrt(Qb[i, i])
      cond.mn <- beta.mn - Qb[i, -i] %*% (beta.time[-i] - beta.mn) / Qb[i, i]

      R <- sum(canll - curll[i, ]) +
        dnorm(canbeta.time, cond.mn, cond.sd, log = TRUE) -
        dnorm(beta.time[i], cond.mn, cond.sd, log = TRUE)

      if (!is.na(exp(R))) { if (log(runif(1)) < R) {
        beta.time[i] <- canbeta.time
        ls[i, ]      <- canls
        curll[i, ]   <- canll
        acc[i]       <- acc[i] + 1
      }}
    }
  }

  SS <- quad.form(Qb, beta.time - beta.mn)
  results <- list(beta.time = beta.time, SS = SS,
                  ls = ls, curll = curll, acc = acc, att = att)
  return(results)
}

updateBeta2 <- function(beta, tau, beta.mn, SS, mu, time,
                        y, theta, theta.xi, ls, xi, thresh,
                        alpha, Qb, curll, acc, att, MH) {
  # beta(ns x nt x 2): spatially varying intercept and time coefficient
  # tau(2): precision for GPs for beta1 and beta2
  # mn(2): overall means for GPs for beta1 and beta2
  # SS(nt x 2): Sum of squares for GP using correlation matrix (no tau)
  # time(nt): vector of standardized times
  acc <- acc + 1
  ns <- dim(beta)[1]
  nt <- dim(beta)[2]

  for (i in 1:ns) {
    VVV     <- Qb[i, i] * tau
    MMM.adj <- Qb[i, -i] / Qb[i, i]
    for (t in 1:nt) {
      canbeta1 <- rnorm(1, beta[i, t, 1], MH[i, t])
      canbeta2 <- rnorm(1, beta[i, t, 2], MH[i, t])
      canls    <- canbeta1 + canbeta2 * time[t]
      # using conditional normal prior
      MMM1 <- beta.mn[1] - MMM.adj %*% (beta[-i, t, 1] - beta.mn[1])
      MMM2 <- beta.mn[2] - MMM.adj %*% (beta[-i, t, 2] - beta.mn[2])

      canll <- loglike(y = y[i, t], theta = theta[i, t],
                       theta.xi = theta.xi[i, t], mu = mu[i, t], ls = canls,
                       xi = xi, thresh = thresh[i, t], alpha = alpha)

      R <- canll - curll[i, t] -
        dnorm(canbeta1, MMM1, 1 / sqrt(VVV[1]), log = TRUE) -
        dnorm(beta[i, t, 1], MMM1, 1 / sqrt(VVV[1]), log = TRUE) +
        dnorm(canbeta2, MMM2, 1 / sqrt(VVV[2]), log = TRUE) -
        dnorm(beta[i, t, 2], MMM2, 1 / sqrt(VVV[2]), log = TRUE)

      if (!any(is.na(exp(R)))) { if (log(runif(1)) < R) {
        beta[i, t, 1] <- canbeta1
        beta[i, t, 2] <- canbeta2
        ls[i, t]      <- canls
        curll[i, t]   <- canll
        acc[i, t]     <- acc[i, t] + 1
      }}
    }
  }

  for (t in 1:nt) {
    SS[t, 1] <- quad.form(Qb, beta[, t, 1] - beta.mn[1])
    SS[t, 2] <- quad.form(Qb, beta[, t, 2] - beta.mn[2])
  }

  results <- list(beta = beta, SS = SS, curll = curll, acc = acc, att = att)
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

updateGPMean <- function(beta.sd, Qb, beta.int, tau.int, beta.time, tau.time) {

  tXQ  <- colSumsC(Qb)  # at the moment, same across all times
  tXQX <- sum(tXQ)
  VVV.int  <- tau.int * tXQX + 1 / beta.sd^2  # same for all times
  VVV.time <- tau.time * tXQX + 1 / beta.sd^2  # same for all times
  MMM.int <- tau.int * tXQ %*% beta.int
  MMM.time <- tau.time * tXQ %*% beta.time

  VVV.int <- 1 / VVV.int
  VVV.time <- 1 / VVV.time

  mn.int  <- MMM.int * VVV.int
  mn.time <- MMM.time * VVV.time
  sd.int  <- sqrt(VVV.int)
  sd.time <- sqrt(VVV.time)

  beta.int.mn  <- rnorm(1, mn.int, sd.int)
  beta.time.mn <- rnorm(1, mn.time, sd.time)

  SS.int  <- quad.form(Qb, beta.int - beta.int.mn)
  SS.time <- quad.form(Qb, beta.time - beta.time.mn)

  results <- list(beta.int.mn = beta.int.mn, SS.int = SS.int,
                  beta.time.mn = beta.time.mn, SS.time = SS.time)
  return(results)
}

updateGPTau <- function(SS.int, SS.time, tau.a, tau.b, ns) {
  nt <- length(SS.int)

  tau.int  <- rgamma(1, ns / 2 + tau.a, SS.int / 2 + tau.b)
  tau.time <- rgamma(1, ns / 2 + tau.a, SS.time / 2 + tau.b)

  results <- list(tau.int = tau.int, tau.time = tau.time)
  return(results)
}

#### TODO: Update updateBW and test
updateBW <- function(bw, bw.min, bw.max, Qb, logdetQb, d,
                     beta.int, tau.int, SS.int, beta.int.mn,
                     beta.time, tau.time, SS.time, beta.time.mn,
                     acc, att, MH) {
  # update the bandwidth term for the gaussian process
  ns <- nrow(beta.int)

  att <- att + 1
  bw.star <- transform$logit(bw, bw.min, bw.max)
  canbw.star <- rnorm(1, bw.star, MH)
  canbw <- transform$inv.logit(canbw.star, bw.min, bw.max)
  canSigma.chol <- chol(exp(-d / canbw))
  canQb <- chol2inv(canSigma.chol)
  canlogdetQb <- -logdet(chol = canSigma.chol)

  canSS.int <- canSS.time <- rep(0, 2)
  for (i in 1:2) {
    canSS.int[i]  <- quad.form(canQb, beta.int[, i] - beta.int.mn[i])
    canSS.time[i] <- quad.form(canQb, beta.time[, i] - beta.time.mn[i])
  }
  # canSS1 <- getGPSS(Qb = canQb, param = , Xb = Xb1)
  # canSS2 <- getGPSS(Qb = canQb, param = ls, Xb = Xb2)

  R <- 0
  for (i in 1:2) {
    R <- R - 0.5 * tau.int[i] * (canSS.int[i] - SS.int[i])
    R <- R - 0.5 * tau.time[i] * (canSS.time[i] - SS.time[i])
  }
  # For R, multiply by 4 diff(logdet) because of 4 GPs
  R <- R + 0.5 * 4 * (canlogdetQb - logdetQb) +
    log(canbw - bw.min) + log(bw.max - canbw) -  # Jacobian of the prior
    log(bw - bw.min) - log(bw.max - bw)

  # print(R)

  if (!is.na(exp(R))) { if (log(runif(1)) < R) {
    acc <- acc + 1
    SS.int <- canSS.int
    SS.time <- canSS.time
    bw <- canbw
    Qb <- canQb
    logdetQb <- canlogdetQb
  }}

  results <- list(bw = bw, Qb = Qb, logdetQb = logdetQb, SS.int = SS.int,
                  SS.time = SS.time, acc = acc, att = att)
  return(results)
}

# updateGPBW <- function(bw, bw.min, bw.max, Qb, logdetQb, d,
#                        mu, Xb1, tau1, SS1, ls, Xb2, tau2, SS2,
#                        acc, att, MH) {
#   # update the bandwidth term for the gaussian process
#   ns <- nrow(Xb1)
#   nt <- ncol(Xb1)
#   att <- att + 1
#   bw.star <- transform$logit(bw, bw.min, bw.max)
#   canbw.star <- rnorm(1, bw.star, MH)
#   canbw <- transform$inv.logit(canbw.star, bw.min, bw.max)
#
#   canQb <- chol2inv(chol(exp(-d / canbw)))
#   canlogdetQb <- logdet(canQb)
#
#   canSS1 <- getGPSS(Qb = canQb, param = mu, Xb = Xb1)
#   canSS2 <- getGPSS(Qb = canQb, param = ls, Xb = Xb2)
#
#   # For R, multiply nt * 2 for 2 spatially varying terms
#   R <- - 0.5 * sum(tau1 * (canSS1 - SS1)) - 0.5 * sum(tau2 * (canSS2 - SS2)) +
#     0.5 * (nt * 2) * (canlogdetQb - logdetQb) +
#     log(canbw - bw.min) + log(bw.max - canbw) -  # Jacobian of the prior
#     log(bw - bw.min) - log(bw.max - bw)
#
#   # print(R)
#
#   if (!is.na(exp(R))) { if (log(runif(1)) < R) {
#     acc <- acc + 1
#     SS1 <- canSS1
#     SS2 <- canSS2
#     bw <- canbw
#     Qb <- canQb
#     logdetQb <- canlogdetQb
#   }}
#
#   results <- list(bw = bw, Qb = Qb, logdetQb = logdetQb, SS1 = SS1, SS2 = SS2,
#                   acc = acc, att = att)
#   return(results)
# }

updateXi <- function(xi, xi.min, xi.max, xi.mn, xi.sd, y, mu, ls, curll, theta,
                     theta.xi, thresh, alpha, acc, att, MH) {
  # update xi term
  # using a slightly more complicated update with truncated normals
  # prior := N(xi.mn, xi.sd)
  # cand  := TN(xi, MH, lower = xi.min, upper = xi.max)
  att <- att + 1

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
  canxi.star  <- runif(1, curlower.U + 1e-6, curupper.U - 1e-6)
  canxi       <- qnorm(canxi.star, xi, MH)
  cantheta.xi <- theta^canxi

  canlower.U <- pnorm(q = this.min, mean = canxi, sd = MH)
  canupper.U <- pnorm(q = this.max, mean = canxi, sd = MH)

  # adjust for asymmetrical candidate: adjusting for truncation
  canlog.cand <- dnorm(canxi, xi, MH, log = TRUE) -
    log(curupper.U - curlower.U)
  curlog.cand <- dnorm(xi, canxi, MH, log = TRUE) -
    log(canupper.U - canlower.U)


  canll  <- loglike(y = y, mu = mu, ls = ls, xi = canxi, theta = theta,
                    theta.xi = cantheta.xi, thresh = thresh, alpha = alpha)
  # We do not need to account for truncation in prior because the
  # scaling cancels out in R, but we do need to account for asymmetrical
  # candidate distribution
  R <- sum(canll - curll) +
    dnorm(canxi, xi.mn, xi.sd, log = TRUE) -
    dnorm(xi, xi.mn, xi.sd, log = TRUE) +
    curlog.cand - canlog.cand  # adjust for asymmetrical candidate


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

# updateXBasisBW <- function(bw, bw.min, bw.max,
#                            beta1, Xb1, mu, tau1, SS1,
#                            beta2, Xb2, ls, tau2, SS2,
#                            Qb, X, dw2, time.interact, acc, att, MH) {
#   # update bandwidth to get the basis functions for X
#   # TODO: adjust for uniform prior
#   # does not impact likelihood
#   att <- att + 1
#   bw.star <- transform$logit(bw, bw.min, bw.max)
#   canbw.star <- rnorm(1, bw.star, MH)
#   canbw <- transform$inv.logit(canbw.star, bw.min, bw.max)
#
#   canB <- makeW(dw2 = dw2, rho = canbw)
#   # canX1  <- rep.basis.X(X = X1, newB = canB, time.interact = time.interact)
#   # canX2  <- rep.basis.X(X = X2, newB = canB, time.interact = time.interact)
#   canX   <- rep.basis.X(X = X, newB = canB, time.interact = time.interact)
#   canXb1 <- getXBeta(X = canX, beta = beta1)
#   canXb2 <- getXBeta(X = canX, beta = beta2)
#
#   # bw.basis only impacts GP prior NOT the likelihood
#   canSS1 <- getGPSS(Qb = Qb, param = mu, Xb = canXb1)
#   canSS2 <- getGPSS(Qb = Qb, param = ls, Xb = canXb2)
#
#   R <- -0.5 * sum(tau1 * (canSS1 - SS1)) - 0.5 * sum(tau2 * (canSS2 - SS2)) +
#     log(canbw - bw.min) + log(bw.max - canbw) -  # Jacobian of the prior
#     log(bw - bw.min) - log(bw.max - bw)
#
#   if (!is.na(R)) { if (log(runif(1)) < R) {
#     acc <- acc + 1
#     bw  <- canbw
#     Xb1 <- canXb1
#     SS1 <- canSS1
#     Xb2 <- canXb2
#     SS2 <- canSS2
#     X   <- canX
#   }}
#
#   results <- list(bw = bw, X = X,
#                   Xb1 = Xb1, SS1 = SS1,
#                   Xb2 = Xb2, SS2 = SS2,
#                   acc = acc, att = att)
#   return(results)
# }
#
# updateMu <- function(mu, tau, Xb, SS, y, theta, theta.xi, ls, xi, thresh, alpha,
#                      Qb, curll, acc, att, MH) {
#   # update mu(s, t)
#   ns <- nrow(mu)
#   nt <- ncol(mu)
#
#   for (t in 1:nt) {
#     att[, t] <- att[, t] + 1
#     canmu.mn <- mu[, t] + MH[, t]^2 / 2 *
#       logpost.mu.grad(mu = mu[, t], Xb = Xb[, t], tau = tau[t], Qb = Qb,
#                       y = y[, t], ls = ls[, t], xi = xi, theta = theta[, t],
#                       theta.xi = theta.xi[, t], thresh = thresh[, t],
#                       alpha = alpha)
#
#     canmu <- rnorm(ns, canmu.mn, MH[, t])
#     if (xi < 0 & any(y[, t] - canmu > -exp(ls[, t]) / xi, na.rm = TRUE)) {
#       R <- -Inf
#     } else if (xi > 0 & any(y[, t] - canmu < -exp(ls[, t]) / xi,
#                             na.rm = TRUE)) {
#       R <- -Inf
#     } else {
#       canll <- loglike(y = y[, t], theta = theta[, t], theta.xi = theta.xi[, t],
#                        mu = canmu, ls = ls[, t], xi = xi, thresh = thresh[, t],
#                        alpha = alpha)
#       canSS <- getGPSS(Qb = Qb, param = canmu, Xb = Xb[, t])
#
#       curmu.mn <- canmu + MH[, t]^2 / 2 *
#         logpost.mu.grad(mu = canmu, Xb = Xb[, t], tau = tau[t], Qb = Qb,
#                         y = y[, t], ls = ls[, t], xi = xi, theta = theta[, t],
#                         theta.xi = theta.xi[, t], thresh = thresh[, t],
#                         alpha = alpha)
#
#       R <- canll - curll[, t] -
#         0.5 * tau[t] * canSS +
#         0.5 * tau[t] * SS[t] +
#         dnorm(mu[, t], curmu.mn, MH[, t], log = TRUE) -
#         dnorm(canmu, canmu.mn, MH[, t], log = TRUE)
#
#       if (!any(is.na(exp(R)))) {
#         keep <- log(runif(ns)) < R
#         mu[keep, t] <- canmu[keep]
#         curll[keep, t] <- canll[keep]
#         acc[keep, t]   <- acc[keep, t] + 1
#       }
#     }
#   }
#
#   SS <- getGPSS(Qb = Qb, param = mu, Xb = Xb)
#
#   results <- list(mu = mu, SS = SS, curll = curll, acc = acc, att = att)
#   return(results)
# }
#
# updateLS <- function(ls, tau, Xb, SS, y, theta, theta.xi, mu, xi, thresh, alpha,
#                      Qb, curll, acc, att, MH) {
#   # update logsig(s, t)
#   ns <- nrow(ls)
#   nt <- ncol(ls)
#
#   for (t in 1:nt) {
#     att[, t] <- att[, t] + 1
#     canls.mn <- ls[, t] + MH[, t]^2 / 2 *
#       logpost.logsig.grad(ls = ls[, t], Xb = Xb[, t], tau = tau[t], Qb = Qb,
#                           y = y[, t], mu = mu[, t], xi = xi, theta = theta[, t],
#                           theta.xi = theta.xi[, t], thresh = thresh[, t],
#                           alpha = alpha)
#     canls <- rnorm(ns, canls.mn, MH[, t])
#     if (xi < 0 & any(y[, t] - mu[, t] > -exp(canls) / xi, na.rm = TRUE)) {
#       R <- -Inf
#     } else if (xi > 0 & any(y[, t] - mu[, t] < -exp(canls) / xi,
#                             na.rm = TRUE)) {
#       R <- -Inf
#     } else {
#       canll <- loglike(y = y[, t], theta = theta[, t], theta.xi = theta.xi[, t],
#                        mu = mu[, t], ls = canls, xi = xi, thresh = thresh[, t],
#                        alpha = alpha)
#       canSS <- getGPSS(Qb = Qb, param = canls, Xb = Xb[, t])
#
#       curls.mn <- canls + MH[, t]^2 / 2 *
#         logpost.logsig.grad(ls = canls, Xb = Xb[, t], tau = tau[t], Qb = Qb,
#                             y = y[, t], mu = mu[, t], xi = xi,
#                             theta = theta[, t], theta.xi = theta.xi[, t],
#                             thresh = thresh[, t],
#                             alpha = alpha)
#
#       R <- canll - curll[, t] -
#         0.5 * tau[t] * canSS +
#         0.5 * tau[t] * SS[t] +
#         dnorm(ls[, t], curls.mn, MH[, t], log = TRUE) -
#         dnorm(canls, canls.mn, MH[, t], log = TRUE)
#
#       if (!any(is.na(exp(R)))) {
#         keep           <- log(runif(ns)) < R
#         ls[keep, t]    <- canls[keep]
#         curll[keep, t] <- canll[keep]
#         acc[keep, t]   <- acc[keep, t] + 1
#       }
#     }
#   }
#
#   SS <- getGPSS(Qb = Qb, param = ls, Xb = Xb)
#
#   results <- list(ls = ls, SS = SS, curll = curll, acc = acc, att = att)
#   return(results)
# }
#
# updateMuTest <- function(mu, tau, Xb, SS, y, ls, xi, Qb, curll, acc, att,
#                          MH) {
#   # update mu(s, t)
#   ns <- nrow(mu)
#   nt <- ncol(mu)
#
#   for (t in 1:nt) {
#     att[, t] <- att[, t] + 1
#     canmu.mn <- mu[, t] + MH[, t]^2 / 2 *
#       logpost.mu.grad.test(mu = mu[, t], Xb = Xb[, t], tau = tau[t], Qb = Qb,
#                            y = y[, t], ls = ls[, t], xi = xi)
#     canmu <- rnorm(ns, canmu.mn, MH[, t])
#     if (xi < 0 & any(y[, t] - canmu > -exp(ls[, t]) / xi, na.rm = TRUE)) {
#       R <- -Inf
#     } else if (xi > 0 & any(y[, t] - canmu < -exp(ls[, t]) / xi,
#                             na.rm = TRUE)) {
#       R <- -Inf
#     } else {
#       canll <- dgev(x = y[, t], loc = canmu, exp(ls[, t]), xi, log = TRUE)
#       canSS <- getGPSS(Qb = Qb, param = canmu, Xb = Xb[, t])
#
#       curmu.mn <- canmu + MH[, t]^2 / 2 *
#         logpost.mu.grad.test(mu = canmu, Xb = Xb[, t], tau = tau[t], Qb = Qb,
#                              y = y[, t], ls = ls[, t], xi = xi)
#
#       R <- canll - curll[, t] -
#         0.5 * tau[t] * canSS +
#         0.5 * tau[t] * SS[t] +
#         dnorm(mu[, t], curmu.mn, MH[, t], log = TRUE) -
#         dnorm(canmu, canmu.mn, MH[, t], log = TRUE)
#
#       if (!any(is.na(exp(R)))) {
#         keep <- log(runif(ns)) < R
#         mu[keep, t] <- canmu[keep]
#         curll[keep, t] <- canll[keep]
#         acc[keep, t]   <- acc[keep, t] + 1
#       }
#     }
#   }
#
#   SS <- getGPSS(Qb = Qb, param = mu, Xb = Xb)
#
#   results <- list(mu = mu, SS = SS, curll = curll, acc = acc, att = att)
#   return(results)
# }
#
# updateLSTest <- function(ls, tau, Xb, SS, y, mu, xi, Qb, curll, acc, att,
#                          MH) {
#   # update mu(s, t)
#   ns <- nrow(mu)
#   nt <- ncol(mu)
#
#   for (t in 1:nt) {
#     att[, t] <- att[, t] + 1
#     canls.mn <- ls[, t] + MH[, t]^2 / 2 *
#       logpost.logsig.grad.test(mu = mu[, t], Xb = Xb[, t], tau = tau[t],
#                                Qb = Qb, y = y[, t], ls = ls[, t], xi = xi)
#     canls <- rnorm(ns, canls.mn, MH[, t])
#     if (xi < 0 & any(y[, t] - mu[, t] > -exp(canls) / xi, na.rm = TRUE)) {
#       R <- -Inf
#     } else if (xi > 0 & any(y[, t] - mu[, t] < -exp(canls) / xi,
#                             na.rm = TRUE)) {
#       R <- -Inf
#     } else {
#       canll <- dgev(x = y[, t], loc = mu[, t], exp(canls), xi, log = TRUE)
#       canSS <- getGPSS(Qb = Qb, param = canls, Xb = Xb[, t])
#
#       curls.mn <- canls + MH[, t]^2 / 2 *
#         logpost.logsig.grad.test(mu = mu[, t], Xb = Xb[, t], tau = tau[t],
#                                  Qb = Qb, y = y[, t], ls = canls, xi = xi)
#
#       R <- canll - curll[, t] -
#         0.5 * tau[t] * canSS +
#         0.5 * tau[t] * SS[t] +
#         dnorm(ls[, t], curls.mn, MH[, t], log = TRUE) -
#         dnorm(canls, canls.mn, MH[, t], log = TRUE)
#
#       if (!any(is.na(exp(R)))) {
#         keep           <- log(runif(ns)) < R
#         ls[keep, t]    <- canls[keep]
#         curll[keep, t] <- canll[keep]
#         acc[keep, t]   <- acc[keep, t] + 1
#       }
#     }
#   }
#
#   SS <- getGPSS(Qb = Qb, param = ls, Xb = Xb)
#
#   results <- list(ls = ls, SS = SS, curll = curll, acc = acc, att = att)
#   return(results)
# }