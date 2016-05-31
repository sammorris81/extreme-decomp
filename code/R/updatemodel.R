updateA <- function(A, cuts, bins, Ba, y, mu, logsig, xi, thresh, alpha, curll,
                    MH) {
  # update positive stable random effects
  nt <- ncol(A)
  L  <- ncol(A)
  
  l1    <- get.level(A, cuts)
  CANA  <- A * exp(MH[l1] * rnorm(nt * L))
  l2    <- get.level(CANA, cuts)
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
  
  results <- list(A = A, theta = theta, curll = curll)
  return(results)
}

updateXBasisBW <- function(bw, bw.min, bw.max,
                           X.mu, beta1, Xb1, mu, tau1,
                           X.sig, beta2, Xb2, logsig, tau2,
                           Qb, dw2, time.interact, acc, att, MH) {
  # update bandwidth to get the basis functions for X
  # does not impact likelihood
  att <- att + 1
  bw.star <- transform$logit(bw, lower = bw.min, upper = bw.max)
  canbw.star <- rnorm(1, bw.star, MH)
  canbw <- transform$inv.logit(canbw.star, lower = bw.min, upper = bw.max)
  
  canB <- makeW(dw2 = dw2, rho = canbw)
  canX.mu  <- rep.basis.X(X = X.mu, newB = canB, time.interact = time.interact)
  canX.sig <- rep.basis.X(X = X.sig, newB = canB, time.interact = time.interact)
  
  canXb1 <- getXBeta(X = canX.mu, beta = beta1)
  canXb2 <- getXBeta(X = canX.sig, beta = beta2)
  
  # need curll - Both mu and sigma
  canll <- curll <- 0
  for (i in 1:nt) {
    canll <- canll - 0.5 * tau1[t] * quad.form(Qb, mu[, t] - canXb1[, t])
    canll <- canll - 0.5 * tau2[t] * quad.form(Qb, logsig[, t] - canXb2[, t])
    
    curll <- curll - 0.5 * tau1[t] * quad.form(Qb, mu[, t] - Xb1[, t])
    curll <- curll - 0.5 * tau2[t] * quad.form(Qb, logsig[, t] - Xb2[, t])
  }
  # need canll - Both mu and sigma
  R <- canll - curll +
    dnorm(canbw.star, 0, 3, log = TRUE) - dnorm(bw.star, 0, 3, log = TRUE)
  
  if (!is.na(R)) { if (log(runif(1)) < R) {
    acc <- acc + 1
    X.mu <- canX.mu
    Xb1 <- canXb1
    X.sig <- canX.sig
    Xb2 <- canXb2
    bw <- canbw
  }}
  
  results <- list(bw = bw, X.mu = X.mu, Xb1 = Xb1, X.sig = X.sig, Xb2 = Xb2,
                  acc = acc, att = att)
  return(results)
}

updateMu <- function(mu, Qb, tau, Xb, logsig, xi, thresh, alpha, curll,
                     acc, att, MH) {
  # update mu(s, t)
  ns <- nrow(mu)
  nt <- ncol(mu)
  
  for (i in 1:ns) { for (t in 1:nt) {
    VVV <- tau[t] * Qb[i, i]
    MMM <- tau[t] * Qb[i, i] * Xb[i, t] -
      tau[t] * sum(Qb[-i, i] * (mu[-i, t] - Xb[-i, t]))
    
    att[i, t] <- att[i, t] + 1
    # canmu <- mu[, t]  # only for the tth day
    canmu <- rnorm(1, mu[i, t], MH[i, t])
    canll <- loglike(y = y[i, t], theta = theta[i, t], mu = canmu,
                     logsig = logsig[i, t], xi = xi, thresh = thresh[i, t],
                     alpha = alpha)
    
    R <- sum(canll - curll[i, t]) +
      dnorm(canmu, MMM / VVV, 1 / sqrt(VVV), log = TRUE) -
      dnorm(mu[i, t], MMM / VVV, 1 / sqrt(VVV), log = TRUE)
    
    if (!is.na(exp(R))) { if (log(runif(1)) < R) {
      mu[i, t] <- canmu
      curll[i, t] <- canll
      acc[i, t] <- acc[i, t] + 1
    }}
  }}
  
  results <- list(mu = mu, curll = curll, acc = acc, att = att)
  return(results)
}

updateLS <- function(logsig, Qb, tau, Xb, mu, xi, thresh, alpha, curll,
                     acc, att, MY) {
  # update logsig(s, t)
  ns <- nrow(logsig)
  nt <- ncol(logsig)
  
  # update mu(s, t)
  for (i in 1:ns) { for (t in 1:nt) {
    VVV <- tau[t] * Qb[i, i]
    MMM <- tau[t] * Qb[i, i] * Xb[i, t] -
      tau[t] * sum(Qb[-i, i] * (logsig[-i, t] - Xb[-i, t]))
    
    att[i, t] <- att[i, t] + 1
    # canmu <- mu[, t]  # only for the tth day
    canlog <- rnorm(1, logsig[i, t], MH[i, t])
    canll  <- loglike(y = y[i, t], theta = theta[i, t], mu = mu[i, t],
                      logsig = canlog, xi = xi, thresh = thresh[i, t],
                      alpha = alpha)
    
    R <- sum(canll - curll[i, t]) +
      dnorm(canlog, MMM / VVV, 1 / sqrt(VVV), log = TRUE) -
      dnorm(logsig[i, t], MMM / VVV, 1 / sqrt(VVV), log = TRUE)
    
    if (!is.na(exp(R))) { if (log(runif(1)) < R) {
      logsig[i, t] <- canlog
      curll[i, t] <- canll
      acc[i, t] <- acc[i, t] + 1
    }}
  }}
  
  results <- list(logsig = logsig, curll = curll, acc = acc, att = att)
  return(results)
}

updateGPBeta <- function(beta, beta.sd, Qb, param, X, tau) {
  # update the beta parameters for the mean of the GPs
  nt <- dim(X)[2]
  np <- length(beta)
  
  VVV <- diag(1 / (beta.sd^2), np)
  MMM <- rep(0, p)
  for (t in 1:nt) {
    X.t  <- X[, t, ]
    tXQ  <- t(X.t) %*% Qb
    tXQX <- tXQ %*% X.t
    VVV  <- VVV + tau[t] * tXQX
    MMM  <- tau[t] * tXQ %*% mu[, t]
  }
  
  VVV <- chol2inv(chol(VVV))
  beta <- VVV %*% MMM + t(chol(VVV)) %*% rnorm(np)
  
  results <- list(beta = beta)
  return(results)
}

updateGPBetaSD <- function(beta, beta.tau.a, beta.tau.b) {
  # update the prior standard deviation on beta
  np <- length(beta)
  beta.var <- 1 / rgamma(1, beta.tau.a + np / 2,
                         beta.tau.b + sum(beta)^2 / 2)
  
  results <- list(beta.sd = sqrt(beta.var))
  return(results)
}

updateGPTau <- function(tau, Qb, param, Xb, tau.a, tau.b) {
  # update the variance parameters for the GPs
  ns <- nrow(param)
  nt <- ncol(param)
  
  for (t in 1:nt) {
    SS <- quad.form(Qb, param[, t] - Xb[, t])
    tau[t] <- rgamma(1, ns / 2 + tau.a, SS / 2 + tau.b)
  }
  
  results <- list(tau = tau)
  return(results)
}

updateGPBW <- function(bw, bw.min, bw.max, Qb, d,
                       mu, Xb1, tau1, logsig, Xb2, tau2,
                       acc, att, MH) {
  # update the bandwidth term for the gaussian process
  att <- att + 1
  bw.star <- log(bw)
  canbw.star <- rnorm(1, bw.star, MH)
  canbw <- exp(canbw.star)
  canQb <- chol2inv(chol(exp(-d / canbw)))
  
  # For R, multiply nt * 2 for 2 spatially varying terms
  R <- 0.5 * (nt * 2) * (logdet(canQb) - logdet(Qb))
  for (t in 1:nt) {
    R <- R - 0.5 * tau1[t] * quad.form(canQb, mu[, t] - Xb1[, t]) -
      0.5 * tau2[t] * quad.form(canQb, logsig[, t] - Xb2[, t]) +
      0.5 * tau1[t] * quad.form(Qb, mu[, t] - Xb1[, t]) +
      0.5 * tau2[t] * quad.form(Qb, logsig[, t] - Xb2[, t])
  }
  
  if (!is.na(exp(R))) { if (log(runif(1)) < R) {
    acc <- acc + 1
    bw <- canbw
    Qb <- canQb
  }}
  
  results <- list(bw = bw, Qb = Qb, acc = acc, att = att)
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