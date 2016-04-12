#############################################################:
###    THE MAIN MCMC FUNCTION TO FIT THE HEATWAVE MODEL   ###:
#############################################################:

#############################################################:
#
# INPUTS:
#
#   y       := ns x nt matrix of data
#   X       := ns x nt x p matrix of covariance for the GPD prob and scale
#   thresh  := ns x nt matrix of GPD thresholds
#   B       := ns x L matrix of PCs
#   alpha   := PS parameter alpha
#   A       := initial A terms
#   beta1   := initial values for the GEV location coefficients
#   beta2   := initial values for the GEV scale coefficients
#   xi      := initial value for the GEV shape parameter
#   beta.sd := prior sd of the elements of beta1 and beta2
#
# OUTPUTS:
#
#   beta1    := Posterior samples of beta1
#   beta2    := Posterior samples of beta2
#   beta.var := Posterior samples of prior variances for beta1 and beta2
#   xi       := Posterior samples of xi
#   logA     := Posterior samples of log(A)
#   theta.mn := Posterior mean of theta
#
################################################################################

ReShMCMC<-function(y, X, X.mu = NULL, X.sig = NULL, thresh, B, alpha,
                   beta1 = NULL, beta1.mn = 0, beta1.sd = 10,
                   beta2 = NULL, beta2.mn = 0, beta2.sd = 1,
                   xi = 0.001, A = NULL,
                   can.mu.sd = 0.1, can.sig.sd = 0.1,
                   beta1.block = FALSE, beta2.block = FALSE,
                   beta1.tau.a = 1, beta1.tau.b = 1, beta1.sd.fix = FALSE,
                   beta2.tau.a = 1, beta2.tau.b = 1, beta2.sd.fix = FALSE,
                   beta1.attempts = 50, beta2.attempts = 50, xi.attempts = 50,
                   keep.burn = FALSE, iters = 5000, burn = 1000, update = 10,
                   iterplot = FALSE){
  require(extRemes)
  # BOOKKEEPING

  ns   <- nrow(y)
  nt   <- ncol(y)
  L    <- ncol(B)

  if (is.null(X.mu)) {
    X.mu <- X
  }
  if (is.null(X.sig)) {
    X.sig <- X
  }

  p.mu    <- dim(X.mu)[3]
  p.sig   <- dim(X.sig)[3]

  miss  <- is.na(y)
  y     <- ifelse(y < thresh, thresh, y)

  # initial missing values to be set at threshold
  # looping over time because threshold is a vector of ns length
  if (any(miss)) {
    missing.times <- which(colSums(miss) > 0)
  }

  # dPS approximation:

  npts      <- 50
  Ubeta     <- qbeta(seq(0, 1, length = npts + 1), 0.5, 0.5)
  MidPoints <- (Ubeta[-1] + Ubeta[-(npts + 1)]) / 2
  BinWidth  <- Ubeta[-1] - Ubeta[-(npts + 1)]
  bins      <- list(npts = npts, MidPoints = MidPoints, BinWidth = BinWidth)

  # INITIAL VALUES:
  if (is.null(beta2)) {
    beta2    <- rep(0, p.sig)
    beta2[1] <- log(sqrt(6) * sd(as.vector(y), na.rm = TRUE) / pi)
  }

  if (is.null(beta1)) {
    beta1    <- rep(0, p.mu)
    beta1[1] <- -0.57722 * beta2[1]
  }

  mu <- logsig <- 0

  for(j in 1:p.mu){
    mu <- mu + X.mu[, , j] * beta1[j]
  }
  for(j in 1:p.sig){
    logsig <- logsig + X.sig[, , j] * beta2[j]
  }

  if (is.null(A)) {
    A <- matrix(1, L, nt)
  } else if (length(A) == 1) {
    A <- matrix(A, L, nt)
  } else if (length(A) != (nt * L)) {
    stop("The length of the initial A must be either 1 or L * nt")
  }

  Ba    <- B^(1 / alpha)
  theta <- (Ba %*% A)^alpha

  curll <- loglike(y, theta, mu, logsig, xi, thresh, alpha)

  # STORAGE:
  keep.beta1   <- matrix(0, iters, p.mu)
  keep.beta2   <- matrix(0, iters, p.sig)
  keep.beta.sd <- matrix(0, iters, 2)  # variance terms for beta priors
  keep.xi      <- rep(0, iters)
  keep.A       <- array(0, dim = c(iters, L, nt))
  if (any(miss)) {
    keep.y <- matrix(0, iters, sum(miss))  # only record the missing data
  } else {
    keep.y <- NULL
  }

  theta.mn   <- 0

  # TUNING:

  cuts <- exp(c(-1, 0, 1, 2, 5, 10))
  MH.a  <- rep(1, 100)
  att.a <- acc.a <- 0 * MH.a
  att.beta1 <- acc.beta1 <- MH.beta1 <- rep(can.mu.sd, p.mu)
  att.beta2 <- acc.beta2 <- MH.beta2 <- rep(can.sig.sd, p.sig)
  att.xi    <- acc.xi    <- MH.xi    <- 0.1

  tic <- proc.time()[3]
  for (iter in 1:iters) {
    ####################################################
    ##############      Random effects A    ############
    ####################################################
    oldA  <- A
    l1    <- get.level(A, cuts)
    CANA  <- A * exp(MH.a[l1] * rnorm(nt * L))
    l2    <- get.level(CANA, cuts)
    q     <- dPS(CANA, alpha, bins) -
             dPS(A, alpha, bins) +
             dlognormal(A, CANA, matrix(MH.a[l2], L, nt)) -
             dlognormal(CANA, A, matrix(MH.a[l1], L, nt))

    for (l in 1:L) {
      canA     <- A
      canA[l, ] <- CANA[l, ]
      cantheta <- (Ba %*% canA)^alpha
      canll    <- loglike(y, cantheta, mu, logsig, xi, thresh, alpha)

      R    <- colSums(canll - curll) + q[l, ]
      keep <- log(runif(nt)) < R

      A[l, keep]     <- canA[l, keep]
      theta[, keep]  <- cantheta[, keep]
      curll[, keep]  <- canll[, keep]
    }

    ####################################################
    ##############      GEV parameters      ############
    ####################################################

    # am splitting out beta1 and beta2 to allow for potentially different
    # covariates for mu and log(sigma)
    if (beta1.block) {
      att.beta1 <- att.beta1 + 1
      canb      <- rnorm(p.mu, beta1, MH.beta1)
      canmu     <- 0
      for (j in 1:p.mu) {
        canmu   <- mu + X.mu[, , j] * canb[j]
      }
      canll     <- loglike(y, theta, canmu, logsig, xi, thresh, alpha)
      R         <- sum(canll - curll) +
                   sum(dnorm(canb, 0, beta1.sd, log = TRUE) -
                       dnorm(beta1, 0, beta1.sd, log = TRUE))
      if (log(runif(1)) < R) {
        acc.beta1 <- acc.beta1 + 1
        beta1     <- canb
        mu        <- canmu
        curll     <- canll
      }
    } else {
      for (j in 1:p.mu) {  # beta1
        att.beta1[j] <- att.beta1[j] + 1
        canb         <- rnorm(1, beta1[j], MH.beta1[j])
        canmu        <- mu + X.mu[, , j] * (canb - beta1[j])
        canll        <- loglike(y, theta, canmu, logsig, xi, thresh, alpha)
        R            <- sum(canll - curll) +
          dnorm(canb, beta1.mn, beta1.sd, log = TRUE) -
          dnorm(beta1[j], beta1.mn, beta1.sd, log = TRUE)
        if (log(runif(1)) < R) {
          acc.beta1[j] <- acc.beta1[j] + 1
          beta1[j]        <- canb
          mu              <- canmu
          curll           <- canll
        }
      }
    }

    if (beta2.block) {
      att.beta2 <- att.beta2 + 1
      canb      <- rnorm(p.sig, beta2, MH.beta2)
      canlogs   <- 0
      for (j in 1:p.sig) { # beta2
        canlogs <- canlogs + X.sig[, , j] * canb[j]
      }
      canll      <- loglike(y, theta, mu, canlogs, xi, thresh, alpha)
      R          <- sum(canll - curll) +
        sum(dnorm(canb, 0, beta2.sd, log = TRUE) -
              dnorm(beta2, 0, beta2.sd, log = TRUE))
      if (log(runif(1)) < R) {
        acc.beta2 <- acc.beta2 + 1
        beta2     <- canb
        logsig    <- canlogs
        curll     <- canll
      }
    } else {
      for (j in 1:p.sig) { # beta2
        att.beta2[j] <- att.beta2[j] + 1
        canb       <- rnorm(1, beta2[j], MH.beta2[j])
        canlogs    <- logsig + X.sig[, , j] * (canb - beta2[j])
        canll      <- loglike(y, theta, mu, canlogs, xi, thresh, alpha)
        R          <- sum(canll - curll) +
          dnorm(canb, beta2.mn, beta2.sd, log = TRUE) -
          dnorm(beta2[j], beta2.mn, beta2.sd, log = TRUE)
        if (log(runif(1)) < R) {
          acc.beta2[j] <- acc.beta2[j] + 1
          beta2[j]     <- canb
          logsig       <- canlogs
          curll        <- canll
        }
      }
    }

    # update prior standard deviations: prior IG(a, b)
    if (!beta1.sd.fix) {
      beta1.var <- 1 / rgamma(1, beta1.tau.a + p.mu / 2,
                              beta1.tau.b + sum((beta1 - beta1.mn)^2) / 2)
      beta1.sd  <- sqrt(beta1.var)
    }

    if (!beta2.sd.fix) {
      beta2.var <- 1 / rgamma(1, beta2.tau.a + p.sig / 2,
                              beta2.tau.b + sum((beta2 - beta2.mn)^2) / 2)

      beta2.sd  <- sqrt(beta2.var)
    }

    # xi
    att.xi <- att.xi + 1
    canxi  <- rnorm(1, xi, MH.xi)
    canll  <- loglike(y, theta, mu, logsig, canxi, thresh, alpha)
    R      <- sum(canll - curll) +
              dnorm(canxi, 0, 0.5, log = TRUE) -
              dnorm(xi, 0, 0.5, log = TRUE)

    if (log(runif(1)) < R) {
      acc.xi <- acc.xi + 1
      xi     <- canxi
      curll  <- canll
    }

    # TUNING

    if (iter < burn / 2) {
      for (j in 1:length(MH.a)) {
        acc.a[j] <- acc.a[j] + sum(oldA[l1 == j] != A[l1 == j])
        att.a[j] <- att.a[j] + sum(l1 == j)
        if (att.a[j] > 1000) {
          if (acc.a[j] / att.a[j] < 0.3) { MH.a[j] <- MH.a[j] * 0.8 }
          if (acc.a[j] / att.a[j] > 0.6) { MH.a[j] <- MH.a[j] * 1.2 }
          acc.a[j] <- att.a[j] <- 0
        }
      }

      for (p in 1:p.mu) { if (att.beta1[p] > beta1.attempts) {
        acc.rate <- acc.beta1[p] / att.beta1[p]
        if (acc.rate < 0.3) { MH.beta1[p] <- MH.beta1[p] * 0.8 }
        if (acc.rate > 0.6) { MH.beta1[p] <- MH.beta1[p] * 1.2 }
        acc.beta1[p] <- att.beta1[p] <- 0
      }}

      for (p in 1:p.sig) { if (att.beta2[p] > beta2.attempts) {
        acc.rate <- acc.beta2[p] / att.beta2[p]
        if (acc.rate < 0.3) { MH.beta2[p] <- MH.beta2[p] * 0.8 }
        if (acc.rate > 0.6) { MH.beta2[p] <- MH.beta2[p] * 1.2 }
        acc.beta2[p] <- att.beta2[p] <- 0
      }}

      if (att.xi > xi.attempts) {
        acc.rate <- acc.xi / att.xi
        if (acc.rate < 0.3) { MH.xi <- MH.xi * 0.8 }
        if (acc.rate > 0.6) { MH.xi <- MH.xi * 1.2 }
        acc.xi <- att.xi <- 0
      }
    }

    ####################################################
    ##############      Impute missing      ############
    ####################################################
    if (any(miss)) {
      y.tmp <- y
      for (t in missing.times) {
        # calculate mu and sigma
        miss.t   <- miss[, t]
        sigma    <- exp(logsig[miss.t, t])
        mu_star  <- mu[miss.t, t] + sigma * (theta[miss.t, t]^xi - 1) / xi
        sig_star <- alpha * sigma * theta[miss.t, t]^xi
        xi_star  <- alpha * xi
        y.tmp[miss.t, t] <- revd(n = sum(miss.t), loc = mu_star,
                                 scale = sig_star, shape = xi_star,
                                 type = "GEV")
      }
    }

    #KEEP TRACK OF STUFF:
    keep.beta1[iter, ]    <- beta1
    keep.beta2[iter, ]    <- beta2
    keep.xi[iter]         <- xi
    keep.A[iter, , ]      <- A
    keep.beta.sd[iter, ] <- c(beta1.sd, beta2.sd)
    if (any(miss)) {
      keep.y[iter, ]      <- y.tmp[miss]
    }
    if (iter > burn) {
      theta.mn <- theta.mn + theta / (iters - burn)
    }


    #DISPLAY CURRENT VALUE:

    if (iter %% update == 0) {
      if (iterplot) {
        acc.rate.mu     <- round(acc.beta1 / att.beta1, 3)
        acc.rate.logsig <- round(acc.beta2 / att.beta2, 3)
        acc.rate.xi     <- round(acc.xi / att.xi, 3)

        if (iter > burn) {
          start <- burn + 1
        } else {
          start <- max(iter - 2000, 1)
        }

        if (!exists("this.plot")) {
          this.plot <- "mu"
        }
        # this.plot <- "sig"

        if (this.plot == "mu") {
          if (L <= 9) {
            par(mfrow = c(3, 3))
          } else {
            par(mfrow = c(3, 4))
          }
          for (plot.idx in 1:min(p.mu, 12)) {
            plot(keep.beta1[start:iter, plot.idx],
                 main = bquote(paste(mu, ": ", beta[.(plot.idx)])),
                 xlab = acc.rate.mu[plot.idx], type = "l",
                 ylab = paste("MH =", round(MH.beta1[plot.idx], 3)))
          }
          this.plot <- "sig"
        } else if (this.plot == "sig") {
          if (L <= 9) {
            par(mfrow = c(3, 3))
          } else {
            par(mfrow = c(3, 4))
          }
          for (plot.idx in 1:min(p.sig, 12)) {
            plot(keep.beta2[start:iter, plot.idx],
                 main = bquote(paste(sigma, ": ", beta[.(plot.idx)])),
                 xlab = acc.rate.logsig[plot.idx], type = "l",
                 ylab = paste("MH =", round(MH.beta2[plot.idx], 3)))
          }
          this.plot <- "all"
        } else {
          par(mfrow = c(3, 4))

          plot(keep.beta.sd[start:iter, 1],
               main = bquote(paste(mu, ": ", sigma[beta])), type = "l")

          plot(keep.beta1[start:iter, 1],
               main = bquote(paste(mu, ": ", beta[0])),
               xlab = acc.rate.mu[1], type = "l",
               ylab = paste("MH =", round(MH.beta1[1], 3)))

          plot(keep.beta1[start:iter, 2],
               main = bquote(paste(mu, ": ", beta[t])),
               xlab = acc.rate.mu[2], type = "l",
               ylab = paste("MH =", round(MH.beta1[2], 3)))

          plot(keep.beta1[start:iter, L],
               main = bquote(paste(mu, ": ", beta[L])),
               xlab = acc.rate.mu[L], type = "l",
               ylab = paste("MH =", round(MH.beta1[L], 3)))

          plot(keep.beta.sd[start:iter, 2],
               main = bquote(paste(sigma, ": ", sigma[sigma])), type = "l")

          plot(keep.beta2[start:iter, 1],
               main = bquote(paste(sigma, ": ", beta[0])),
               xlab = acc.rate.logsig[1], type = "l",
               ylab = paste("MH =", round(MH.beta2[1], 3)))

          plot(keep.beta2[start:iter, 2],
               main = bquote(paste(sigma, ": ", beta[t])),
               xlab = acc.rate.logsig[2], type = "l",
               ylab = paste("MH =", round(MH.beta2[2], 3)))

          plot(keep.beta2[start:iter, L],
               main = bquote(paste(sigma, ": ", beta[L])),
               xlab = acc.rate.logsig[L], type = "l",
               ylab = paste("MH =", round(MH.beta2[L], 3)))

          plot(keep.xi[start:iter], main = bquote(xi),
               xlab = acc.rate.xi, type = "l",
               ylab = paste("MH =", round(MH.xi, 3)))

          plot(log(keep.A[start:iter, 1, 1]), main = "log(A[1, 1])", type = "l")

          plot(log(keep.A[start:iter, 2, 1]), main = "log(A[2, 1])", type = "l")

          plot(log(keep.A[start:iter, L, 1]), main = "log(A[L, 1])", type = "l")

          this.plot <- "mu"
        }
      }
      cat("    Finished fit:", iter, "of", iters, "iters \n")
    }

  }#end iter

  toc <- proc.time()[3]

  if (keep.burn) {
    return.iters <- 1:iters
  } else {
    return.iters <- (burn + 1):iters
  }

  if (any(miss)) {
    keep.y <- keep.y[return.iters, , drop = FALSE]
  }

  list(beta1 = keep.beta1[return.iters, , drop = FALSE],
       beta2 = keep.beta2[return.iters, , drop = FALSE],
       xi = keep.xi[return.iters],
       betasd = keep.beta.sd[return.iters, , drop = FALSE],
       theta.mn = theta.mn,
       A = keep.A[return.iters, , , drop = FALSE],
       y.pred = keep.y,
       timing = toc - tic)
}

#############################################################:
###                PREDICT AT NEW LOCATIONS               ###:
#############################################################:

#############################################################:
#
# INPUTS:
#
#   mcmcoutput := list of output from mcmc
#   X          := npred x nt x p matrix of covariance for the GPD prob and scale
#   thresh     := npred x nt matrix of GPD thresholds
#   B          := npred x L matrix of PCs
#   alpha      := PS parameter alpha
#
# OUTPUTS:
#
#   ypred      := npred x nt matrix of predicted y value
#
################################################################################

pred.ReShMCMC <- function (mcmcoutput, X.pred, B, alpha, start = 1, end = NULL,
                           thin = 1, update = NULL) {
  require(extRemes)
  if (is.null(end)) {
    end <- length(mcmcoutput$xi)
  }

  if (length(dim(X.pred)) != 3) {
    stop("X.pred must be an array with dimensions npred x nt x p")
  }

  npred <- dim(X.pred)[1]
  nt    <- dim(X.pred)[2]
  L     <- ncol(B)
  p     <- dim(X.pred)[3]

  # stays the same for all iterations
  Ba    <- B^(1 / alpha)

  # make sure we are iterating over the post burnin samples
  niters <- length(start:end)
  beta1  <- matrix(mcmcoutput$beta1[start:end, , drop = F], niters, p)
  beta2  <- matrix(mcmcoutput$beta2[start:end, , drop = F], niters, p)
  xi     <- mcmcoutput$xi[start:end]
  A      <- mcmcoutput$A[start:end, , , drop = F]

  # storage for predictions
  y.pred <- array(-99999, dim = c(niters, npred, nt))
  iters <- length(start:end)
  for (iter in 1:iters) {
    xi.i <- xi[iter]

    # calculate mu and sigma
    for (t in 1:nt) {
      theta.i <- (Ba %*% A[iter, , t])^alpha
      mu.i  <- X.pred[, t, ] %*% beta1[iter, ]
      sig.i <- exp(X.pred[, t, ] %*% beta2[iter, ])

      mu.star  <- mu.i + sig.i * (theta.i^xi.i - 1) / xi.i
      sig.star <- alpha * sig.i * theta.i^xi.i
      xi.star  <- alpha * xi.i

      y.pred[iter, , t] <- revd(n = npred, loc = mu.star, scale = sig.star,
                                shape = xi.star, type = "GEV")
    }

    if (iter %% update == 0) {
      cat("    Finished pred:", iter, "of", iters, "iters \n")
    }
  }

  return(y.pred)
}


#############################################################:
###           OTHER FUNCTIONS USED IN THE MCMC            ###:
#############################################################:

loglike <- function(y, theta, mu, logsig, xi, thresh, alpha){

  sigma    <- exp(logsig)
  mu_star  <- mu + sigma * ((theta^xi) - 1) / xi
  sig_star <- alpha * sigma * (theta^xi)
  xi_star  <- alpha * xi

  tx       <- (1 + xi_star * (y - mu_star) / sig_star)^(-1 / xi_star)
  ll       <- -tx + (y > thresh) * ((xi_star + 1) * log(tx) - log(sig_star))
  ll       <- ifelse(is.na(y), 0, ll)  # maybe this handles the missing data
  ll       <- ifelse(is.na(ll), -Inf, ll)

  return(ll)
}

grad_loglike_betamu <- function(beta1, X.mu, y, theta, logsig, xi, thresh,
                                alpha) {
  mu <- 0
  p.mu <- dim(X.mu)[3]
  for (j in 1:p.mu) {
    mu <- mu + X.mu[, , j] * beta1[j]
  }

  sigma    <- exp(logsig)
  mu_star  <- mu + sigma * ((theta^xi) - 1) / xi
  sig_star <- alpha * sigma * (theta^xi)
  xi_star  <- alpha * xi

  tx_star  <- (1 + xi_star * (y - mu_star) / sig_star)

  grad <- rep(0, p.mu)
  for (j in 1:p.mu) {
    this.j <- -(tx_star^(-1 / xi_star - 1)) / sig_star +
      (y > thresh) * (xi_star + 1) / (sig_star * tx_star)
    grad[j] <- sum(this.j * X.mu[, , j])
  }
  return(grad)
}

hess_loglike_betamu <- function(beta1, X.mu, y, theta, logsig, xi, thresh,
                                alpha) {
  mu <- 0
  p.mu <- dim(X.mu)[3]
  for (j in 1:p.mu) {
    mu <- mu + X.mu[, , j] * beta1[j]
  }

  sigma    <- exp(logsig)
  mu_star  <- mu + sigma * ((theta^xi) - 1) / xi
  sig_star <- alpha * sigma * (theta^xi)
  xi_star  <- alpha * xi

  tx_star  <- (1 + xi_star * (y - mu_star) / sig_star)

  hess <- matrix(0, p.mu, p.mu)
  for (j in 1:p.mu) {
    this.j <- -(xi_star + 1) * tx_star^(-1 / xi_star - 2) / sig_star^2 +
      (y > thresh) * (xi_star + 1) * xi_star / (sig_star^2 * tx_star^2)
    for (i in j:p.mu) {
      hess[i, j] <- hess[j, i] <- sum(this.j * X.mu[, , i] * X.mu[, , j])
    }
  }

  return(hess)
}


loglike_mu <- function(beta1, X.mu, y, theta, logsig, xi, thresh, alpha) {
  mu <- 0
  p.mu <- dim(X.mu)[3]
  for (j in 1:p.mu) {
    mu <- mu + X.mu[, , j] * beta1[j]
  }

  ll <- loglike(y = y, theta = theta, mu = mu, logsig = logsig, xi = xi,
                thresh = thresh, alpha = alpha)

  return(sum(ll))
}

grad_loglike_betasig <- function(beta2, X.sig, y, theta, mu, xi, thresh,
                                 alpha) {
  logsig <- 0
  p.sig <- dim(X.sig)[3]
  for (j in 1:p.sig) {
    logsig <- logsig + X.sig[, , j] * beta2[j]
  }

  sigma    <- exp(logsig)
  mu_star  <- mu + sigma * ((theta^xi) - 1) / xi
  sig_star <- alpha * sigma * (theta^xi)
  xi_star  <- alpha * xi

  tx_star  <- (1 + xi_star * (y - mu_star) / sig_star)

  grad <- rep(0, p.sig)
  for (j in 1:p.sig) {
    this.j <- (-tx_star^(-1 / xi_star - 1) *
                 xi * (y - mu) / (xi_star * sigma * theta^xi) +
                 (y > thresh) * (-1 + ((xi_star + 1) * xi * (y - mu)) /
                                   (xi_star * sigma * tx_star * theta^xi)))
    grad[j] <- sum(this.j * X.sig[, , j])
  }

  return(grad)
}

hess_loglike_betasig <- function(beta2, X.sig, y, theta, mu, xi, thresh,
                                 alpha) {
  logsig <- 0
  p.sig <- dim(X.sig)[3]
  for (j in 1:p.sig) {
    logsig <- logsig + X.sig[, , j] * beta2[j]
  }

  sigma    <- exp(logsig)
  mu_star  <- mu + sigma * ((theta^xi) - 1) / xi
  sig_star <- alpha * sigma * (theta^xi)
  xi_star  <- alpha * xi

  tx_star  <- (1 + xi_star * (y - mu_star) / sig_star)

  hess <- matrix(0, p.sig, p.sig)
  dtx <- -xi * (y - mu) / (theta^xi * sigma)
  for (i in 1:p.sig) {
    d12 <- (-1 / xi_star - 1) * tx_star^(-1 / xi_star - 2) * dtx^2
    d21 <- tx_star^(-1 / xi_star - 1) * (-dtx)
    this.j <- (d12 + d21) / xi_star
    this.j <- this.j + (y > thresh) * (xi_star + 1) / alpha * (y - mu) / theta^xi *
      (-1 / (sigma * tx_star) - dtx / (tx_star^2 * sigma))
    for (j in i:p.sig) {
      hess[i, j] <- hess[j, i] <- sum(this.j * X.sig[, , i] * X.sig[, , j])
    }
  }

  return(hess)
}

loglike_sig <- function(beta2, X.sig, y, theta, mu, xi, thresh, alpha) {
  logsig <- 0
  p.sig <- dim(X.sig)[3]
  for (j in 1:p.sig) {
    logsig <- logsig + X.sig[, , j] * beta2[j]
  }

  ll <- loglike(y = y, theta = theta, mu = mu, logsig = logsig, xi = xi,
                thresh = thresh, alpha = alpha)

  return(sum(ll))
}


###########  PS functions  ############

ld <- function(u, A, alpha){
  psi <- pi * u
  c <- (sin(alpha * psi) / sin(psi))^(1 / (1 - alpha))
  c <- c * sin((1 - alpha) * psi) / sin(alpha * psi)
  logd <- log(alpha) - log(1 - alpha) - (1 / (1 - alpha)) * log(A) +
          log(c) - c * (1/A^(alpha / (1 - alpha)))
  exp(logd)
}

dPS <- function(A, alpha, bins){
  l <- 0
  for(j in 1:bins$npts){
    l <- l + bins$BinWidth[j] * ld(bins$MidPoints[j], A, alpha)
  }
  l <- ifelse(A > 0, log(l), -Inf)
  return(l)
}

get.level <- function(A, cuts){
  lev <- A * 0 + 1
  for (j in 1:length(cuts)) {
    lev <- ifelse(A > cuts[j], j + 1, lev)
  }
  return(lev)
}

dlognormal <- function(x, mu, sig){
  dnorm(log(x), log(mu), sig, log = T) - log(x)
}


######################################################
#######            DUMB EXAMPLE            ###########
######################################################
if (FALSE) {
  ns       <- 100
  nt       <- 50
  L        <- 2
  p        <- 2
  alpha    <- 0.5
  B        <- matrix(rgamma(ns * L, 2, 3), ns, L)
  B        <- sweep(B, 1, rowSums(B), "/")
  X        <- array(rnorm(ns * nt * p), c(ns, nt, p))
  X[, , 1] <- 1
  y        <- 0.5 * X[, , 2] + 10 * matrix(rnorm(ns * nt), ns, nt)
  thresh   <- matrix(0, ns, nt)

  fit    <- ReShMCMC(y = y, X = X, thresh = thresh, B = B, alpha = alpha)
}