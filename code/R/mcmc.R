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

ReShMCMC<-function(y, X, X.mu = NULL, X.sig = NULL, s, knots, thresh, B, alpha,
                   beta1 = NULL, beta1.mu = 0, beta1.sd = 10, mu1.sd = 100,
                   beta2 = NULL, beta2.mu = 0, beta2.sd = 1, mu2.sd = 10,
                   xi = 0.001, A = NULL,
                   can.mu.sd = 0.1, can.sig.sd = 0.1,
                   beta1.block = FALSE, beta2.block = FALSE,
                   beta1.tau.a = 0.1, beta1.tau.b = 0.1,
                   beta2.tau.a = 0.1, beta2.tau.b = 0.1,
                   mu.tau.a = 0.1, mu.tau.b = 0.1, mu.sd.fix = FALSE,
                   logsig.tau.a = 1, logsig.tau.b = 1, logsig.sd.fix = FALSE,
                   mu.attempts = 50, logsig.attempts = 50, xi.attempts = 50,
                   bw.basis.init = NULL, bw.basis.random = TRUE,
                   bw.basis.attempts = 50,
                   bw.gp.init = NULL, bw.gp.attempts = 50,
                   time.interact = FALSE,
                   keep.burn = FALSE, iters = 5000, burn = 1000, update = 10,
                   iterplot = FALSE){
  require(extRemes)
  require(fields)
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

  # initialize distance matrix from site to knot
  dw2 <- rdist(s, knots)^2
  dw2[dw2 < 1e-4] <- 0

  # if (bw.basis.random) {
    # initialize basis functions
    if (is.null(bw.basis.init)) {
      bw.basis <- quantile(dw2, 0.1)
    } else {
      bw.basis <- bw.basis.init
    }

    bw.basis.min <- quantile(dw2, 0.025)
    bw.basis.max <- quantile(dw2, 0.975)

    # B.X <- makeW(dw2 = dw2, rho = bw.basis)
    B.X  <- getW(rho = bw.basis, dw2 = dw2)
    X.mu <- add.basis.X(X = X.mu, B = B.X, time.interact = time.interact)
    X.sig <- add.basis.X(X = X.sig, B = B.X, time.interact = time.interact)
  # }

  p.mu    <- dim(X.mu)[3]
  p.sig   <- dim(X.sig)[3]

  miss  <- is.na(y)
  if (length(thresh) == 1) {
    thresh <- matrix(thresh, ns, nt)
  }
  y     <- ifelse(y < thresh, thresh, y)

  # initial missing values to be set at threshold
  # looping over time because threshold is a vector of ns length
  if (any(miss)) {
    missing.times <- which(colSums(miss) > 0)
  }

  # INITIAL VALUES:
  if (is.null(beta2)) {
    beta2    <- rep(0, p.sig)
    beta2[1] <- log(sqrt(6) * sd(as.vector(y), na.rm = TRUE) / pi)
  }

  if (is.null(beta1)) {
    beta1    <- rep(0, p.mu)
    beta1[1] <- -0.57722 * beta2[1]
  }

  # set the initial mu and logsig to the mean
  mu <- Xb1 <- getXBeta(X = X.mu, beta = beta1)
  logsig <- Xb2 <- getXBeta(X = X.sig, beta = beta2)

  # get the initial covariance matrix for the gaussian process
  d <- rdist(s)
  diag(d) <- 0
  bw.gp.min <- quantile(d[upper.tri(d)], 0.001)
  bw.gp.max <- quantile(d[upper.tri(d)], 0.999)
  if (!is.null(bw.gp.init)) {
    bw.gp <- bw.gp.init
  } else {
    bw.gp <- min(d[upper.tri(d)]) * 4
    # bw.gp <- 0.7
  }
  Sigma <- exp(-d / bw.gp)
  Qb    <- chol2inv(chol(Sigma))

  SS1 <- SS2 <- rep(0, nt)

  tau1 <- tau2 <- rep(0, nt)
  for (t in 1:nt) {
    tau1[t] <- 1 / var(mu[, t])
    tau2[t] <- 1 / var(logsig[, t])
  }

  # adjust xi until it is valid for all y, mu, and sig
  if (xi < 0) {
    while (any(y - mu > -exp(logsig) / xi, na.rm = TRUE)) {
      xi <- xi * 0.8
    }
  } else if (xi > 0) {
    while (any(y - mu < -exp(logsig) / xi, na.rm = TRUE)) {
      xi <- xi * 0.8
    }
  }

  # dPS approximation:
  npts      <- 50
  Ubeta     <- qbeta(seq(0, 1, length = npts + 1), 0.5, 0.5)
  MidPoints <- (Ubeta[-1] + Ubeta[-(npts + 1)]) / 2
  BinWidth  <- Ubeta[-1] - Ubeta[-(npts + 1)]
  bins      <- list(npts = npts, MidPoints = MidPoints, BinWidth = BinWidth)

  if (is.null(A)) {
    A <- matrix(1, L, nt)
  } else if (length(A) == 1) {
    A <- matrix(A, L, nt)
  } else if (length(A) != (nt * L)) {
    stop("The length of the initial A must be either 1 or L * nt")
  }
  oldA <- A

  Ba    <- B^(1 / alpha)
  theta <- (Ba %*% A)^alpha  # should be ns x nt
  theta.xi <- theta^xi

  # initial values for loglikelihood and gradients
  curll <- loglike(y, theta, mu, logsig, xi, thresh, alpha)

  # beta1.grad.cur <- grad_logpost_betamu(beta1 = beta1, beta.mu = beta1.mu,
  #                                       beta.sd = beta1.sd, X.mu = X.mu, y = y,
  #                                       theta = theta, logsig = logsig, xi = xi,
  #                                       thresh = thresh, alpha = alpha)
  #
  # beta2.grad.cur <- grad_logpost_betasig(beta2 = beta2, beta.mu = beta2.mu,
  #                                        beta.sd = beta2.sd, X.sig = X.sig,
  #                                        y = y, theta = theta, mu = mu, xi = xi,
  #                                        thresh = thresh, alpha = alpha)
  #
  # beta1.hess.cur <- hess_logpost_betamu(beta1 = beta1, beta.mu = beta1.mu,
  #                                       beta.sd = beta1.sd, X.mu = X.mu, y = y,
  #                                       theta = theta, logsig = logsig,
  #                                       xi = xi, thresh = thresh, alpha = alpha)
  #
  # beta2.hess.cur <- hess_logpost_betasig(beta2 = beta2, beta.mu = beta2.mu,
  #                                        beta.sd = beta2.sd, X.sig = X.sig, y = y,
  #                                        theta = theta, mu = mu, xi = xi,
  #                                        thresh = thresh, alpha = alpha)

  # STORAGE:
  keep.beta1   <- matrix(0, iters, p.mu)
  keep.beta2   <- matrix(0, iters, p.sig)
  keep.tau1    <- matrix(0, iters, nt)
  keep.tau2    <- matrix(0, iters, nt)
  keep.mu      <- array(0, dim = c(iters, 5, 5))
  keep.logsig  <- array(0, dim = c(iters, 5, 5))
  these.sites <- sample(1:ns, 5)
  these.days <- sample(1:nt, 5)
  keep.beta.mu <- matrix(0, iters, 2)  # mean terms for beta priors
  keep.beta.sd <- matrix(0, iters, 2)  # variance terms for beta priors
  keep.xi      <- rep(0, iters)
  keep.bw.basis <- rep(0, iters)  # bandwidth term for basis functions
  keep.bw.gp    <- rep(0, iters)  # bandwidth for GP
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
  att.mu     <- acc.mu     <- MH.mu     <- matrix(0.1, ns, nt)
  att.logsig <- acc.logsig <- MH.logsig <- matrix(0.1, ns, nt)
  att.xi    <- acc.xi    <- MH.xi    <- 0.1
  att.bw.basis <- acc.bw.basis <- MH.bw.basis <- 0.001
  att.bw.gp    <- acc.bw.gp    <- MH.bw.gp    <- 0.05

  tic <- proc.time()[3]
  for (iter in 1:iters) {
    ####################################################
    ##############      Random effects A    ############
    ####################################################
    oldA  <- A
    # print(mean(A))
    this.update <- updateA(A = A, cuts = cuts, bins = bins, Ba = Ba,
                           theta = theta, y = y, mu = mu, logsig = logsig,
                           xi = xi, thresh = thresh, alpha = alpha,
                           curll = curll, MH = MH.a)

    A     <- this.update$A
    l1    <- this.update$l1
    theta <- this.update$theta
    curll <- this.update$curll
    # print(mean(A))
    ####################################################
    ##########      bandwidth for kernels      #########
    ####################################################
    if (bw.basis.random) {
      this.update <- updateXBasisBW(bw = bw.basis, bw.min = bw.basis.min,
                                    bw.max = bw.basis.max, bw.mn = -2,
                                    bw.sd = 1, X.mu = X.mu,
                                    beta1 = beta1, Xb1 = Xb1, mu = mu,
                                    tau1 = tau1, SS1 = SS1, X.sig = X.sig,
                                    beta2 = beta2, Xb2 = Xb2, logsig = logsig,
                                    tau2 = tau2, SS2 = SS2, Qb = Qb, dw2 = dw2,
                                    time.interact = time.interact,
                                    acc = acc.bw.basis, att = att.bw.basis,
                                    MH = MH.bw.basis)

      bw.basis     <- this.update$bw
      X.mu         <- this.update$X.mu
      Xb1          <- this.update$Xb1
      SS1          <- this.update$SS1
      X.sig        <- this.update$X.sig
      Xb2          <- this.update$Xb2
      SS2          <- this.update$SS2
      acc.bw.basis <- this.update$acc
      att.bw.basis <- this.update$att
    }

    ####################################################
    ##############      GEV parameters      ############
    ####################################################
    # am trying a Gaussian process here
    # mu
    this.update <- updateMu(mu = mu, Qb = Qb, tau = tau1, Xb = Xb1,
                            y = y, theta = theta, logsig = logsig, xi = xi,
                            thresh = thresh, alpha = alpha, SS = SS1,
                            curll = curll, acc = acc.mu, att = att.mu,
                            MH = MH.mu)
    mu     <- this.update$mu
    SS1    <- this.update$SS
    curll  <- this.update$curll
    acc.mu <- this.update$acc
    att.mu <- this.update$att

    # logsig
    this.update <- updateLS(logsig = logsig, Qb = Qb, tau = tau2, Xb = Xb2,
                            y = y, theta = theta, mu = mu, xi = xi,
                            thresh = thresh, alpha = alpha, SS = SS2,
                            curll = curll, acc = acc.logsig, att = att.logsig,
                            MH = MH.logsig)
    logsig     <- this.update$logsig
    SS2        <- this.update$SS
    curll      <- this.update$curll
    acc.logsig <- this.update$acc
    att.logsig <- this.update$att

    # xi
    this.update <- updateXi(xi = xi, y = y, mu = mu, logsig = logsig,
                            curll = curll, theta = theta, thresh = thresh,
                            alpha = alpha, acc = acc.xi, att = att.xi,
                            MH = MH.xi)
    xi     <- this.update$xi
    curll  <- this.update$curll
    acc.xi <- this.update$acc
    att.xi <- this.update$att

    ####################################################
    ###########      GEV hyper parameters      #########
    ####################################################
    # None of these terms should impact the log likelihood because they only
    # impact the gaussian process prior terms for mu and logsig

    # beta terms
    this.update <- updateGPBeta(beta = beta1, beta.sd = beta1.sd, Qb = Qb,
                                param = mu, X = X.mu, SS = SS1, tau = tau1)
    beta1 <- this.update$beta
    Xb1   <- this.update$Xb
    SS1   <- this.update$SS

    this.update <- updateGPBeta(beta = beta2, beta.sd = beta2.sd, Qb = Qb,
                                param = logsig, X = X.sig, SS = SS2, tau = tau2)
    beta2 <- this.update$beta
    Xb2   <- this.update$Xb
    SS2   <- this.update$SS

    # beta sd
    this.update <- updateGPBetaSD(beta = beta1, tau.a = beta1.tau.a,
                                  tau.b = beta1.tau.b)
    beta1.sd <- this.update$beta.sd

    this.update <- updateGPBetaSD(beta = beta2, tau.a = beta2.tau.a,
                                  tau.b = beta2.tau.b)
    beta2.sd <- this.update$beta.sd

    # variances
    this.update <- updateGPTau(tau = tau1, SS = SS1,
                               tau.a = mu.tau.a, tau.b = mu.tau.b, ns = ns)
    tau1 <- this.update$tau

    this.update <- updateGPTau(tau = tau2, SS = SS2,
                               tau.a = mu.tau.a, tau.b = mu.tau.b, ns = ns)
    tau2 <- this.update$tau

    # spatial range
    # if (iter> 100) {
    this.update <- updateGPBW(bw = bw.gp, bw.min = bw.gp.min,
                              bw.mn = -2, bw.sd = 1, bw.max = bw.gp.max,
                              Qb = Qb, d = d, mu = mu, Xb1 = Xb1, tau1 = tau1,
                              SS1 = SS1, logsig = logsig, Xb2 = Xb2,
                              tau2 = tau2, SS2 = SS2,
                              acc = acc.bw.gp, att = att.bw.gp, MH = MH.bw.gp)
    bw.gp     <- this.update$bw
    Qb        <- this.update$Qb
    SS1       <- this.update$SS1
    SS2       <- this.update$SS2
    acc.bw.gp <- this.update$acc
    att.bw.gp <- this.update$att
    # }


    # ####################################################
    # ##############      GEV parameters      ############
    # ####################################################
    # # am splitting out beta1 and beta2 to allow for potentially different
    # # covariates for mu and log(sigma)
    # # beta1
    # if (FALSE) {  # random walk proposal block
    #   att.beta1 <- att.beta1 + 1
    #   canb      <- rnorm(p.mu, beta1, MH.beta1)
    #   canmu     <- 0
    #   for (j in 1:p.mu) {
    #     canmu   <- canmu + X.mu[, , j] * canb[j]
    #   }
    #   canll     <- loglike(y, theta, canmu, logsig, xi, thresh, alpha)
    #   # canll     <- loglike(y, theta.xi, canmu, logsig, xi, thresh, alpha)
    #   R         <- sum(canll - curll) +
    #     sum(dnorm(canb, beta1.mu, beta1.sd, log = TRUE) -
    #           dnorm(beta1, beta1.mu, beta1.sd, log = TRUE))
    #   if (log(runif(1)) < R) {
    #     acc.beta1 <- acc.beta1 + 1
    #     beta1     <- canb
    #     mu        <- canmu
    #     curll     <- canll
    #   }
    # } else if (FALSE) {  # sequential random walk proposal
    #   for (j in 1:p.mu) {
    #     att.beta1[j] <- att.beta1[j] + 1
    #     canb         <- rnorm(1, beta1[j], MH.beta1[j])  # simple rw proposal
    #     canmu        <- mu + X.mu[, , j] * (canb - beta1[j])
    #     canll        <- loglike(y, theta, canmu, logsig, xi, thresh, alpha)
    #     # canll        <- loglike(y, theta.xi, canmu, logsig, xi, thresh, alpha)
    #     R            <- sum(canll - curll) +
    #       dnorm(canb, beta1.mu, beta1.sd, log = TRUE) -
    #       dnorm(beta1[j], beta1.mu, beta1.sd, log = TRUE)
    #     if (log(runif(1)) < R) {
    #       acc.beta1[j] <- acc.beta1[j] + 1
    #       beta1[j]        <- canb
    #       mu              <- canmu
    #       curll           <- canll
    #     }
    #   }
    # } else if (iter < 2000) {  # sequential Langevin update
    #   for (j in 1:p.mu) {
    #     att.beta1[j] <- att.beta1[j] + 1
    #     canb     <- beta1
    #     mean.can <- beta1[j] + 0.5 * MH.beta1[j]^2 * beta1.grad.cur[j]
    #     canb[j]  <- rnorm(1, mean.can, MH.beta1[j])
    #     canmu    <- mu + X.mu[, , j] * (canb[j] - beta1[j])
    #     if (xi < 0 & any(y - canmu > -exp(logsig) / xi, na.rm = TRUE)) {
    #       R <- -Inf
    #     } else if (xi > 0 & any(y - canmu < -exp(logsig) / xi, na.rm = TRUE)) {
    #       R <- -Inf
    #     } else {
    #       canll    <- loglike(y, theta, canmu, logsig, xi, thresh, alpha)
    #       # canll    <- loglike(y, theta.xi, canmu, logsig, xi, thresh, alpha)
    #
    #       # langevin update ratio
    #       beta1.grad.can <- grad_logpost_betamu(beta1 = canb, beta.mu = beta1.mu,
    #                                             beta.sd = beta1.sd, X.mu = X.mu,
    #                                             y = y, theta = theta,
    #                                             logsig = logsig, xi = xi,
    #                                             thresh = thresh, alpha = alpha)
    #
    #       mean.cur <- canb[j] + 0.5 * MH.beta1[j]^2 * beta1.grad.can[j]
    #
    #       R <- sum(canll - curll) +
    #         dnorm(canb[j], beta1.mu, beta1.sd, log = TRUE) -
    #         dnorm(beta1[j], beta1.mu, beta1.sd, log = TRUE) +
    #         dnorm(beta1[j], mean.cur, MH.beta1[j], log = TRUE) -
    #         dnorm(canb[j], mean.can, MH.beta1[j], log = TRUE)
    #       if (!is.nan(R)) { if (log(runif(1)) < R) {
    #         acc.beta1[j]   <- acc.beta1[j] + 1
    #         beta1[j]       <- canb[j]
    #         beta1.grad.cur <- beta1.grad.can
    #         mu             <- canmu
    #         curll          <- canll
    #       } }
    #     }
    #   }
    # } else {  # block Langevin update
    #   if (iter == 1) {
    #     print("block gradient update")
    #   }
    #   # MH.beta1 <- rep(mean(MH.beta1), p.mu)
    #   att.beta1 <- att.beta1 + 1
    #   canb      <- beta1
    #   mean.can  <- beta1 + 0.5 * MH.beta1^2 * beta1.grad.cur
    #   canb      <- rnorm(p.mu, mean.can, MH.beta1)
    #   canmu <- 0
    #   for (j in 1:p.mu) {
    #     canmu   <- canmu + X.mu[, , j] * canb[j]
    #   }
    #
    #   if (xi < 0 & any(y - canmu > -exp(logsig) / xi, na.rm = TRUE)) {
    #     R <- -Inf
    #   } else if (xi > 0 & any(y - canmu < -exp(logsig) / xi, na.rm = TRUE)) {
    #     R <- -Inf
    #   } else {
    #     canll     <- loglike(y, theta, canmu, logsig, xi, thresh, alpha)
    #     # canll     <- loglike(y, theta.xi, canmu, logsig, xi, thresh, alpha)
    #
    #     # langevin update ratio
    #     beta1.grad.can <- grad_logpost_betamu(beta1 = canb, beta.mu = beta1.mu,
    #                                           beta.sd = beta1.sd, X.mu = X.mu,
    #                                           y = y, theta = theta,
    #                                           logsig = logsig, xi = xi,
    #                                           thresh = thresh, alpha = alpha)
    #
    #     mean.cur <- canb + 0.5 * MH.beta1^2 * beta1.grad.can
    #
    #     R <- sum(canll - curll) +
    #       sum(dnorm(canb, beta1.mu, beta1.sd, log = TRUE) -
    #             dnorm(beta1, beta1.mu, beta1.sd, log = TRUE)) +
    #       sum(dnorm(beta1, mean.cur, MH.beta1, log = TRUE) -
    #             dnorm(canb, mean.can, MH.beta1, log = TRUE))
    #
    #     if (!is.nan(R)) { if (log(runif(1)) < R) {
    #       acc.beta1   <- acc.beta1 + 1
    #       beta1       <- canb
    #       beta1.grad.cur <- beta1.grad.can
    #       mu             <- canmu
    #       curll          <- canll
    #     } }
    #   }
    # }
    #
    # # beta2
    # if (FALSE) {  # random walk proposal block
    #   att.beta2 <- att.beta2 + 1
    #   canb      <- rnorm(p.sig, beta2, MH.beta2)
    #   canlogs   <- 0
    #   for (j in 1:p.sig) { # beta2
    #     canlogs <- canlogs + X.sig[, , j] * canb[j]
    #   }
    #   canll      <- loglike(y, theta, mu, canlogs, xi, thresh, alpha)
    #   # canll      <- loglike(y, theta.xi, mu, canlogs, xi, thresh, alpha)
    #   R          <- sum(canll - curll) +
    #     sum(dnorm(canb, beta2.mu, beta2.sd, log = TRUE) -
    #           dnorm(beta2, beta2.mu, beta2.sd, log = TRUE))
    #   if (log(runif(1)) < R) {
    #     acc.beta2 <- acc.beta2 + 1
    #     beta2     <- canb
    #     logsig    <- canlogs
    #     curll     <- canll
    #   }
    # } else if (FALSE) {  # sequential random walk proposal
    #   for (j in 1:p.sig) {
    #     att.beta2[j] <- att.beta2[j] + 1
    #     canb       <- rnorm(1, beta2[j], MH.beta2[j])
    #     canlogs    <- logsig + X.sig[, , j] * (canb - beta2[j])
    #     canll      <- loglike(y, theta, mu, canlogs, xi, thresh, alpha)
    #     # canll      <- loglike(y, theta.xi, mu, canlogs, xi, thresh, alpha)
    #     R          <- sum(canll - curll) +
    #       dnorm(canb, beta2.mu, beta2.sd, log = TRUE) -
    #       dnorm(beta2[j], beta2.mu, beta2.sd, log = TRUE)
    #     if (log(runif(1)) < R) {
    #       acc.beta2[j] <- acc.beta2[j] + 1
    #       beta2[j]     <- canb
    #       logsig       <- canlogs
    #       curll        <- canll
    #     }
    #   }
    # } else if (iter < 2000) {  # sequential Langevin update
    #   for (j in 1:p.sig) {
    #     att.beta2[j] <- att.beta2[j] + 1
    #     canb     <- beta2
    #     mean.can <- beta2[j] + 0.5 * MH.beta2[j]^2 * beta2.grad.cur[j]
    #     canb[j]  <- rnorm(1, mean.can, MH.beta2[j])
    #     canlogs  <- logsig + X.sig[, , j] * (canb[j] - beta2[j])
    #
    #     # if we're outside of the parameter space, automatically reject
    #     if (xi < 0 & any(y - mu > -exp(canlogs) / xi, na.rm = TRUE)) {
    #       R <- -Inf
    #     } else if (xi > 0 & any(y - mu < -exp(canlogs) / xi, na.rm = TRUE)) {
    #       R <- -Inf
    #     } else {
    #       canll    <- loglike(y, theta, mu, canlogs, xi, thresh, alpha)
    #       # canll    <- loglike(y, theta.xi, mu, canlogs, xi, thresh, alpha)
    #
    #       # langevin update ratio
    #       beta2.grad.can <- grad_logpost_betasig(beta2 = canb,
    #                                              beta.mu = beta2.mu,
    #                                              beta.sd = beta2.sd,
    #                                              X.sig = X.sig, y = y,
    #                                              theta = theta, mu = mu,
    #                                              xi = xi, thresh = thresh,
    #                                              alpha = alpha)
    #
    #       mean.cur <- canb[j] + 0.5 * MH.beta2[j]^2 * beta2.grad.can[j]
    #
    #       R <- sum(canll - curll) +
    #         dnorm(canb[j], beta2.mu, beta2.sd, log = TRUE) -
    #         dnorm(beta2[j], beta2.mu, beta2.sd, log = TRUE) +
    #         dnorm(beta2[j], mean.cur, MH.beta2[j], log = TRUE) -
    #         dnorm(canb[j], mean.can, MH.beta2[j], log = TRUE)
    #       if (!is.nan(R)) { if (log(runif(1)) < R) {
    #         # print(paste("updating beta", j, "to", canb[j]))
    #         acc.beta2[j]   <- acc.beta2[j] + 1
    #         beta2[j]       <- canb[j]
    #         beta2.grad.cur <- beta2.grad.can
    #         logsig         <- canlogs
    #         curll          <- canll
    #       } }
    #     }
    #
    #     # print(paste("iter", iter, mean(mu)))
    #   }
    # } else {  # block Langevin update
    #   if (iter == 1) {
    #     print("block gradient update")
    #   }
    #   # MH.beta2 <- rep(mean(MH.beta2), p.sig)
    #   att.beta2 <- att.beta2 + 1
    #   mean.can <- beta2 + 0.5 * MH.beta2^2 * beta2.grad.cur
    #   canb     <- rnorm(p.sig, mean.can, MH.beta2)
    #   canlogs <- 0
    #   for (j in 1:p.sig) {
    #     canlogs <- canlogs + X.sig[, , j] * canb[j]
    #   }
    #   if (xi < 0 & any(y - mu > -exp(canlogs) / xi, na.rm = TRUE)) {
    #     R <- -Inf
    #   } else if (xi > 0 & any(y - mu < -exp(canlogs) / xi, na.rm = TRUE)) {
    #     R <- -Inf
    #   } else {
    #     canll <- loglike(y, theta, mu, canlogs, xi, thresh, alpha)
    #     # canll <- loglike(y, theta.xi, mu, canlogs, xi, thresh, alpha)
    #
    #     # langevin update ratio
    #     beta2.grad.can <- grad_logpost_betasig(beta2 = canb, beta.mu = beta2.mu,
    #                                            beta.sd = beta2.sd,
    #                                            X.sig = X.sig, y = y,
    #                                            theta = theta, mu = mu,
    #                                            xi = xi, thresh = thresh,
    #                                            alpha = alpha)
    #
    #     mean.cur <- canb + 0.5 * MH.beta2^2 * beta2.grad.can
    #
    #     R <- sum(canll - curll) +
    #       sum(dnorm(canb, beta2.mu, beta2.sd, log = TRUE) -
    #             dnorm(beta2, beta2.mu, beta2.sd, log = TRUE)) +
    #       sum(dnorm(beta2, mean.cur, MH.beta2, log = TRUE) -
    #             dnorm(canb, mean.can, MH.beta2, log = TRUE))
    #     # print(R)
    #     # print(sum(canll))
    #     # print(sum(curll))
    #     # print(as.numeric(canb))
    #     # print(as.numeric(mean.can))
    #     # print(as.numeric(beta2))
    #     # print(as.numeric(mean.cur))
    #     # stop()
    #     if (!is.nan(R)) { if (log(runif(1)) < R) {
    #       # print(paste("updating beta", j, "to", canb[j]))
    #       acc.beta2      <- acc.beta2 + 1
    #       beta2          <- canb
    #       beta2.grad.cur <- beta2.grad.can
    #       logsig         <- canlogs
    #       curll          <- canll
    #     } }
    #   }
    # }
    #
    # # update beta1.mu: prior N(0, mu1.sd)
    # # mu1.post.var <- 1 / (1 / mu1.sd^2 + p.mu / beta1.sd^2)
    # # mu1.post.mn  <- sum(beta1)/beta1.sd^2 * mu1.post.var
    # # beta1.mu <- rnorm(1, mu1.post.mn, sqrt(mu1.post.var))
    # # beta1.mu <- 0
    #
    # # update beta2.mu: prior N(0, mu2.sd)
    # # mu2.post.var <- 1 / (1 / mu2.sd^2 + p.sig / beta2.sd^2)
    # # mu2.post.mn  <- sum(beta2)/beta2.sd^2 * mu2.post.var
    # # beta2.mu <- rnorm(1, mu2.post.mn, sqrt(mu2.post.var))
    # # beta2.mu <- 0
    #
    # # update prior standard deviations: prior IG(a, b)
    # if (!beta1.sd.fix) {
    #   beta1.var <- 1 / rgamma(1, beta1.tau.a + p.mu / 2,
    #                           beta1.tau.b + sum((beta1 - beta1.mu)^2) / 2)
    #   beta1.sd  <- sqrt(beta1.var)
    # }
    #
    # if (!beta2.sd.fix) {
    #   beta2.var <- 1 / rgamma(1, beta2.tau.a + p.sig / 2,
    #                           beta2.tau.b + sum((beta2 - beta2.mu)^2) / 2)
    #
    #   beta2.sd  <- sqrt(beta2.var)
    # }

    # TUNING

    if (iter < burn / 2) {
      # if (iter > burn * 0.10) {
      for (j in 1:length(MH.a)) {
        acc.a[j] <- acc.a[j] + sum(oldA[l1 == j] != A[l1 == j])
        att.a[j] <- att.a[j] + sum(l1 == j)
        if (att.a[j] > 1000) {
          if (acc.a[j] / att.a[j] < 0.3) { MH.a[j] <- MH.a[j] * 0.8 }
          if (acc.a[j] / att.a[j] > 0.6) { MH.a[j] <- MH.a[j] * 1.2 }
          acc.a[j] <- att.a[j] <- 0
        }
      }
      # }
      this.update <- mhUpdate(acc = acc.bw.basis, att = att.bw.basis,
                              MH = MH.bw.basis, nattempts = bw.basis.attempts,
                              target.min = 0.3, target.max = 0.6,
                              lower = 0.8, higher = 1.2)
      acc.bw.basis <- this.update$acc
      att.bw.basis <- this.update$att
      MH.bw.basis  <- this.update$MH

      this.update <- mhUpdate(acc = acc.mu, att = att.mu, nattempts = 50,
                              MH = MH.mu, target.min = 0.3, target.max = 0.6,
                              lower = 0.8, higher = 1.2)
      acc.mu <- this.update$acc
      att.mu <- this.update$att
      MH.mu  <- this.update$MH

      this.update <- mhUpdate(acc = acc.logsig, att = att.logsig,
                              MH = MH.logsig, nattempts = 50,
                              target.min = 0.3, target.max = 0.6,
                              lower = 0.8, higher = 1.2)
      acc.logsig <- this.update$acc
      att.logsig <- this.update$att
      MH.logsig  <- this.update$MH


      this.update <- mhUpdate(acc = acc.xi, att = att.xi,
                              nattempts = xi.attempts,
                              MH = MH.xi, target.min = 0.3, target.max = 0.6,
                              lower = 0.8, higher = 1.2)
      acc.xi <- this.update$acc
      att.xi <- this.update$att
      MH.xi  <- this.update$MH

      this.update <- mhUpdate(acc = acc.bw.gp, att = att.bw.gp,
                              nattempts = bw.gp.attempts,
                              MH = MH.bw.gp, target.min = 0.3, target.max = 0.6,
                              lower = 0.8, higher = 1.2)
      acc.bw.gp <- this.update$acc
      att.bw.gp <- this.update$att
      MH.bw.gp  <- this.update$MH
    }

    ####################################################
    ##############      Impute missing      ############
    ####################################################
    # if (any(miss)) {
    #   y.tmp <- y
    #   for (t in missing.times) {
    #     # calculate mu and sigma
    #     miss.t   <- miss[, t]
    #     sigma    <- exp(logsig[miss.t, t])
    #     mu_star  <- mu[miss.t, t] + sigma * (theta[miss.t, t]^xi - 1) / xi
    #     sig_star <- alpha * sigma * theta[miss.t, t]^xi
    #     xi_star  <- alpha * xi
    #     y.tmp[miss.t, t] <- revd(n = sum(miss.t), loc = mu_star,
    #                              scale = sig_star, shape = xi_star,
    #                              type = "GEV")
    #   }
    # }

    #KEEP TRACK OF STUFF:
    keep.bw.basis[iter]  <- bw.basis
    # print(bw.basis)
    keep.mu[iter, , ]    <- mu[these.sites, these.days]
    keep.logsig[iter, , ] <- logsig[these.sites, these.days]
    keep.tau1[iter, ]    <- tau1
    keep.tau2[iter, ]    <- tau2
    keep.beta1[iter, ]   <- beta1
    keep.beta2[iter, ]   <- beta2
    keep.xi[iter]        <- xi
    keep.A[iter, , ]     <- A
    # keep.beta.mu[iter, ] <- c(beta1.mu, beta2.mu)
    keep.beta.sd[iter, ] <- c(beta1.sd, beta2.sd)
    keep.bw.gp[iter]     <- bw.gp
    # if (any(miss)) {
    #   keep.y[iter, ]     <- y.tmp[miss]
    # }
    if (iter > burn) {
      theta.mn <- theta.mn + theta / (iters - burn)
    }


    #DISPLAY CURRENT VALUE:

    if (iter %% update == 0) {
      if (iterplot) {
        acc.rate.mu     <- round(acc.mu / att.mu, 3)
        acc.rate.logsig <- round(acc.logsig / att.logsig, 3)
        acc.rate.xi     <- round(acc.xi / att.xi, 3)
        acc.rate.bw.basis <- round(acc.bw.basis / att.bw.basis, 3)
        acc.rate.bw.gp    <- round(acc.bw.gp / att.bw.gp, 3)

        if (iter > burn) {
          start <- burn + 1
        } else {
          start <- max(iter - 2000, 1)
        }

        if (!exists("this.plot")) {
          this.plot <- "all"
        }
        # this.plot <- "sig"

        if (this.plot == "mu") {
          if (L <= 9) {
            par(mfrow = c(3, 3))
            for (plot.idx in 1:min(p.mu, 7)) {
              plot(keep.beta1[start:iter, plot.idx],
                   main = bquote(paste(mu, ": ", beta[.(plot.idx)])),
                   xlab = acc.rate.mu[plot.idx], type = "l")#,
                   # ylab = paste("MH =", round(MH.beta1[plot.idx], 3)))
            }
          } else {
            par(mfrow = c(3, 4))
            for (plot.idx in 1:min(p.mu, 10)) {
              plot(keep.beta1[start:iter, plot.idx],
                   main = bquote(paste(mu, ": ", beta[.(plot.idx)])),
                   xlab = acc.rate.mu[plot.idx], type = "l")#,
                   #ylab = paste("MH =", round(MH.beta1[plot.idx], 3)))
            }
          }

          # plot(keep.beta.mu[start:iter, 1],
          #      main = bquote(paste(mu, ": ", mu[beta])), type = "l")

          plot(keep.beta.sd[start:iter, 1],
               main = bquote(paste(mu, ": ", sigma[beta])), type = "l")
          this.plot <- "sig"
        } else if (this.plot == "sig") {
          if (L <= 9) {
            par(mfrow = c(3, 3))
            for (plot.idx in 1:min(p.sig, 7)) {
              plot(keep.beta2[start:iter, plot.idx],
                   main = bquote(paste(sigma, ": ", beta[.(plot.idx)])),
                   xlab = acc.rate.logsig[plot.idx], type = "l")#,
                   #ylab = paste("MH =", round(MH.beta2[plot.idx], 3)))
            }
          } else {
            par(mfrow = c(3, 4))
            for (plot.idx in 1:min(p.sig, 10)) {
              plot(keep.beta2[start:iter, plot.idx],
                   main = bquote(paste(sigma, ": ", beta[.(plot.idx)])),
                   xlab = acc.rate.logsig[plot.idx], type = "l")#,
                   #ylab = paste("MH =", round(MH.beta2[plot.idx], 3)))
            }
          }

          # plot(keep.beta.mu[start:iter, 2],
          #      main = bquote(paste(sigma, ": ", mu[beta])), type = "l")

          plot(keep.beta.sd[start:iter, 2],
               main = bquote(paste(sigma, ": ", sigma[beta])), type = "l")
          this.plot <- "all"
          # this.plot <- "mu"
        } else {
          par(mfrow = c(3, 5))

          for (i in 1:2) {
            plot(keep.beta1[start:iter, i],
                 main = bquote(paste(mu, ": ", beta[.(i)])),
                 xlab = acc.rate.mu[i], type = "l")#,
                 #ylab = paste("MH =", round(MH.beta1[i], 3)))
          }

          for (i in 1:2) {
            plot(keep.mu[start:iter, i, i], type = "l",
                 main = bquote(paste(mu[.(i)])))
          }

          # plot(keep.beta1[start:iter, p.mu],
          #      main = bquote(paste(mu, ": ", beta[.(p.mu)])),
          #      xlab = acc.rate.mu[p.mu], type = "l",
          #      ylab = paste("MH =", round(MH.beta1[p.mu], 3)))

          plot(keep.xi[start:iter], main = bquote(xi),
               xlab = acc.rate.xi, type = "l",
               ylab = paste("MH =", round(MH.xi, 3)))

          for (i in 1:2) {
            plot(keep.beta2[start:iter, i],
                 main = bquote(paste(sigma, ": ", beta[.(i)])),
                 xlab = acc.rate.logsig[i], type = "l")#,
                 # ylab = paste("MH =", round(MH.beta2[i], 3)))
          }

          for (i in 1:2) {
            plot(keep.logsig[start:iter, i, i], type = "l",
                 main = bquote(paste("log(", sigma[.(i)], ")")))
          }

          # plot(keep.beta2[start:iter, p.sig],
          #      main = bquote(paste(sigma, ": ", beta[.(p.sig)])),
          #      xlab = acc.rate.logsig[p.sig], type = "l",
          #      ylab = paste("MH =", round(MH.beta2[p.sig], 3)))

          plot(keep.bw.basis[start:iter], main = "kernel bandwidth",
               xlab = acc.rate.bw.basis, type = "l",
               ylab = paste("MH =", round(MH.bw.basis, 3)))
          plot(keep.tau1[start:iter, 1], main = "tau1[1]", type = "l")
          plot(keep.tau1[start:iter, 2], main = "tau1[2]", type = "l")
          plot(keep.tau2[start:iter, 1], main = "tau2[1]", type = "l")
          plot(keep.tau2[start:iter, 2], main = "tau2[2]", type = "l")
          # plot(log(keep.A[start:iter, 1, 1]), main = "log(A[1, 1])", type = "l")
          #
          # plot(log(keep.A[start:iter, 2, 1]), main = "log(A[2, 1])", type = "l")
          #
          # plot(log(keep.A[start:iter, L, 1]), main = "log(A[L, 1])", type = "l")
          #
          # plot(log(keep.A[start:iter, L, 2]), main = "log(A[L, 2])", type = "l")

          # plot(log(keep.A[start:iter, L, nt]),
          #      main = "log(A[L, nt])", type = "l")

          plot(keep.bw.gp[start:iter], main = "gaussian process bandwidth",
               xlab = acc.rate.bw.gp, type = "l",
               ylab = paste("MH =", round(MH.bw.gp, 3)))

          # this.plot <- "mu"
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
       bw.basis = keep.bw.basis[return.iters],
       xi = keep.xi[return.iters],
       # betamu = keep.beta.mu[return.iters, , drop = FALSE],
       # betasd = keep.beta.sd[return.iters, , drop = FALSE],
       bw.gp = keep.bw.gp[return.iters],
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

  theta.xi  <- theta^xi
  sigma    <- exp(logsig)
  mu_star  <- mu + sigma * ((theta.xi) - 1) / xi
  sig_star <- alpha * sigma * (theta.xi)
  xi_star  <- alpha * xi

  tx <- (1 + xi_star * (y - mu_star) / sig_star)^(-1 / xi_star)
  ll <- -tx + (y > thresh) * ((xi_star + 1) * log(tx) - log(sig_star))
  ll[is.na(y)] <- 0
  ll[is.na(ll)] <- -Inf

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

  this.grad <- -(tx_star^(-1 / xi_star - 1)) / sig_star +
    (y > thresh) * (xi_star + 1) / (sig_star * tx_star)
  this.grad[is.na(y)] <- 0  # for the missing data

  grad <- rep(0, p.mu)
  for (j in 1:p.mu) {
    grad[j] <- sum(this.grad * X.mu[, , j])
  }

  return(grad)
}

grad_logpost_betamu <- function(beta1, beta.mu, beta.sd, X.mu, y, theta, logsig,
                                xi, thresh, alpha) {
  mu <- 0
  p.mu <- dim(X.mu)[3]
  for (j in 1:p.mu) {
    mu <- mu + X.mu[, , j] * beta1[j]
  }

  sig <- exp(logsig)

  tx  <- (1 + xi * (y - mu) / sig)

  this.grad <- - (tx^(-1 / (alpha * xi) - 1)) * theta^(1 / alpha) /
    (alpha * sig) +
    (y > thresh) * (alpha * xi + 1) / (alpha * tx * sig)
  this.grad[is.na(y)] <- 0  # for the missing data

  grad <- rep(0, p.mu)
  for (j in 1:p.mu) {
    grad[j] <- sum(this.grad * X.mu[, , j]) - (beta1[j] - beta.mu) / beta.sd^2
  }

  return(grad)
}

hess_logpost_betamu <- function(beta1, beta.mu, beta.sd, X.mu, y, theta, logsig,
                                xi, thresh, alpha) {
  mu <- 0
  p.mu <- dim(X.mu)[3]
  for (j in 1:p.mu) {
    mu <- mu + X.mu[, , j] * beta1[j]
  }

  sigma <- exp(logsig)
  tx <- (1 + xi * (y - mu) / sigma)

  this.hess <- - (tx^(-1 / (alpha * xi) - 2) * (alpha * xi + 1) *
                    theta^(1 / alpha)) / (alpha^2 * sigma^2) +
    (y > thresh) * (alpha * xi + 1) * xi / (alpha * tx^2 * sigma^2)
  this.hess[is.na(y)] <- 0  # for missing data

  hess <- matrix(0, p.mu, p.mu)
  for (j in 1:p.mu) {
    for (i in j:p.mu) {
      hess[i, j] <- hess[j, i] <- sum(this.hess * X.mu[, , i] * X.mu[, , j])
    }
  }

  hess <- hess - diag(1 / beta.sd^2, p.mu)  # account for prior

  return(hess)
}

hess_loglike_betamu <- function(beta1, X.mu, y, theta, logsig, xi, thresh,
                                alpha) {
  mu <- 0
  p.mu <- dim(X.mu)[3]
  for (j in 1:p.mu) {
    mu <- mu + X.mu[, , j] * beta1[j]
  }

  # print(mean(mu))
  sigma    <- exp(logsig)
  mu_star  <- mu + sigma * ((theta^xi) - 1) / xi
  sig_star <- alpha * sigma * (theta^xi)
  xi_star  <- alpha * xi

  tx_star  <- (1 + xi_star * (y - mu_star) / sig_star)
  this.hess <- -(xi_star + 1) * tx_star^(-1 / xi_star - 2) / sig_star^2 +
    (y > thresh) * (xi_star + 1) * xi_star / (sig_star^2 * tx_star^2)
  this.hess[is.na(y)] <- 0

  hess <- matrix(0, p.mu, p.mu)
  for (i in 1:p.mu) {
    for (j in i:p.mu) {
      hess[i, j] <- hess[j, i] <- sum(this.hess * X.mu[, , i] * X.mu[, , j])
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

  return(sum(ll[!is.na(y)]))
}

logpost_mu <- function(beta1, beta.mu, beta.sd, X.mu, y, theta, logsig, xi,
                       thresh, alpha) {
  mu <- 0
  p.mu <- dim(X.mu)[3]
  for (j in 1:p.mu) {
    mu <- mu + X.mu[, , j] * beta1[j]
  }

  ll <- loglike(y = y, theta = theta, mu = mu, logsig = logsig, xi = xi,
                thresh = thresh, alpha = alpha)

  lp <- sum(dnorm(beta1, beta.mu, beta.sd, log = TRUE))

  return(sum(ll[!is.na(y)]) + lp)
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
  this.grad <- (-tx_star^(-1 / xi_star - 1) *
                  xi * (y - mu) / (xi_star * sigma * theta^xi) +
                  (y > thresh) * (-1 + ((xi_star + 1) * xi * (y - mu)) /
                                    (xi_star * sigma * tx_star * theta^xi)))
  this.grad[is.na(y)] <- 0

  grad <- rep(0, p.sig)
  for (j in 1:p.sig) {
    grad[j] <- sum(this.grad * X.sig[, , j])
  }

  return(grad)
}

grad_logpost_betasig <- function(beta2, beta.mu, beta.sd, X.sig, y, theta, mu,
                                 xi, thresh, alpha) {
  logsig <- 0
  p.sig <- dim(X.sig)[3]
  for (j in 1:p.sig) {
    logsig <- logsig + X.sig[, , j] * beta2[j]
  }

  sigma <- exp(logsig)
  res <- y - mu

  tx <- (1 + xi * res / sigma)

  this.grad <- - theta^(1 / alpha) * res * tx^(-1 / (alpha * xi) - 1) /
    (alpha * sigma) +
    (y > thresh) * ((alpha * xi + 1) * res / (alpha * tx * sigma) - 1)

  this.grad[is.na(y)] <- 0

  grad <- rep(0, p.sig)
  for (j in 1:p.sig) {
    grad[j] <- sum(this.grad * X.sig[, , j]) - (beta2[j] - beta.mu) / beta.sd^2
  }

  return(grad)
}

hess_logpost_betasig <- function(beta2, beta.mu, beta.sd, X.sig, y, theta, mu,
                                 xi, thresh, alpha) {
  logsig <- 0
  p.sig <- dim(X.sig)[3]
  for (j in 1:p.sig) {
    logsig <- logsig + X.sig[, , j] * beta2[j]
  }

  sigma <- exp(logsig)
  res <- y - mu

  tx <- (1 + xi * res / sigma)

  this.hess <- - theta^(1 / alpha) * res * tx^(-1 / (alpha * xi) - 1) * (
    (alpha * xi + 1) * xi * res / (alpha * xi * tx * sigma) - 1
  ) / (alpha * sigma) -
    (y > thresh) * (alpha * xi + 1) * res * sigma /
    (alpha * (sigma + xi * res)^2)

  hess <- matrix(0, p.sig, p.sig)
  for (i in 1:p.sig) {
    for (j in i:p.sig) {
      hess[i, j] <- hess[j, i] <- sum(this.hess * X.sig[, , i] * X.sig[, , j])
    }
  }

  hess <- hess - diag(1 / beta.sd^2, p.sig)  # account for prior

  return(hess)
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
  d12 <- (-1 / xi_star - 1) * tx_star^(-1 / xi_star - 2) * dtx^2
  d21 <- tx_star^(-1 / xi_star - 1) * (-dtx)
  this.hess <- (d12 + d21) / xi_star
  this.hess <- this.hess + (y > thresh) * (xi_star + 1) / alpha *
    (y - mu) / theta^xi * (-1 / (sigma * tx_star) - dtx / (tx_star^2 * sigma))
  this.hess[is.na(y)] <- 0
  for (i in 1:p.sig) {
    for (j in i:p.sig) {
      hess[i, j] <- hess[j, i] <- sum(this.hess * X.sig[, , i] * X.sig[, , j])
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

logpost_sig <- function(beta2, beta.mu, beta.sd, X.sig, y, theta, mu, xi,
                        thresh, alpha) {
  logsig <- 0
  p.sig <- dim(X.sig)[3]
  for (j in 1:p.sig) {
    logsig <- logsig + X.sig[, , j] * beta2[j]
  }

  ll <- loglike(y = y, theta = theta, mu = mu, logsig = logsig, xi = xi,
                thresh = thresh, alpha = alpha)

  lp <- sum(dnorm(beta2, beta.mu, beta.sd, log = TRUE))
  return(sum(ll[!is.na(y)]) + lp)
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



# # trying a Langevin block update
# cur.grad <- grad_loglike_betamu(beta1 = beta1, X.mu = X.mu, y = y,
#                                 theta = theta, logsig = logsig, xi = xi,
#                                 thresh = thresh, alpha = alpha)
# cur.hess <- hess_loglike_betamu(beta1 = beta1, X.mu = X.mu, y = y,
#                                 theta = theta, logsig = logsig, xi = xi,
#                                 thresh = thresh, alpha = alpha)
#
# cat(cur.hess)
# cur.mn  <- beta1 + 0.5 * MH.beta1^2 * cur.grad
# cur.var <- tryCatch(expr = solve(cur.hess), error = function(e) 0)
# cur.var <- diag(MH.beta1^2)
# canb <- cur.mn + t(chol(cur.var)) %*% rnorm(p.mu, 0, 1)
#
# can.grad <- grad_loglike_betamu(beta1 = canb, X.mu = X.mu, y = y,
#                                 theta = theta, logsig = logsig, xi = xi,
#                                 thresh = thresh, alpha = alpha)
# can.hess <- hess_loglike_betamu(beta1 = canb, X.mu = X.mu, y = y,
#                                 theta = theta, logsig = logsig, xi = xi,
#                                 thresh = thresh, alpha = alpha)
#
#
#
# R <- sum(canll - curll) +
#       dnorm(canb, 0, beta1.sd, log = TRUE) -
#       dnorm(beta1, 0, beta1.sd, log = TRUE) +
#       dnorm(beta1, can.mn, MH.beta1[j], log = TRUE) -
#       dnorm(canb, cur.mn, MH.beta1[j], log = TRUE)


# for (j in 1:p.mu) {
#   att.beta1[j] <- att.beta1[j] + 1
#   cur.mn  <- beta1[j] + 0.5 * MH.beta1[j]^2 * cur.grad[j]
#   # cur.hess <- hess_loglike_betamu(beta1 = beta1, X.mu = X.mu, y = Y,
#   #                                 theta = theta, logsig = logsig, xi = xi,
#   #                                 thresh = thresh, alpha = alpha)
#
#   canb    <- beta1
#   canb[j] <- rnorm(1, cur.mn, MH.beta1[j])
#
#   # print(paste("canb:", canb[j]))
#   canmu <- mu + X.mu[, , j] * (canb[j] - beta1[j])
#   canll <- loglike(y, theta, canmu, logsig, xi, thresh, alpha)
#
#   # need to adjust the candidate distribution for acc.ratio
#   can.grad <- grad_loglike_betamu(beta1 = canb, X.mu = X.mu, y = y,
#                                   theta = theta, logsig = logsig, xi = xi,
#                                   thresh = thresh, alpha = alpha)
#   can.mn <- canb[j] + 0.5 * MH.beta1[j]^2 * can.grad[j]
#
#   R <- sum(canll - curll) +
#     dnorm(canb[j], 0, beta1.sd, log = TRUE) -
#     dnorm(beta1[j], 0, beta1.sd, log = TRUE) +
#     dnorm(beta1[j], can.mn, MH.beta1[j], log = TRUE) -
#     dnorm(canb[j], cur.mn, MH.beta1[j], log = TRUE)
#
#   if (!is.na(R)) { if (log(runif(1)) < R) {
#     acc.beta1[j] <- acc.beta1[j] + 1
#     beta1[j]     <- canb[j]
#     mu           <- canmu
#     cur.grad[j]  <- can.grad[j]
#     curll        <- canll
#   }}
#
# }

# # trying a Langevin block update
# cur.grad <- grad_loglike_betasig(beta2 = beta2, X.sig = X.sig, y = y,
#                                  theta = theta, mu = mu, xi = xi,
#                                  thresh = thresh, alpha = alpha)
# for (j in 1:p.mu) {
#   att.beta2[j] <- att.beta2[j] + 1
#   cur.mn  <- beta2[j] + 0.5 * MH.beta2[j]^2 * cur.grad[j]
#   # cur.hess <- hess_loglike_betamu(beta1 = beta1, X.mu = X.mu, y = Y,
#   #                                 theta = theta, logsig = logsig, xi = xi,
#   #                                 thresh = thresh, alpha = alpha)
#
#   canb    <- beta2
#   canb[j] <- rnorm(1, cur.mn, MH.beta2[j])
#
#   # print(paste("canb:", canb[j]))
#   canlogs <- logsig + X.sig[, , j] * (canb[j] - beta2[j])
#   canll   <- loglike(y, theta, mu, canlogs, xi, thresh, alpha)
#
#   # need to adjust the candidate distribution for acc.ratio
#   can.grad <- grad_loglike_betasig(beta2 = canb, X.sig = X.sig, y = y,
#                                    theta = theta, mu = mu, xi = xi,
#                                    thresh = thresh, alpha = alpha)
#   can.mn <- canb[j] + 0.5 * MH.beta2[j]^2 * can.grad[j]
#
#   R <- sum(canll - curll) +
#     dnorm(canb[j], 0, beta2.sd, log = TRUE) -
#     dnorm(beta2[j], 0, beta2.sd, log = TRUE) +
#     dnorm(beta2[j], can.mn, MH.beta2[j], log = TRUE) -
#     dnorm(canb[j], cur.mn, MH.beta2[j], log = TRUE)
#
#   if (!is.na(R)) { if (log(runif(1)) < R) {
#     acc.beta2[j] <- acc.beta2[j] + 1
#     beta2[j]     <- canb[j]
#     logsig       <- canlogs
#     cur.grad[j]  <- can.grad[j]
#     curll        <- canll
#   }}
#
# }


#
# else if (FALSE) {  # block update with normal approximation
#   att.beta1 <- att.beta1 + 1
#   VVV.cur <- solve(-beta1.hess.cur)
#   mean.cur <- beta1 + VVV.cur %*% beta1.grad.cur
#
#   canb  <- mean.cur + t(chol(VVV.cur)) %*% rnorm(p.mu)
#   canmu <- 0
#   for (j in 1:p.mu) {
#     canmu   <- canmu + X.mu[, , j] * canb[j]
#   }
#   if (xi < 0 & any(y - canmu > -exp(logsig) / xi, na.rm = TRUE)) {
#     R <- -Inf
#   } else if (xi > 0 & any(y - canmu < -exp(logsig) / xi, na.rm = TRUE)) {
#     R <- -Inf
#   } else {
#     canll    <- loglike(y, theta, canmu, logsig, xi, thresh, alpha)
#     # canll    <- loglike(y, theta.xi, canmu, logsig, xi, thresh, alpha)
#
#     # get the adjustments for the ratio of the update
#     beta1.grad.can <- grad_logpost_betamu(beta1 = canb, beta.mu = beta1.mu,
#                                           beta.sd = beta1.sd, X.mu = X.mu,
#                                           y = y, theta = theta,
#                                           logsig = logsig, xi = xi,
#                                           thresh = thresh, alpha = alpha)
#     beta1.hess.can <- hess_logpost_betamu(beta1 = canb, beta.mu = beta1.mu,
#                                           beta.sd = beta1.sd, X.mu = X.mu,
#                                           y = y, theta = theta,
#                                           logsig = logsig, xi = xi,
#                                           thresh = thresh, alpha = alpha)
#
#     VVV.can <- solve(-beta1.hess.can)
#     mean.can <- canb + VVV.can %*% beta1.grad.can
#
#     prop.cur <- 0.5 * log(det(-beta1.hess.can)) -
#       0.5 * t(beta1 - mean.can) %*% (-beta1.hess.can) %*% (beta1 - mean.can)
#     prop.can <- 0.5 * log(det(-beta1.hess.cur)) -
#       0.5 * t(canb - mean.cur) %*% (-beta1.hess.cur) %*% (canb - mean.cur)
#
#     R <- sum(canll - curll) +
#       sum(dnorm(canb, beta1.mu, beta1.sd, log = TRUE)) -
#       sum(dnorm(beta1, beta1.mu, beta1.sd, log = TRUE)) +
#       prop.cur - prop.can
#   }
#   if (!is.nan(R)) { if (log(runif(1)) < R) {
#     acc.beta1   <- acc.beta1 + 1
#     beta1       <- canb
#     beta1.grad.cur <- beta1.grad.can
#     beta1.hess.cur <- beta1.hess.can
#     mu             <- canmu
#     curll          <- canll
#   } }
# }