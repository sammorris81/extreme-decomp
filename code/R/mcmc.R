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

ReShMCMC<-function(y, X, X1 = NULL, X2 = NULL, s, knots, thresh, B, alpha,
                   beta1 = NULL, beta1.mu = 0, beta1.sd = 10,
                   beta2 = NULL, beta2.mu = 0, beta2.sd = 1,
                   xi = 0.001, A = NULL,
                   can.mu.sd = 0.1, can.ls.sd = 0.1,
                   beta1.tau.a = 0.1, beta1.tau.b = 0.1,
                   beta2.tau.a = 0.1, beta2.tau.b = 0.1,
                   tau1.a = 0.1, tau1.b = 0.1, tau1.fix = FALSE,
                   tau2.a = 0.1, tau2.b = 0.1, tau2.fix = FALSE,
                   mu.attempts = 50, ls.attempts = 50, xi.attempts = 50,
                   bw.basis.init = NULL, bw.basis.random = TRUE,
                   bw.basis.attempts = 50,
                   bw.gp.init = NULL, bw.gp.attempts = 50,
                   time.interact = FALSE,
                   keep.sites = NULL, keep.days = NULL, keep.knots = NULL,
                   keep.burn = FALSE, iters = 5000, burn = 1000, update = 10,
                   iterplot = FALSE){
  require(extRemes)
  require(fields)
  # BOOKKEEPING

  ns   <- nrow(y)
  nt   <- ncol(y)
  L    <- ncol(B)
  if (is.null(keep.sites)) {
    keep.sites <- ns
  }
  if (is.null(keep.days)) {
    keep.days <- nt
  }
  if (is.null(keep.knots)) {
    keep.knots <- L
  } else {
    keep.knots <- min(keep.knots, L)
  }

  if (is.null(X1)) {
    X1 <- X
  }
  if (is.null(X2)) {
    X2 <- X
  }

  # initialize distance matrix from site to knot
  dw2 <- rdist(s, knots)^2
  dw2[dw2 < 1e-4] <- 0

  if (is.null(bw.basis.init)) {
    bw.basis <- quantile(dw2, 0.1)
  } else {
    bw.basis <- bw.basis.init
  }

  bw.basis.min <- quantile(dw2, 0.005)
  bw.basis.max <- quantile(dw2, 0.995)

  B.X <- makeW(dw2 = dw2, rho = bw.basis)
  X1 <- add.basis.X(X = X1, B = B.X, time.interact = time.interact)
  X2 <- add.basis.X(X = X2, B = B.X, time.interact = time.interact)

  p1 <- dim(X1)[3]
  p2 <- dim(X2)[3]

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
    beta2    <- rep(0, p2)
    beta2[1] <- log(sqrt(6) * sd(as.vector(y), na.rm = TRUE) / pi)
  }

  if (is.null(beta1)) {
    beta1    <- rep(0, p1)
    beta1[1] <- -0.57722 * beta2[1]
  }

  # set the initial mu and logsig to the mean
  mu <- Xb1 <- getXBeta(X = X1, beta = beta1)
  ls <- Xb2 <- getXBeta(X = X2, beta = beta2)

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
    tau2[t] <- 1 / var(ls[, t])
  }

  # adjust xi until it is valid for all y, mu, and sig
  if (xi < 0) {
    while (any(y - mu > -exp(ls) / xi, na.rm = TRUE)) {
      xi <- xi * 0.8
    }
  } else if (xi > 0) {
    while (any(y - mu < -exp(ls) / xi, na.rm = TRUE)) {
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
  curll <- loglike(y, theta, mu, ls, xi, thresh, alpha)

  # STORAGE:
  these.sites   <- sort(sample(1:ns, keep.sites))  # which sites to keep GP
  these.days    <- sort(sample(1:nt, keep.days))   # which days to keep GP and A
  these.knots   <- sort(sample(1:L, keep.knots))   # which knots to keep A
  keep.beta1    <- matrix(0, iters, p1)
  keep.beta2    <- matrix(0, iters, p2)
  keep.xi       <- rep(0, iters)
  keep.beta.sd  <- matrix(0, iters, 2)  # sd terms for beta priors
  keep.tau1     <- matrix(0, iters, nt)
  keep.tau2     <- matrix(0, iters, nt)
  keep.mu       <- array(0, dim = c(iters, keep.sites, keep.days))
  keep.ls       <- array(0, dim = c(iters, keep.sites, keep.days))
  keep.bw.basis <- rep(0, iters)  # bandwidth term for basis functions
  keep.bw.gp    <- rep(0, iters)  # bandwidth for GP
  keep.A        <- array(0, dim = c(iters, keep.knots, keep.days))
  if (any(miss)) {
    keep.y <- matrix(0, iters, sum(miss))  # only record the missing data
  } else {
    keep.y <- NULL
  }

  # TUNING:

  cuts <- exp(c(-1, 0, 1, 2, 5, 10))
  MH.a  <- rep(1, 100)
  att.a <- acc.a <- 0 * MH.a
  att.mu <- acc.mu <- MH.mu <- matrix(can.mu.sd, ns, nt)
  att.ls <- acc.ls <- MH.ls <- matrix(can.ls.sd, ns, nt)
  att.xi <- acc.xi    <- MH.xi    <- 0.1
  att.bw.basis <- acc.bw.basis <- MH.bw.basis <- 0.001
  att.bw.gp    <- acc.bw.gp    <- MH.bw.gp    <- 0.05

  tic <- proc.time()[3]
  for (iter in 1:iters) {
    ####################################################
    ##############      Random effects A    ############
    ####################################################
    oldA  <- A
    this.update <- updateA(A = A, cuts = cuts, bins = bins, Ba = Ba,
                           theta = theta, y = y, mu = mu, ls = ls,
                           xi = xi, thresh = thresh, alpha = alpha,
                           curll = curll, MH = MH.a)

    A     <- this.update$A
    l1    <- this.update$l1
    theta <- this.update$theta
    curll <- this.update$curll

    ####################################################
    ##########      bandwidth for kernels      #########
    ####################################################
    if (bw.basis.random) {
      this.update <- updateXBasisBW(bw = bw.basis, bw.min = bw.basis.min,
                                    bw.max = bw.basis.max, bw.mn = -2,
                                    bw.sd = 1,
                                    X1 = X1, beta1 = beta1, Xb1 = Xb1, mu = mu,
                                    tau1 = tau1, SS1 = SS1,
                                    X2 = X2, beta2 = beta2, Xb2 = Xb2, ls = ls,
                                    tau2 = tau2, SS2 = SS2,
                                    Qb = Qb, dw2 = dw2,
                                    time.interact = time.interact,
                                    acc = acc.bw.basis, att = att.bw.basis,
                                    MH = MH.bw.basis)

      bw.basis     <- this.update$bw
      X1           <- this.update$X1
      Xb1          <- this.update$Xb1
      SS1          <- this.update$SS1
      X2           <- this.update$X2
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
                            y = y, theta = theta, ls = ls, xi = xi,
                            thresh = thresh, alpha = alpha, SS = SS1,
                            curll = curll, acc = acc.mu, att = att.mu,
                            MH = MH.mu)
    mu     <- this.update$mu
    SS1    <- this.update$SS
    curll  <- this.update$curll
    acc.mu <- this.update$acc
    att.mu <- this.update$att

    # logsig
    this.update <- updateLS(ls = ls, Qb = Qb, tau = tau2, Xb = Xb2,
                            y = y, theta = theta, mu = mu, xi = xi,
                            thresh = thresh, alpha = alpha, SS = SS2,
                            curll = curll, acc = acc.ls, att = att.ls,
                            MH = MH.ls)
    ls     <- this.update$ls
    SS2    <- this.update$SS
    curll  <- this.update$curll
    acc.ls <- this.update$acc
    att.ls <- this.update$att

    # xi
    this.update <- updateXi(xi = xi, y = y, mu = mu, ls = ls,
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
                                param = mu, X = X1, SS = SS1, tau = tau1)
    beta1 <- this.update$beta
    Xb1   <- this.update$Xb
    SS1   <- this.update$SS

    this.update <- updateGPBeta(beta = beta2, beta.sd = beta2.sd, Qb = Qb,
                                param = ls, X = X2, SS = SS2, tau = tau2)
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
                               tau.a = tau1.a, tau.b = tau1.b, ns = ns)
    tau1 <- this.update$tau

    this.update <- updateGPTau(tau = tau2, SS = SS2,
                               tau.a = tau2.a, tau.b = tau2.b, ns = ns)
    tau2 <- this.update$tau

    # spatial range for GP
    this.update <- updateGPBW(bw = bw.gp, bw.min = bw.gp.min,
                              bw.mn = -2, bw.sd = 1, bw.max = bw.gp.max,
                              Qb = Qb, d = d,
                              mu = mu, Xb1 = Xb1, tau1 = tau1, SS1 = SS1,
                              ls = ls, Xb2 = Xb2, tau2 = tau2, SS2 = SS2,
                              acc = acc.bw.gp, att = att.bw.gp, MH = MH.bw.gp)
    bw.gp     <- this.update$bw
    Qb        <- this.update$Qb
    SS1       <- this.update$SS1
    SS2       <- this.update$SS2
    acc.bw.gp <- this.update$acc
    att.bw.gp <- this.update$att

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

      this.update <- mhUpdate(acc = acc.mu, att = att.mu, MH = MH.mu,
                              nattempts = 50,
                              target.min = 0.3, target.max = 0.6,
                              lower = 0.8, higher = 1.2)
      acc.mu <- this.update$acc
      att.mu <- this.update$att
      MH.mu  <- this.update$MH

      this.update <- mhUpdate(acc = acc.ls, att = att.ls, MH = MH.ls,
                              nattempts = 50,
                              target.min = 0.3, target.max = 0.6,
                              lower = 0.8, higher = 1.2)
      acc.ls <- this.update$acc
      att.ls <- this.update$att
      MH.ls  <- this.update$MH

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
    keep.beta1[iter, ]   <- beta1
    keep.beta2[iter, ]   <- beta2
    keep.xi[iter]        <- xi
    keep.beta.sd[iter, ] <- c(beta1.sd, beta2.sd)
    keep.tau1[iter, ]    <- tau1
    keep.tau2[iter, ]    <- tau2
    keep.mu[iter, , ]    <- mu[these.sites, these.days]
    keep.ls[iter, , ]    <- ls[these.sites, these.days]
    keep.bw.basis[iter]  <- bw.basis
    keep.bw.gp[iter]     <- bw.gp
    keep.A[iter, , ]     <- A[these.knots, these.days]
    # if (any(miss)) {
    #   keep.y[iter, ]     <- y.tmp[miss]
    # }

    #DISPLAY CURRENT VALUE:

    if (iter %% update == 0) {
      if (iterplot) {
        acc.rate.mu     <- round(acc.mu / att.mu, 3)
        acc.rate.ls     <- round(acc.ls / att.ls, 3)
        acc.rate.xi     <- round(acc.xi / att.xi, 3)
        acc.rate.bw.basis <- round(acc.bw.basis / att.bw.basis, 3)
        acc.rate.bw.gp    <- round(acc.bw.gp / att.bw.gp, 3)

        if (iter > burn) {
          start <- burn + 1
        } else {
          start <- max(iter - 2000, 1)
        }

        par(mfrow = c(3, 5))

        for (i in 1:2) {
          plot(keep.beta1[start:iter, i],
               main = bquote(paste(mu, ": ", beta[.(i)])))
        }

        for (i in 1:2) {
          this.site <- these.sites[i]
          this.day  <- these.days[i]
          plot(keep.mu[start:iter, this.site, this.day], type = "l",
               main = bquote(paste(mu[.(i)])))
        }

        plot(keep.xi[start:iter], main = bquote(xi),
             xlab = acc.rate.xi, type = "l",
             ylab = paste("MH =", round(MH.xi, 3)))

        for (i in 1:2) {
          plot(keep.beta2[start:iter, i],
               main = bquote(paste(sigma, ": ", beta[.(i)])))
        }

        for (i in 1:2) {
          this.site <- these.sites[i]
          this.day  <- these.days[i]
          plot(keep.ls[start:iter, this.site, this.day], type = "l",
               main = bquote(paste("log(", sigma[.(i)], ")")))
        }

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

      }
    }
    cat("    Finished fit:", iter, "of", iters, "iters \n")
  } #end iter

  toc <- proc.time()[3]

  if (keep.burn) {
    return.iters <- 1:iters
  } else {
    return.iters <- (burn + 1):iters
  }

  if (any(miss)) {
    keep.y <- keep.y[return.iters, , drop = FALSE]
  }

  keep.beta1[iter, ]   <- beta1
  keep.beta2[iter, ]   <- beta2
  keep.xi[iter]        <- xi
  keep.beta.sd[iter, ] <- c(beta1.sd, beta2.sd)
  keep.tau1[iter, ]    <- tau1
  keep.tau2[iter, ]    <- tau2
  keep.mu[iter, , ]    <- mu[these.sites, these.days]
  keep.ls[iter, , ]    <- ls[these.sites, these.days]
  keep.bw.basis[iter]  <- bw.basis
  keep.bw.gp[iter]     <- bw.gp
  keep.A[iter, , ]     <- A[these.knots, these.days]
  # if (any(miss)) {
  #   keep.y[iter, ]     <- y.tmp[miss]
  # }

  list(beta1   = keep.beta1[return.iters, , drop = FALSE],
       beta2   = keep.beta2[return.iters, , drop = FALSE],
       xi      = keep.xi[return.iters],
       beta.sd = keep.beta.sd[return.iters, ],
       tau1    = keep.tau1[return.iters, ],
       tau2    = keep.tau2[return.iters, ],
       mu      = keep.mu[return.iters, , , drop = FALSE],
       ls      = keep.ls[return.iters, , , drop = FALSE],
       bw.basis = keep.bw.basis[return.iters],
       bw.gp = keep.bw.gp[return.iters],
       A = keep.A[return.iters, , , drop = FALSE],
       # y.pred = keep.y,
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