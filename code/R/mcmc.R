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

ReShMCMC<-function(y, s, thresh, B, alpha,
                   # details for the GP prior on mu and logsig
                   beta.int = NULL, beta.int.mn = NULL,
                   beta.int.attempts = 100, canbeta.int.sd = 0.1,
                   beta.time = NULL, beta.time.mn = NULL,
                   beta.time.attempts = 100, canbeta.time.sd = 0.1,
                   # this is the prior standard deviation on the mean of the
                   # GPs
                   mu.beta.pri.sd = 100, ls.beta.pri.sd = 10,
                   xi = 0.001, xi.min = -0.5, xi.max = 0.5,
                   xi.mn = 0, xi.sd = 0.5, xi.attempts = 50, canxi.sd = 0.1,
                   tau.int = NULL, tau.int.a = 0.1, tau.int.b = 0.1,
                   tau.time = NULL, tau.time.a = 0.1, tau.time.b = 0.1,
                   bw.init = NULL, canbw.sd = 0.05, bw.attempts = 50,
                   # starting value for PS
                   A = NULL,
                   # how many sites, days, and knots to keep. default is all
                   keep.sites = NULL, keep.days = NULL, keep.knots = NULL,
                   iters = 5000, burn = 1000, update = 10, keep.burn = FALSE,
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

  time <- rep(0, nt)
  for (t in 1:nt) {
    time[t] <- (t - nt / 2) / nt
  }

  miss <- is.na(y)
  if (length(thresh) == 1) {
    thresh <- matrix(thresh, ns, nt)
  }
  y <- ifelse(y < thresh, thresh, y)

  # initial missing values to be set at threshold
  # looping over time because threshold is a vector of ns length
  if (any(miss)) {
    missing.times <- which(colSums(miss) > 0)
  }

  # INITIAL VALUES:
  if (is.null(beta.time)) {
    beta.time <- matrix(0, ns, 2)
  }
  if (is.null(beta.int)) {
    beta.int <- matrix(0, ns, 2)
    beta.int[, 2] <- log(sqrt(6) * sd(as.vector(y), na.rm = TRUE) / pi)
    beta.int[, 1] <- -0.5772 * beta.int[, 2]
  }

  beta.int.mn <- beta.time.mn <- rep(0, 2)
  beta.pri.sd <- c(mu.beta.pri.sd, ls.beta.pri.sd)

  mu <- ls <- matrix(0, ns, nt)
  for (t in 1:nt) {
    mu[, t] <- beta.int[, 1] + beta.time[, 1] * time[t]
    ls[, t] <- beta.int[, 2] + beta.time[, 2] * time[t]
  }

  # get the initial covariance matrix for the gaussian process
  d <- rdist(s)
  diag(d) <- 0
  bw.min <- 1e-4
  bw.max <- max(d[upper.tri(d)])
  if (!is.null(bw.init)) {
    bw <- bw.init
  } else {
    bw <- (median(d[upper.tri(d)]) - bw.min) / 2
  }
  Sigma <- exp(-d / bw)
  Qb    <- chol2inv(chol(Sigma))
  logdetQb <- logdet(Qb)

  SS.int  <- SS.time  <- rep(0, 2)
  tau.int <- tau.time <- rep(1, 2)

  # adjust xi until it is valid for all y, mu, and sig
  if (xi < 0) {
    while (any(xi * (y - mu) / exp(ls) < -1, na.rm = TRUE)) {
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
  curll <- loglike(y = y, mu = mu, ls = ls, xi = xi, theta = theta,
                   theta.xi = theta.xi, thresh = thresh, alpha = alpha)
  # print(curll)

  # STORAGE:
  these.sites <- sort(sample(1:ns, keep.sites))  # which sites to keep GP
  these.days  <- sort(sample(1:nt, keep.days))   # which days to keep GP and A
  these.knots <- sort(sample(1:L, keep.knots))   # which knots to keep A
  keep.beta.int     <- array(0, dim = c(iters, ns, 2))
  keep.beta.time    <- array(0, dim = c(iters, ns, 2))
  keep.beta.int.mn  <- matrix(0, iters, 2)
  keep.beta.time.mn <- matrix(0, iters, 2)
  keep.xi           <- rep(0, iters)
  keep.tau.int      <- matrix(0, iters, 2)
  keep.tau.time     <- matrix(0, iters, 2)
  keep.bw           <- rep(0, iters)  # bandwidth for GP
  keep.A            <- array(0, dim = c(iters, keep.knots, keep.days))
  if (any(miss)) {
    # only record the missing data after burnin finishes
    keep.y <- matrix(0, iters - burn, sum(miss))
  } else {
    keep.y <- NULL
  }

  # TUNING:

  cuts <- exp(c(-1, 0, 1, 2, 5, 10))
  MH.a  <- rep(1, 100)
  att.a <- acc.a <- 0 * MH.a
  att.beta.int <- acc.beta.int <- MH.beta.int <- matrix(canbeta.int.sd, ns, 2)
  att.beta.time <- acc.beta.time <- MH.beta.time <- matrix(canbeta.time.sd, ns, 2)
  att.xi <- acc.xi <- MH.xi <- canxi.sd
  att.bw <- acc.bw <- MH.bw <- canbw.sd

  tic <- proc.time()[3]
  for (iter in 1:iters) {
    ####################################################
    ##############      Random effects A    ############
    ####################################################
    oldA  <- A
    this.update <- updateA(A = A, cuts = cuts, bins = bins, Ba = Ba,
                           theta = theta, theta.xi = theta.xi,
                           y = y, mu = mu, ls = ls, xi = xi, thresh = thresh,
                           alpha = alpha, curll = curll, MH = MH.a)

    A        <- this.update$A
    l1       <- this.update$l1
    theta    <- this.update$theta
    theta.xi <- this.update$theta.xi
    curll    <- this.update$curll

    ####################################################
    ##############      GEV parameters      ############
    ####################################################
    #### betas for mu ####
    this.update <- updateBeta1.int(beta.int = beta.int[, 1],
                                   beta.mn = beta.int.mn[1], SS = SS.int[1],
                                   tau = tau.int[1], beta.time = beta.time[, 1],
                                   time = time, y = y, theta = theta,
                                   theta.xi = theta.xi, mu = mu, ls = ls,
                                   xi = xi, thresh = thresh, alpha = alpha,
                                   Qb = Qb, curll = curll,
                                   acc = acc.beta.int[, 1],
                                   att = att.beta.int[, 1],
                                   MH = MH.beta.int[, 1])
    beta.int[, 1]     <- this.update$beta.int
    SS.int[1]         <- this.update$SS
    mu                <- this.update$mu
    curll             <- this.update$curll
    acc.beta.int[, 1] <- this.update$acc
    att.beta.int[, 1] <- this.update$att

    this.update <- updateBeta1.time(beta.time = beta.time[, 1],
                                    beta.mn = beta.time.mn[1],
                                    SS = SS.time[1], tau = tau.time[1],
                                    beta.int = beta.int[, 1], time = time,
                                    y = y, theta = theta,
                                    theta.xi = theta.xi, mu = mu, ls = ls,
                                    xi = xi, thresh = thresh, alpha = alpha,
                                    Qb = Qb, curll = curll,
                                    acc = acc.beta.time[, 1],
                                    att = att.beta.time[, 1],
                                    MH = MH.beta.time[, 1])
    beta.time[, 1]     <- this.update$beta.time
    SS.time[1]         <- this.update$SS
    mu                 <- this.update$mu
    curll              <- this.update$curll
    acc.beta.time[, 1] <- this.update$acc
    att.beta.time[, 1] <- this.update$att

    #### betas for ls ####
    this.update <- updateBeta2.int(beta.int = beta.int[, 2],
                                   beta.mn = beta.int.mn[2],
                                   SS = SS.int[2], tau = tau.int[2],
                                   beta.time = beta.time[, 2], time = time,
                                   y = y, theta = theta,
                                   theta.xi = theta.xi, mu = mu, ls = ls,
                                   xi = xi, thresh = thresh, alpha = alpha,
                                   Qb = Qb, curll = curll,
                                   acc = acc.beta.int[, 2],
                                   att = att.beta.int[, 2],
                                   MH = MH.beta.int[, 2])
    beta.int[, 2]     <- this.update$beta.int
    SS.int[2]         <- this.update$SS
    ls                <- this.update$ls
    curll             <- this.update$curll
    acc.beta.int[, 2] <- this.update$acc
    att.beta.int[, 2] <- this.update$att

    this.update <- updateBeta2.time(beta.time = beta.time[, 2],
                                    beta.mn = beta.time.mn[2],
                                    SS = SS.time[2], tau = tau.time[2],
                                    beta.int = beta.int[, 2], time = time,
                                    y = y, theta = theta,
                                    theta.xi = theta.xi, mu = mu, ls = ls,
                                    xi = xi, thresh = thresh, alpha = alpha,
                                    Qb = Qb, curll = curll,
                                    acc = acc.beta.time[, 2],
                                    att = att.beta.time[, 2],
                                    MH = MH.beta.time[, 2])
    beta.time[, 2]     <- this.update$beta.time
    SS.time[2]         <- this.update$SS
    ls                 <- this.update$ls
    curll              <- this.update$curll
    acc.beta.time[, 2] <- this.update$acc
    att.beta.time[, 2] <- this.update$att

    #### xi ####
    this.update <- updateXi(xi = xi, xi.min = xi.min, xi.max = xi.max,
                            xi.mn = xi.mn, xi.sd = xi.sd,
                            y = y, mu = mu, ls = ls, curll = curll,
                            theta = theta, theta.xi = theta.xi, thresh = thresh,
                            alpha = alpha, acc = acc.xi, att = att.xi,
                            MH = MH.xi)
    xi       <- this.update$xi
    theta.xi <- this.update$theta.xi
    curll    <- this.update$curll
    acc.xi   <- this.update$acc
    att.xi   <- this.update$att

    ####################################################
    ###########      GEV hyper parameters      #########
    ####################################################
    #### means ####
    for (p in 1:2) {
      # mu ~ N(0, beta.pri.sd)
      this.update <- updateGPMean(beta.sd = beta.pri.sd[p], Qb = Qb,
                                  beta.int = beta.int[, p],
                                  tau.int = tau.int[p],
                                  beta.time = beta.time[, p],
                                  tau.time = tau.time[p])
      beta.int.mn[p]  <- this.update$beta.int.mn
      SS.int[p]       <- this.update$SS.int
      beta.time.mn[p] <- this.update$beta.time.mn
      SS.time[p]      <- this.update$SS.time
    }

    # variances
    for (p in 1:2) {
      this.update <- updateGPTau(SS.int = SS.int[p],
                                 SS.time = SS.time[p],
                                 tau.a = 0.1, tau.b = 0.1, ns = ns)
      tau.int[p]  <- this.update$tau.int
      tau.time[p] <- this.update$tau.time
    }

    #### bandwidth for GP ####
    this.update <- updateBW(bw = bw, bw.min = bw.min, bw.max = bw.max,
                            Qb = Qb, logdetQb = logdetQb, d = d,
                            beta.int = beta.int, tau.int = tau.int,
                            SS.int = SS.int, beta.int.mn = beta.int.mn,
                            beta.time = beta.time, tau.time = tau.time,
                            SS.time = SS.time, beta.time.mn = beta.time.mn,
                            acc = acc.bw, att = att.bw, MH = MH.bw)
    bw       <- this.update$bw
    Qb       <- this.update$Qb
    logdetQb <- this.update$logdetQb
    SS.int   <- this.update$SS.int
    SS.time  <- this.update$SS.time
    acc.bw   <- this.update$acc
    att.bw   <- this.update$att

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

      this.update <- mhUpdate(acc = acc.beta.int[, p], att = att.beta.int[, p],
                              MH = MH.beta.int[, 1], target.min = 0.3,
                              target.max = 0.6, nattempts = beta.int.attempts)
      acc.beta.int[, p] <- this.update$acc
      att.beta.int[, p] <- this.update$att
      MH.beta.int[, p]  <- this.update$MH

      this.update <- mhUpdate(acc = acc.beta.time[, p],
                              att = att.beta.time[, p],
                              MH = MH.beta.time[, p], target.min = 0.3,
                              target.max = 0.6, nattempts = beta.time.attempts)
      acc.beta.time[, p] <- this.update$acc
      att.beta.time[, p] <- this.update$att
      MH.beta.time[, p]  <- this.update$MH

      this.update <- mhUpdate(acc = acc.xi, att = att.xi, MH = MH.xi,
                              target.min = 0.3, target.max = 0.6,
                              nattempts = xi.attempts)
      acc.xi <- this.update$acc
      att.xi <- this.update$att
      MH.xi  <- this.update$MH

      this.update <- mhUpdate(acc = acc.bw, att = att.bw, MH = MH.bw,
                              target.min = 0.3, target.max = 0.6,
                              nattempts = bw.attempts)
      acc.bw <- this.update$acc
      att.bw <- this.update$att
      MH.bw  <- this.update$MH
    }

    ####################################################
    ##############      Impute missing      ############
    ####################################################
    if (iter > burn) {
      if (any(miss)) {
        y.tmp <- y
        for (t in missing.times) {
          # calculate mu and sigma
          miss.t     <- miss[, t]  # logical vector of whether it's missing
          theta.xi.t <- theta.xi[miss.t, t]
          mu.t       <- mu[miss.t, t]
          sig.t      <- exp(ls[miss.t, t])
          mu.star.t  <- mu.t + sig.t * (theta.xi.t - 1) / xi
          sig.star.t <- alpha * sig.t * theta.xi.t
          xi.star    <- alpha * xi
          # get unit frechet
          these.miss <- -1 / log(runif(sum(miss.t)))
          # transform to correct marginals
          y.tmp[miss.t, t] <- mu.star.t + sig.star.t *
            (these.miss^xi.star - 1) / xi.star
        }
        keep.y[(iter - burn), ] <- y.tmp[miss]
      }

    }

    #KEEP TRACK OF STUFF:
    keep.beta.int[iter, , ]   <- beta.int
    keep.beta.time[iter, , ]  <- beta.time
    keep.beta.int.mn[iter, ]  <- beta.int.mn
    keep.beta.time.mn[iter, ] <- beta.time.mn
    keep.xi[iter]             <- xi
    keep.tau.int[iter, ]      <- tau.int
    keep.tau.time[iter, ]     <- tau.time
    keep.bw[iter]             <- bw
    keep.A[iter, , ]          <- A[these.knots, these.days]

    #DISPLAY CURRENT VALUE:

    if (iter %% update == 0) {
      cat("    Finished fit:", iter, "of", iters, "iters \n")
      params <- c("mu", "ls")
      if (iterplot) {
        acc.rate.beta.int  <- round(acc.beta.int / att.beta.int, 3)
        acc.rate.beta.time <- round(acc.beta.time / att.beta.time, 3)
        acc.rate.xi        <- round(acc.xi / att.xi, 3)
        acc.rate.bw        <- round(acc.bw / att.bw, 3)

        if (iter > burn) {
          start <- burn + 1
        } else {
          start <- max(iter - 3000, 1)
        }

        par(mfrow = c(5, 4))

        for (i in 1:2) {
          for (p in 1:2) {
            plot(keep.beta.int[start:iter, i, p], type = "l",
                 main = paste(params[p], " beta int ", i, sep = ""),
                 ylab = paste("MH: ", round(MH.beta.int[i, p], 3)),
                 xlab = acc.rate.beta.int[i, p])
            plot(keep.beta.time[start:iter, i, p], type = "l",
                 main = paste(params[p], " beta time ", i, sep = ""),
                 ylab = paste("MH: ", round(MH.beta.time[i, p], 3)),
                 xlab = acc.rate.beta.time[i, p])
          }
        }

        for (p in 1:2) {
          plot(keep.tau.int[start:iter, p], type = "l",
               main = paste(params[p], " tau int", sep = ""))
          plot(keep.tau.time[start:iter, p], type = "l",
               main = paste(params[p], " tau time", sep = ""))
        }
        for (p in 1:2) {
          plot(keep.beta.int.mn[start:iter, p], type = "l",
               main = paste(params[p], " beta int mean", sep = ""))
          plot(keep.beta.time.mn[start:iter, p], type = "l",
               main = paste(params[p], " beta time mean", sep = ""))
        }

        plot(keep.xi[start:iter], main = bquote(xi),
             xlab = acc.rate.xi, type = "l",
             ylab = paste("MH =", round(MH.xi, 3)))

        plot(keep.bw[start:iter], main = "gaussian process bandwidth",
             xlab = acc.rate.bw, type = "l",
             ylab = paste("MH =", round(MH.bw, 3)))

        plot(log(keep.A[start:iter, 1, 1]), main = "log(A[1, 1])", type = "l")
        plot(log(keep.A[start:iter, 2, 1]), main = "log(A[2, 1])", type = "l")
      }
    }

  } #end iter

  toc <- proc.time()[3]

  if (keep.burn) {
    return.iters <- 1:iters
  } else {
    return.iters <- (burn + 1):iters
  }

  list(beta.int     = keep.beta.int[return.iters, , , drop = FALSE],
       beta.time    = keep.beta.time[return.iters, , , drop = FALSE],
       beta.int.mn  = keep.beta.int.mn[return.iters, , drop = FALSE],
       beta.time.mn = keep.beta.time.mn[return.iters, , drop = FALSE],
       xi           = keep.xi[return.iters],
       tau.int      = keep.tau.int[return.iters, , drop = FALSE],
       tau.time     = keep.tau.time[return.iters, , drop = FALSE],
       bw           = keep.bw[return.iters],
       A            = keep.A[return.iters, , , drop = FALSE],
       y.pred       = keep.y,  # only stores for post-burnin
       timing       = toc - tic)
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
  niters    <- length(start:end)
  beta.int  <- matrix(mcmcoutput$beta.int[start:end, , drop = F], niters, p)
  beta.time <- matrix(mcmcoutput$beta.time[start:end, , drop = F], niters, p)
  xi        <- mcmcoutput$xi[start:end]
  A         <- mcmcoutput$A[start:end, , , drop = F]

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