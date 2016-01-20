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
#   beta1   := initial values for the GEV location coefficients
#   beta2   := initial values for the GEV scale coefficients
#   xi      := initial value for the GEV shape parameter
#   beta.sd := prior sd of the elements of beta1 and beta2
#
# OUTPUTS:
#
#   beta1    := Posterior samples of beta1
#   beta2    := Posterior samples of beta2
#   xi       := Posterior samples of xi
#   logA     := Posterior samples of log(A)
#   theta.mn := Posterior mean of theta
#
################################################################################

ReShMCMC<-function(y, X, thresh, B, alpha,
                   beta1 = NULL, beta2 = NULL, xi = 0.001, beta.sd = 10,
                   keep.burn = FALSE, iters = 5000, burn = 1000, update = 10,
                   iterplot = FALSE){
  # BOOKKEEPING
  
  ns   <- nrow(y)
  nt   <- ncol(y)
  L    <- ncol(B)
  p    <- dim(X)[3]
  
  y    <- ifelse(y < thresh, thresh, y)
  
  # dPS approximation:
  
  npts      <- 50
  Ubeta     <- qbeta(seq(0, 1, length = npts + 1), 0.5, 0.5)
  MidPoints <- (Ubeta[-1] + Ubeta[-(npts + 1)]) / 2
  BinWidth  <- Ubeta[-1] - Ubeta[-(npts + 1)]
  bins      <- list(npts = npts, MidPoints = MidPoints, BinWidth = BinWidth)
  
  # INITIAL VALUES:
  
  beta1    <- rep(0, p)
  beta2    <- rep(0, p)
  beta2[1] <- log(sqrt(6) * sd(as.vector(y)) / pi)
  beta1[1] <- mean(y) - exp(beta2[1]) * 0.577      
  
  
  mu <- logsig <- 0
  for(j in 1:p){
    mu     <- mu     + X[, , j] * beta1[j]
    logsig <- logsig + X[, , j] * beta2[j]
  }
  A     <- matrix(1, L, nt)
  
  Ba    <- B^(1 / alpha)
  theta <- (Ba %*% A)^alpha
  curll <- loglike(y, theta, mu, logsig, xi, thresh, alpha)
  
  # STORAGE:
  
  keep.beta1 <- matrix(0, iters, p)
  keep.beta2 <- matrix(0, iters, p)
  keep.xi    <- rep(0, iters)
  keep.A     <- array(0, dim = c(iters, L, nt))
  theta.mn   <- 0 
  
  # TUNING:
  
  cuts <- exp(c(-1, 0, 1, 2, 5, 10))
  MHa  <- rep(1,100)
  atta <- acca <- 0*MHa
  attb <- accb <- MHb <- matrix(.1,p,3)
  
  tic <- proc.time()[3]
  for (iter in 1:iters) {
    
    ####################################################
    ##############      Random effects A    ############
    ####################################################
    oldA  <- A
    l1    <- get.level(A, cuts)
    CANA  <- A * exp(MHa[l1] * rnorm(nt * L))
    l2    <- get.level(CANA, cuts)
    q     <- dPS(CANA, alpha, bins) - 
             dPS(A, alpha, bins) +
             dlognormal(A, CANA, matrix(MHa[l2], L, nt)) -
             dlognormal(CANA, A, matrix(MHa[l1], L, nt))        
    
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
    
    for (j in 1:p) {
      
      # beta1
      attb[j, 1] <- attb[j, 1] + 1
      canb       <- rnorm(1, beta1[j], MHb[j, 1])
      canmu      <- mu + X[, , j] * (canb - beta1[j])
      canll      <- loglike(y, theta, canmu, logsig, xi, thresh, alpha)
      R          <- sum(canll - curll) + 
                    dnorm(canb, 0, beta.sd, log = TRUE) -
                    dnorm(beta1[j], 0, beta.sd, log = TRUE)
      if (log(runif(1)) < R) {
        accb[j, 1] <- accb[j, 1] + 1
        beta1[j]   <- canb
        mu         <- canmu
        curll      <- canll
      }
      
      # beta2
      attb[j, 2] <- attb[j, 2] + 1
      canb       <- rnorm(1, beta2[j], MHb[j, 2])
      canlogs    <- logsig + X[, , j] * (canb - beta2[j])
      canll      <- loglike(y, theta, mu, canlogs, xi, thresh, alpha)
      R          <- sum(canll - curll) + 
                    dnorm(canb, 0, beta.sd, log = TRUE) -
                    dnorm(beta2[j], 0, beta.sd, log = TRUE)
      if (log(runif(1)) < R) {
        accb[j, 2] <- accb[j, 2] + 1
        beta2[j]   <- canb
        logsig     <- canlogs
        curll      <- canll
      }
    }
    
    # xi
    attb[1, 3] <- attb[1, 3] + 1
    canxi      <- rnorm(1, xi, MHb[1, 3])
    canll      <- loglike(y, theta, mu, logsig, canxi, thresh, alpha)
    R          <- sum(canll - curll) + 
                  dnorm(canxi, 0, 0.5, log = TRUE) -
                  dnorm(xi, 0, 0.5, log = TRUE)
    if (log(runif(1)) < R) {
      accb[1, 3] <- accb[1, 3] + 1
      xi         <- canxi
      curll      <- canll
    }
    
    # TUNING
    
    for (j in 1:length(MHa)) {
      acca[j] <- acca[j] + sum(oldA[l1 == j] != A[l1 == j])
      atta[j] <- atta[j] + sum(l1 == j)
      if (iter < burn / 2 & atta[j] > 100) {
        if (acca[j] / atta[j] < 0.3) { MHa[j] <- MHa[j] * 0.9 }
        if (acca[j] / atta[j] > 0.6) { MHa[j] <- MHa[j] * 1.1 }
        acca[j] <- atta[j] <- 0
      }
    }
    
    for (k in 1:p) { for (j in 1:3) { if (iter < burn / 2 & attb[k, j] > 50) {
      if (accb[k, j] / attb[k, j] < 0.3) { MHb[k, j] <- MHb[k, j] * 0.9}
      if (accb[k, j] / attb[k, j] > 0.6) { MHb[k, j] <- MHb[k, j] * 1.1}
      accb[k, j] <- attb[k, j] <- 0
    }}}
    
    
    #KEEP TRACK OF STUFF:
    
    keep.beta1[iter, ] <- beta1
    keep.beta2[iter, ] <- beta2
    keep.xi[iter]      <- xi
    keep.A[iter, , ]   <- A
    if (iter > burn) {
      theta.mn <- theta.mn + theta / (iters - burn)
    }
    
    
    #DISPLAY CURRENT VALUE:
    
    if (iter %% update == 0) {
      if (iterplot) {
        par(mfrow = c(3, 2))
        plot(keep.beta1[1:iter, 1], main = "beta1[1]", type = "l")
        plot(keep.beta1[1:iter, p], main = "beta1[p]", type = "l")
        plot(keep.beta2[1:iter, 1], main = "beta2[1]", type = "l")
        plot(keep.beta2[1:iter, p], main = "beta2[p]", type = "l")
        plot(log(keep.A[1:iter, 1, 1]), main = "log(A[1, 1])", type = "l")
        plot(log(keep.A[1:iter, L, 1]), main = "log(A[L, 1])", type = "l")
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
  
  list(beta1 = keep.beta1[return.iters, , drop = FALSE],
       beta2 = keep.beta2[return.iters, , drop = FALSE],
       xi = keep.xi[return.iters],
       theta.mn = theta.mn,
       A = keep.A[return.iters, , , drop = FALSE], 
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
  mu_star  <- mu + sigma*((theta^xi) - 1) / xi
  sig_star <- alpha * sigma * (theta^xi)
  xi_star  <- alpha * xi
  
  tx       <- (1 + xi_star * (y - mu_star) / sig_star)^(-1 / xi_star)
  ll       <- -tx + (y > thresh) * ((xi_star + 1) * log(tx) - log(sig_star))
  ll       <- ifelse(is.na(y), 0, ll)
  ll       <- ifelse(is.na(ll), -Inf, ll)
  
  return(ll)
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