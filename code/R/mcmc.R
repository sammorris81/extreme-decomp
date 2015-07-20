######################################################
########      THE MAIN MCMC FUNCTION       ###########
######################################################

Bayes_GEV <- function(y, S, knots = NULL,
                      X = NULL,
                      Sp = NULL, Xp = NULL,
                      init.beta = c(0, 0, 0.0001),
                      init.alpha = 0.25,
                      init.range = 1,
                      init.bw = 1,
                      vary = c(T, F, F),
                      pri.mn.range = -2,
                      pri.sd.range = 1,
                      pri.mn.bw = 0,
                      pri.sd.bw = 1,
                      pri.var.a = 0.1,
                      pri.var.b = 0.1,
                      pri.alpha.a = 1,
                      pri.alpha.b = 1,
                      pri.sd.beta = 10,
                      pri.mn.gev = c(0, 0, 0),
                      pri.sd.gev = c(10, 1, 0.25),
                      keep.samples = T,
                      iters = 200, burn = 50, update = 50, nthin = 1){
  
  start.time <- proc.time()
  library(fields)
  library(emulator)
  
  if (is.null(knots)) {knots <- S}
  if (is.null(X)) {X <- cbind(1, S)}
  
  #BOOKEEPING
  n <- nrow(y)
  nt <- ncol(y)
  years <- 3 * ((1:nt) - ((nt + 1) / 2)) / ((nt + 1) / 2)  # put years in (-3, 3)
  nF <- nrow(knots)
  p <- ncol(X)
  
  d <- rdist(S,S)
  diag(d) <- 0
  dw2 <- rdist(S, knots)^2
  dw2[dw2 < 0.0001] <- 0
  
  #INITIAL VALUES:
  beta <- matrix(init.beta, n, 3, byrow = T)
  if (vary[1]) {
    for (j in 1:n) {
      beta[j, 1] <- get.inits.mu(y[j, ], beta[j, 2], beta[j, 3])
    }
  }
  mnb <- matrix(0, p, 3)
  mnb[1,] <- colMeans(beta)
  taub <- rep(1, 3)
  lograngeb <- log(init.range)
  logrange <- log(init.bw)
  logs <- matrix(2, nF, nt)
  u <- matrix(0.5, nF, nt)
  alpha <- init.alpha
  
  #COMPUTE QUANTITIES USED IN THE POSTERIOR
  Qb <- solve(exp(-d / exp(lograngeb)))
  fac <- make.kern(dw2, logrange)
  FAC <- stdKern(fac)
  cur.ll <- matrix(0, n, nt)
  theta <- make.theta(FAC, logs, alpha)
  for (j in 1:n) {
    cur.ll[j, ] <- loglike(y[j, ], beta[j, 1], beta[j, 2], beta[j, 3], 
                          theta[j, ], alpha)
  }
  
  #CREATE PLACES TO STORE THE OUTPUT
  keepers <- matrix(0, iters, 7)
  colnames(keepers) <- c("alpha", "spatial range of the GEV params",
                         "kernel bandwidth",
                         "GEV loc site 1", "GEV log scale site 1",
                         "GEV shape site 1", "log likelihood")
  GEVs <- c("GEV location", "GEV log scale", "GEV shape")
  SITEs <- paste("Site", 1:n)
  beta.mn <- beta.var <- 0 * beta
  colnames(beta.mn) <- colnames(beta.var) <- GEVs
  rownames(beta.mn) <- rownames(beta.var) <- SITEs
  locs <- NULL
  if (keep.samples) {
    locs <- array(0, c(iters, n, 3))
    dimnames(locs) <- list(paste("sample", 1:iters), SITEs, GEVs)
  }
  Yp <- Yp1 <- Yp2 <- locsp <- beta.mn.p <- beta.var.p <- NULL
  if (!is.null(Sp)) {
    np <- nrow(Sp)
    SITEs <- paste("Pred site", 1:np)
    TIMEs <- paste("Replication", 1:nt)
    Yp1 <- Yp2 <- matrix(0, np, nt)
    colnames(Yp1) <- colnames(Yp2) <- TIMEs
    rownames(Yp1) <- rownames(Yp2) <- SITEs
    beta.mn.p <- beta.var.p <- matrix(0, np, 3)   
    colnames(beta.mn.p) <- colnames(beta.var.p) <- GEVs
    rownames(beta.mn.p) <- rownames(beta.var.p) <- SITEs
    if (keep.samples) {
      locsp <- array(0, c(iters, np, 3))
      dimnames(locsp) <- list(paste("sample", 1:iters), SITEs, GEVs)
      Yp <- array(0, c(iters, np, nt))
      dimnames(Yp) <- list(paste("sample", 1:iters), SITEs, TIMEs)
    }
    dw2p <- rdist(Sp, knots)^2
    dw2p[dw2p < 0.0001] <- 0
    if (sum(vary) > 0) {
      d12 <- rdist(Sp, S)
      d22 <- rdist(Sp, Sp)
      d12[d12 < 0.0001] <- 0
      diag(d22) <- 0
    }
  }
  
  
  #SETTINGS FOR THE M-H CANDIDATE DISTRIBUTION
  att.a <- acc.a <- mh.a <- rep(1, nF)  # spatial for location
  att.b <- acc.b <- mh.b <- rep(1, nF)  # time for location
  att.c <- acc.c <- mh.c <- rep(1, nF)  # spatial for scale
  att.d <- acc.d <- mh.d <- 0.3         # shape
  
  att.alpha <- acc.alpha <- mh.alpha <- 0.01
  att.u     <- acc.u <- mh.u <- matrix(0.1, nF, nt)  # maybe not needed
  att.logs  <- acc.logs <- mh.logs <- c(3, 2, 2, 2, rep(1, 10))
  cuts <- seq(0, 15, 2)
  
  
  #START SAMPLING!
  for (i in 1:iters) { for (rep in 1:nthin) {
    
    ##########################################################
    ##############      Random effects S and U    ############
    ##########################################################
    level <- olds <- logs
    
    logs.update <- updateLogs(logs = logs, y = y, mu = mu, sigma = sigma, 
                              xi = xi, theta = theta, alpha = alpha, u = u, 
                              cur.ll = cur.ll, FAC = FAC, cuts = cuts, 
                              mh = mh.logs)
    logs   <- a.update$logs
    theta  <- a.update$theta
    cur.ll <- a.update$cur.ll
    
    for (j in 1:length(MHs)) {
      acc.logs[j] <- acc.logs[j] + sum(olds[level == j] != logs[level == j])
      att.logs[j] <- att.logs[j] + sum(level == j)
    }
    for (j in 1:length(att.logs)) { if ((i < burn / 2) & (att.logs[j] > 50)) {
      if (acc.logs[j] / att.logs[j] < 0.3) { mh.logs[j] <- mh.logs[j] * 0.9 }
      if (acc.logs[j] / att.logs[j] > 0.6) { mh.logs[j] <- mh.logs[j] * 1.1 }
      acc.logs[j] <- att.logs[j] <- 0
    }}
    
    canu <- rtnorm(u)
    R <- h1(logs, canu, alpha, log = TRUE) -
         h1(logs, u, alpha, log = TRUE) +
         dtnorm(u, canu) -
         dtnorm(canu, u)
    acc <- matrix(runif(nt * nF), nF, nt)
    u <- ifelse(acc < exp(R), canu, u)
    
    ##########################################################
    ##############              alpha             ############
    ##########################################################
    alpha.update <- updateAlpha(y = y, alpha = alpha, logs = logs, FAC = FAC, 
                                theta = theta, mu = mu, sigma = sigma, xi = xi, 
                                u = u, cur.ll = cur.ll, acc = acc.beta, 
                                att = att.beta, mh = mh.beta)
    alpha    <- alpha.update$alpha
    theta    <- alpha.update$theta
    cur.ll   <- alpha.update$cur.ll
    acc.beta <- alpha.update$acc
    att.beta <- alpha.update$att
    
    if (iter < burn / 2) {
      mh.update <- mhUpdate(acc = acc.alpha, att = att.alpha, mh = mh.alpha,
                            nattempts = alpha.attempts)
      acc.alpha <- mh.update$acc
      att.alpha <- mh.update$att
      mh.alpha  <- mh.update$mh
    }
    
    ##########################################################
    ##############       KERNEL BANDWIDTH         ############
    ##########################################################
    
#     #### Not needed
#     attb[5] <- attb[5] + 1
#     canlogrange <- rnorm(1, logrange, MHb[5])
#     canfac <- make.kern(dw2, canlogrange)
#     canFAC <- stdKern(canfac)
#     canll <- curll
#     cantheta <- make.theta(canFAC, logs, alpha)
#     for (t in 1:nt) {
#       canll[, t] <- loglike(y[, t], beta[, 1], beta[, 2], beta[, 3],
#                             cantheta[, t], alpha)
#     }
#     R <- sum(canll - curll) +
#          dnorm(canlogrange, pri.mn.bw, pri.sd.bw, log = TRUE) -
#          dnorm(logrange, pri.mn.bw, pri.sd.bw, log = TRUE)
#     if (!is.na(exp(R))) { if (runif(1) < exp(R)) {
#       logrange <- canlogrange
#       fac <- canfac
#       FAC <- canFAC
#       theta <- cantheta
#       curll <- canll
#       accb[5] <- accb[5] + 1
#     }}
    
    ##########################################################
    ##############          GEV PARAMETERS        ############
    ##########################################################
    
    ## SPATIALLY VARYING PARAMETERS
    for (l in 1:3) { if (i > 50 & vary[l]) {
      Xb <- X %*% mnb[, l]
      
      for (j in 1:n) {
        VVV <- taub[l] * Qb[j, j]
        MMM <- taub[l] * Qb[j, j] * Xb[j] -
               taub[l] * sum(Qb[-j, j] * (beta[-j, l] - Xb[-j]))
        
        attb[l] <- attb[l] + 1
        canb <- beta[j, ]
        canb[l] <- rnorm(1, beta[j, l], MHb[l])
        canll <- loglike(y[j, ], canb[1], canb[2], canb[3], theta[j, ], alpha)
        
        R <- sum(canll - curll[j, ]) +
             dnorm(canb[l], MMM / VVV, 1 / sqrt(VVV), log = TRUE) -
             dnorm(beta[j, l], MMM / VVV, 1 / sqrt(VVV), log = TRUE)
        
        if (!is.na(exp(R))) { if (runif(1) < exp(R)) {
          beta[j, ] <- canb
          curll[j, ] <- canll
          accb[l] <- accb[l] + 1
        }}
      }
      
    }}
    
    ## SPATIALLY CONSTANT PARAMETERS
    for (l in 1:3) { if (i > 50 & !vary[l]) {
      attb[l] <- attb[l] + 1
      canb <- beta
      canb[, l] <- rnorm(1, beta[1, l], MHb[l])         
      canb[, l] <- beta[1, l] + MHb[l] * rt(1, df = 5)
      canll <- curll
      for (t in 1:nt) {
        canll[, t] <- loglike(y[, t], canb[, 1], canb[, 2], canb[, 3], 
                              theta[, t], alpha)
      }
      R <- sum(canll - curll) +
           dnorm(canb[1, l], pri.mn.gev[l], pri.sd.gev[l], log = TRUE) -
           dnorm(beta[1, l], pri.mn.gev[l], pri.sd.gev[l], log = TRUE)
      if (!is.na(exp(R))) { if (runif(1) < exp(R)) {
        beta <- canb
        curll <- canll
        accb[l] <- accb[l] + 1
      }}
    }}
    
    ##########################################################
    ##############    Spatial hyperparameters     ############
    ##########################################################
    
    #### not needed?
    if (sum(vary) > 0) {
      tXQ <- t(X) %*% Qb
      tXQX <- tXQ %*% X
    }
    
    for (l in 1:3) { if (vary[l]) {
      #MEAN 
      VVV <- solve(taub[l] * tXQX + (1 / pri.sd.beta^2) * diag(p))
      MMM <- taub[l] * tXQ %*% beta[, l]
      mnb[, l] <- VVV %*% MMM + t(chol(VVV)) %*% rnorm(p)
      
      #VARIANCE
      SS <- quad.form(Qb, beta[, l] - X %*% mnb[, l])
      taub[l] <- rgamma(1, n / 2 + pri.var.a, SS / 2 + pri.var.b)
    }}
    
    #SPATIAL RANGE
    if (sum(vary) > 0) {
      attb[6] <- attb[6] + 1
      canlograngeb <- rnorm(1, lograngeb, MHb[6])
      canQb <- solve(exp(-d / exp(canlograngeb)))
      R<- 0.5 * sum(vary) * logdet(canQb) -
          0.5 * sum(vary) * logdet(Qb) +
          dnorm(canlograngeb, pri.mn.range, pri.sd.range, log = T) -
          dnorm(lograngeb, pri.mn.range, pri.sd.range, log = T)
      for (l in 1:3) { if (vary[l]) {
        R <- R - 0.5 * taub[l] * quad.form(canQb, beta[, l] - X %*% mnb[, l])
        R <- R + 0.5 * taub[l] * quad.form(Qb, beta[, l] - X %*% mnb[, l])
      }}
      if (!is.na(exp(R))) { if (runif(1) < exp(R)) {
        lograngeb <- canlograngeb
        Qb <- canQb
        accb[6] <- accb[6] + 1
      }}
    }
    
    #########  TUNE THE CANDIDATE DISTRIBUTION  #######
    for (j in 1:length(accb)) { if (i < burn / 2 & attb[j] > 50) {
      if (accb[j] / attb[j] < 0.3) {MHb[j] <- MHb[j] * 0.9}
      if (accb[j] / attb[j] > 0.6) {MHb[j] <- MHb[j] * 1.1}
      accb[j] <- attb[j] <- 0
    }}
    
  }#end nthin
    
    
    
    #KEEP TRACK OF STUFF:
    keepers[i, ] <- c(alpha, exp(lograngeb),
                      exp(logrange), beta[1, ], sum(curll))
    if (i > burn) {
      nnn <- iters - burn
      beta.mn <- beta.mn + beta / nnn
      beta.var <- beta.var + beta * beta/nnn
    }
    if (keep.samples) { locs[i, , ] <- beta }
    
    
    #MAKE PREDICTIONS AT NEW LOCATIONS
    if (!is.null(Sp)) {
      YYY <- matrix(0, np, nt)
      facp <- make.kern(dw2p, logrange)
      FACp <- stdKern(facp)
      thetap <- make.theta(FACp, logs, alpha)^alpha
      
      bp <- matrix(beta[1, ], np, 3, byrow = TRUE)
      for (j in 1:3) { if (vary[j]) {
        RRR <- beta[, j] - X %*% mnb[, j]
        bp[, j] <- Xp %*% mnb[, j] +
          proj.beta(RRR, d12, d22, Qb, taub[j], lograngeb)
      }}
      for (t in 1:nt) {
        res <- rGEV(np, thetap[, t], alpha * thetap[, t], alpha)
        YYY[, t] <- bp[, 1] + exp(bp[, 2]) * (res^(bp[, 3]) - 1) / bp[, 3]
      }
      
      if (i > burn) {
        Yp1 <- Yp1 + YYY / (iters - burn)
        Yp2 <- Yp2 + YYY * YYY / (iters - burn)
        beta.mn.p <- beta.mn.p + bp / (iters - burn)
        beta.var.p <- beta.mn.p + bp * bp / (iters - burn)
      }
      if (keep.samples) {
        Yp[i, , ] <- YYY
        locsp[i, , ] <- bp
      }
    }
    
    #DISPLAY CURRENT VALUE:
    if (i %% update == 0) {
      par(mfrow = c(3, 2))
      plot(keepers[1:i, 1], ylab = "alpha", xlab = "iteration", type = "l")
      plot(keepers[1:i, 3], ylab = "bandwidth", xlab = "iteration", type = "l")
      plot(keepers[1:i, 4], ylab = "GEV loc site 1", xlab = "iteration",
           type = "l")
      plot(keepers[1:i, 5], ylab = "GEV log scale site 1", xlab = "iteration",
           type = "l")
      plot(keepers[1:i, 6], ylab = "GEV shape site 1", xlab = "iteration",
           type = "l")
      plot(keepers[1:i, 7], ylab = "Log likelihood", xlab = "iteration",
           type = "l")
    }
    
  }
  
  stop.time <- proc.time()
  
  #Return output:
  
  results <- list(time = stop.time - start.time,
                  beta.var = beta.var - beta.mn^2,
                  beta.mn = beta.mn,
                  beta.samples = locs,
                  beta.var.pred = beta.var.p - beta.mn.p^2,
                  beta.mn.pred = beta.mn.p,
                  beta.samples.pred = locsp,
                  Y.mn.pred = Yp1,
                  Y.var.pred = Yp2 - Yp1^2,
                  Y.samples.pred = Yp,
                  parameters = keepers)
  
  return(results)
}





