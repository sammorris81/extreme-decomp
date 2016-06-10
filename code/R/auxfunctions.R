# Rcpp functions
if (!exists("cppaux.load")) {
  sourceCpp(file = "auxfunctions.cpp")
  cppaux.load <- TRUE
}

#############################################################
########   FUNCTIONS TO COMPUTE INITIAL VALUES    ###########
#############################################################

get.inits.mu <- function(y, ls = 0, xi = 0.1){
  m <- median(y, na.rm = TRUE)
  mu <- m - exp(ls) * (log(2)^(-xi) - 1) / xi
  return(mu)
}


######################################################
########   THE POSITIVE STABLE DENSITY    ###########
######################################################

h1 <- function(logs, u, alpha, log = TRUE){
  s <- exp(logs)
  psi <- pi*u
  c <- (sin(alpha * psi) / sin(psi))^(1 / (1 - alpha))
  c <- c * sin((1 - alpha) * psi) / sin(alpha * psi)

  logd <- log(alpha) - log(1 - alpha) - (1 / (1 - alpha)) * logs +
    log(c) - c * (1 / s^(alpha / (1 - alpha))) +
    logs

  if (!log) {
    logd <- exp(logd)
  }

  return(logd)
}

######################################################
##########     Densities and posteriors    ###########
######################################################
loglike <- function(y, mu, ls, xi, theta, theta.xi, thresh, alpha){
  sigma    <- exp(ls)

  these <- (y > thresh) & !is.na(y) # likelihood changes if y > thresh
  if (abs(xi) <= 1e-4) {
    # theta.xi  <- theta^xi
    sig.star <- alpha * sigma  # theta.xi = 1 when xi = 0
    t.y <- (theta * exp(-(y - mu) / sigma))^(1 / alpha)
    ll <- -t.y
    ll[these] <- ll[these] - log(sig.star[these]) + log(t.y[these])
  } else {
    # theta.xi  <- theta^xi
    mu.star  <- mu + sigma * ((theta.xi) - 1) / xi
    sig.star <- alpha * sigma * (theta.xi)
    xi.star  <- alpha * xi
    t.y <- (1 + xi.star * (y - mu.star) / sig.star)^(-1 / xi.star)
    ll <- -t.y
    ll[these] <- ll[these] + (xi.star + 1) * log(t.y[these]) -
      log(sig.star[these])
  }
  ll[is.na(y)] <- 0
  ll[is.na(ll)] <- -Inf

  return(ll)
}

logd <- function(theta, v){
  sum(log(theta) - theta * v)
}


rgevspatial <- function(nreps, S, knots, mu = 1, sig = 1, xi = 1, alpha = 0.5,
                        bw = 1){
  library(evd)
  library(BMAmevt)

  n      <- nrow(S)
  nknots <- nrow(knots)

  d             <- rdist(S, knots)
  d[d < 0.0001] <- 0
  w             <- make.kern(d^2, log(bw))
  K             <- stdKern(w)^(1 / alpha)

  y <- matrix(0, n, nreps)
  for (t in 1:nreps) {
    A     <- rep(0, nknots)
    for (j in 1:nknots) {
      A[j] <- rstable.posit(alpha)
    }
    theta <- (K %*% A)^alpha

    xi_star  <- alpha*xi
    mu_star  <- mu+sig*(theta^xi-1)/xi
    sig_star <- alpha*sig*theta^xi

    y[,t]    <- rgev(n,mu_star,sig_star,xi_star)
  }

  return(y)
}

rX <- function(ns, nt, np) {
  X.sp <- matrix(rnorm(ns * (np - 2), 0, 0.3), ns, np - 2)
  X.t <- rnorm(nt, 0, 0.3)
  X <- array(1, dim = c(ns, nt, np))
  for (t in 1:nt) {
    X[, t, 2:np] <- cbind(rep(X.t[t], ns), X.sp)
  }
  return(X)
}

ll.ind <- function(beta, X, y, xi = -0.1) {
  # nt <- dim(X)[2]
  # np <- dim(X)[3]

  np <- dim(X)[2]
  beta1 <- beta[1:np]
  beta2 <- beta[(np + 1):(2 * np)]
  ll <- 0

  mu <- X %*% beta1
  sig <- exp(X %*% beta2)
  # if (any(sig[!is.na(y[, t])] == 0)) {
  #   return(Inf)
  # }
  ll <- sum(
    dgev(y, loc = mu, scale = sig, shape = xi, log = TRUE),
    na.rm = TRUE)

  return(-ll)
}

ll.ind.xi <- function(beta, X, y) {
  # nt <- dim(X)[2]
  # np <- dim(X)[3]

  np <- dim(X)[2]
  beta1 <- beta[1:np]
  beta2 <- beta[(np + 1):(2 * np)]
  xi <- tail(beta, 1)
  ll <- 0

  mu <- X %*% beta1
  sig <- exp(X %*% beta2)
  if (any(sig[!is.na(y)] == 0)) {
    return(Inf)
  }
  # print(sig)
  ll <- sum(
    dgev(y, loc = mu, scale = sig, shape = xi, log = TRUE),
    na.rm = TRUE)

  return(-ll)
}


#### Logposteriors and gradients - needs to be done for each timepoint
logpost.beta1.int <- function(beta.int, beta.int.mn, tau, Qb,
                              beta.time, time, y, ls, xi, theta,
                              theta.xi, thresh, alpha) {

  lp1 <- -0.5 * tau * quad.form(Qb, beta.int - beta.int.mn)

  mu  <- beta.int + beta.time * time
  lp2 <- loglike(y = y, mu = mu, ls = ls, xi = xi, theta = theta,
                 theta.xi = theta.xi, thresh = thresh, alpha = alpha)

  logpost <- lp1 + sum(lp2)

  return(logpost)
}

logpost.beta1.int.grad <- function(beta.int, beta.int.mn, tau, Qb,
                                   beta.time, time, y, ls, xi,
                                   theta, theta.xi, thresh, alpha) {

  d1dmu <- -tau * crossprod(Qb, (beta.int - beta.int.mn))  # crossprod is faster
  # d1dmu <- -tau * (Qb %*% (beta.int - beta.int.mn))

  mu <- beta.int + beta.time * time
  sig <- exp(ls)

  these <- (y > thresh) & !is.na(y) # likelihood changes if y > thresh
  if (abs(xi) <= 1e-4) {
    sig.star <- alpha * sig  # when xi = 0, theta.xi = 1
    t.y <- (theta * exp(-(y - mu) / sig))^(1 / alpha)
    d2dmu <- -t.y / (sig.star)
    d2dmu[these] <- d2dmu[these] + 1 / (sig.star[these])
  } else {
    # theta.xi <- theta^xi
    mu.star  <- mu + sig * ((theta.xi) - 1) / xi
    sig.star <- alpha * sig * (theta.xi)
    xi.star  <- alpha * xi
    t.y <- 1 + xi.star * (y - mu.star) / sig.star
    d2dmu <- - t.y^(-1 / xi.star - 1) / sig.star
    d2dmu[these] <- d2dmu[these] + (xi.star + 1) /
      (sig.star[these] * t.y[these])
  }

  d2dmu[is.na(y)] <- 0

  grad <- d1dmu + d2dmu

  return(grad)
}

logpost.beta1.time <- function(beta.time, beta.time.mn, time, tau, Qb,
                               beta.int, y, ls, xi, theta,
                               theta.xi, thresh, alpha) {
  sig <- exp(ls)
  lp1 <- -0.5 * tau * quad.form(Qb, beta.time - beta.time.mn)

  mu <- beta.int + beta.time * time
  lp2 <- loglike(y = y, mu = mu, ls = ls, xi = xi, theta = theta,
                 theta.xi = theta.xi, thresh = thresh, alpha = alpha)

  logpost <- lp1 + sum(lp2)

  return(logpost)
}

logpost.beta1.time.grad <- function(beta.time, beta.time.mn, time, tau, Qb,
                                    beta.int, y, ls, xi,
                                    theta, theta.xi, thresh, alpha) {
  sig <- exp(ls)
  d1dmu <- -tau * crossprod(Qb, (beta.time - beta.time.mn))  # faster than %*%
  # d1dmu <- -tau * (Qb %*% (beta.time - beta.time.mn))

  mu <- beta.int + beta.time * time

  these <- (y > thresh) & !is.na(y) # likelihood changes if y > thresh
  if (abs(xi) <= 1e-4) {
    sig.star <- alpha * sig  # when xi = 0, theta.xi = 1
    t.y <- (theta * exp(-(y - mu) / sig))^(1 / alpha)
    d2dmu <- -t.y / (sig.star)
    d2dmu[these] <- d2dmu[these] + 1 / (sig.star[these])
  } else {
    # theta.xi <- theta^xi
    mu.star  <- mu + sig * ((theta.xi) - 1) / xi
    sig.star <- alpha * sig * (theta.xi)
    xi.star  <- alpha * xi
    t.y <- 1 + xi.star * (y - mu.star) / sig.star
    d2dmu <- - t.y^(-1 / xi.star - 1) / sig.star
    d2dmu[these] <- d2dmu[these] + (xi.star + 1) /
      (sig.star[these] * t.y[these])
  }

  d2dmu[is.na(y)] <- 0

  grad <- d1dmu + d2dmu * time

  return(grad)
}

logpost.beta1.grad <- function(beta.int, beta.int.mn, tau.int,
                               beta.time, beta.time.mn, time, tau.time,
                               Qb, y, ls, xi, theta, theta.xi, thresh,
                               alpha) {
  sig <- exp(ls)
  # crossprod is generally a little faster than %*%
  d1dbeta.int  <- -tau.int * crossprod(Qb, (beta.int - beta.int.mn))
  d1dbeta.time <- -tau.time * crossprod(Qb, (beta.time - beta.time.mn))
  # d1dmu <- -tau * (Qb %*% (beta.time - beta.time.mn))

  mu <- beta.int + beta.time * time
  these <- (y > thresh) & !is.na(y) # likelihood changes if y > thresh
  if (abs(xi) <= 1e-4) {
    sig.star <- alpha * sig  # when xi = 0, theta.xi = 1
    t.y <- (theta * exp(-(y - mu) / sig))^(1 / alpha)
    d2dmu <- -t.y / (sig.star)
    d2dmu[these] <- d2dmu[these] + 1 / (sig.star[these])
  } else {
    # theta.xi <- theta^xi
    mu.star  <- mu + sig * ((theta.xi) - 1) / xi
    sig.star <- alpha * sig * (theta.xi)
    xi.star  <- alpha * xi
    t.y <- 1 + xi.star * (y - mu.star) / sig.star
    d2dmu <- - t.y^(-1 / xi.star - 1) / sig.star
    d2dmu[these] <- d2dmu[these] + (xi.star + 1) /
      (sig.star[these] * t.y[these])
  }

  d2dmu[is.na(y)] <- 0

  grad.beta.int  <- d1dbeta.int + d2dmu
  grad.beta.time <- d1dbeta.time + d2dmu * time

  results <- list(grad.beta.int = grad.beta.int,
                  grad.beta.time = grad.beta.time)

  return(results)
}

logpost.beta2.int <- function(beta.int, beta.int.mn, tau, beta.time, time,
                              Qb, y, mu, xi, theta, theta.xi, thresh, alpha) {


  lp1 <- -0.5 * tau * quad.form(Qb, beta.int - beta.int.mn)

  ls <- beta.int + beta.time * time
  sig <- exp(ls)
  lp2 <- loglike(y = y, mu = mu, ls = ls, xi = xi, theta = theta,
                 theta.xi = theta.xi, thresh = thresh, alpha = alpha)

  logpost <- lp1 + sum(lp2)

  return(logpost)
}

logpost.beta2.int.grad <- function(beta.int, beta.int.mn, tau,
                                   beta.time, time, Qb, beta, y, mu,
                                   xi, theta, theta.xi, thresh, alpha) {

  d1dls <- -tau * crossprod(Qb, (beta.int - beta.int.mn))  # faster than %*%

  ls <- beta.int + beta.time * time

  sig <- exp(ls)
  these <- (y > thresh) & !is.na(y) # likelihood changes if y > thresh
  if (abs(xi) <= 1e-4) {
    sig.star <- alpha * sig  # when xi = 0, theta.xi = 1
    t.y <- (theta * exp(-(y - mu) / sig))^(1 / alpha)
    d2dls <- -t.y * (y - mu) / (sig.star)
    d2dls[these] <- d2dls[these] - 1 + (y[these] - mu[these]) /
      (sig.star[these])
  } else {
    # theta.xi <- theta^xi
    mu.star  <- mu + sig * ((theta.xi) - 1) / xi
    sig.star <- alpha * sig * (theta.xi)
    xi.star  <- alpha * xi
    y.star <- (y - mu) / sig
    t.y <- (1 + xi * y.star) / theta.xi
    d2dls <- -y.star * t.y^(-1 / xi.star - 1) / (alpha * theta.xi)

    d2dls[these] <- d2dls[these] - 1 + y.star[these] * (xi.star + 1) /
      (t.y[these] * alpha * theta.xi[these])
  }
  d2dls[is.na(y)] <- 0

  grad <- d1dls + d2dls
  return(grad)
}

logpost.beta2.time <- function(beta.time, beta.time.mn, time, tau, beta.int,
                               Qb, y, mu, xi, theta, theta.xi, thresh, alpha) {


  lp1 <- -0.5 * tau * quad.form(Qb, beta.time - beta.time.mn)

  ls <- beta.int + beta.time * time
  sig <- exp(ls)
  lp2 <- loglike(y = y, mu = mu, ls = ls, xi = xi, theta = theta,
                 theta.xi = theta.xi, thresh = thresh, alpha = alpha)

  logpost <- lp1 + sum(lp2)

  return(logpost)
}

logpost.beta2.time.grad <- function(beta.time, beta.time.mn, time, tau,
                                    beta.int, Qb, beta, y, mu,
                                    xi, theta, theta.xi, thresh, alpha) {

  d1dls <- -tau * crossprod(Qb, (beta.time - beta.time.mn))  # faster than %*%

  ls <- beta.int + beta.time * time

  sig <- exp(ls)
  these <- (y > thresh) & !is.na(y) # likelihood changes if y > thresh
  if (abs(xi) <= 1e-4) {
    sig.star <- alpha * sig  # when xi = 0, theta.xi = 1
    t.y <- (theta * exp(-(y - mu) / sig))^(1 / alpha)
    d2dls <- -t.y * (y - mu) / (sig.star)
    d2dls[these] <- d2dls[these] - 1 + (y[these] - mu[these]) /
      (sig.star[these])
  } else {
    # theta.xi <- theta^xi
    mu.star  <- mu + sig * ((theta.xi) - 1) / xi
    sig.star <- alpha * sig * (theta.xi)
    xi.star  <- alpha * xi
    y.star <- (y - mu) / sig
    t.y <- (1 + xi * y.star) / theta.xi
    d2dls <- -y.star * t.y^(-1 / xi.star - 1) / (alpha * theta.xi)

    d2dls[these] <- d2dls[these] - 1 + y.star[these] * (xi.star + 1) /
      (t.y[these] * alpha * theta.xi[these])
  }
  d2dls[is.na(y)] <- 0

  grad <- d1dls + d2dls * time
  return(grad)
}

logpost.beta2.grad <- function(beta.time, beta.time.mn, time, tau,
                                    beta.int, Qb, beta, y, mu,
                                    xi, theta, theta.xi, thresh, alpha) {

  d1dbeta.int  <- -tau * crossprod(Qb, (beta.int - beta.int.mn))
  d1dbeta.time <- -tau * crossprod(Qb, (beta.time - beta.time.mn))

  ls <- beta.int + beta.time * time

  sig <- exp(ls)
  these <- (y > thresh) & !is.na(y) # likelihood changes if y > thresh
  if (abs(xi) <= 1e-4) {
    sig.star <- alpha * sig  # when xi = 0, theta.xi = 1
    t.y <- (theta * exp(-(y - mu) / sig))^(1 / alpha)
    d2dls <- -t.y * (y - mu) / (sig.star)
    d2dls[these] <- d2dls[these] - 1 + (y[these] - mu[these]) /
      (sig.star[these])
  } else {
    # theta.xi <- theta^xi
    mu.star  <- mu + sig * ((theta.xi) - 1) / xi
    sig.star <- alpha * sig * (theta.xi)
    xi.star  <- alpha * xi
    y.star <- (y - mu) / sig
    t.y <- (1 + xi * y.star) / theta.xi
    d2dls <- -y.star * t.y^(-1 / xi.star - 1) / (alpha * theta.xi)

    d2dls[these] <- d2dls[these] - 1 + y.star[these] * (xi.star + 1) /
      (t.y[these] * alpha * theta.xi[these])
  }
  d2dls[is.na(y)] <- 0

  gradbeta.int  <- d1dbeta.int + d2dls
  gradbeta.time <- d1dbeta.time + d2dls * time
  results <- list(grad.beta.int = grad.beta.int,
                  grad.beta.time = grad.beta.time)
  return(results)
}

# if GP prior on mu
logpost.mu <- function(mu, Xb, tau, Qb, y, ls, xi, theta, theta.xi, thresh,
                       alpha) {
  sig <- exp(ls)
  lp1 <- -0.5 * tau * quad.form(Qb, mu - Xb)

  # lp2 <- loglike(y = y, mu = mu, ls = ls, xi = xi, theta = theta,
  #                thresh = thresh, alpha = alpha)
  lp2 <- loglike(y = y, mu = mu, ls = ls, xi = xi, theta = theta,
                 theta.xi = theta.xi, thresh = thresh, alpha = alpha)

  logpost <- lp1 + sum(lp2)

  return(logpost)
}

logpost.mu.grad <- function(mu, Xb, tau, Qb, y, ls, xi, theta, theta.xi,
                            thresh, alpha) {
  sig <- exp(ls)
  d1dmu <- -tau * crossprod(Qb, (mu - Xb))  # crossprod is actually a bit faster

  these <- (y > thresh) & !is.na(y) # likelihood changes if y > thresh
  if (abs(xi) <= 1e-4) {
    sig.star <- alpha * sig  # when xi = 0, theta.xi = 1
    t.y <- (theta * exp(-(y - mu) / sig))^(1 / alpha)
    d2dmu <- -t.y / (sig.star)
    d2dmu[these] <- d2dmu[these] + 1 / (sig.star[these])
  } else {
    # theta.xi <- theta^xi
    mu.star  <- mu + sig * ((theta.xi) - 1) / xi
    sig.star <- alpha * sig * (theta.xi)
    xi.star  <- alpha * xi
    t.y <- 1 + xi.star * (y - mu.star) / sig.star
    d2dmu <- - t.y^(-1 / xi.star - 1) / sig.star
    d2dmu[these] <- d2dmu[these] + (xi.star + 1) /
      (sig.star[these] * t.y[these])
  }

  d2dmu[is.na(y)] <- 0

  grad <- d1dmu + d2dmu

  return(grad)
}

logpost.logsig <- function(ls, Xb, tau, Qb, y, mu, xi, theta, theta.xi, thresh,
                           alpha)
  {
  sig <- exp(ls)
  lp1 <- -0.5 * tau * quad.form(Qb, ls - Xb)

  lp2 <- loglike(y = y, mu = mu, ls = ls, xi = xi, theta = theta,
                 theta.xi = theta.xi, thresh = thresh, alpha = alpha)

  logpost <- lp1 + sum(lp2)

  return(logpost)
}

logpost.logsig.grad <- function(ls, Xb, tau, Qb, y, mu, xi, theta, theta.xi,
                                thresh, alpha) {
  sig <- exp(ls)
  d1dls <- -tau * crossprod(Qb, (ls - Xb))  # crossprod is actually a bit faster

  these <- (y > thresh) & !is.na(y) # likelihood changes if y > thresh
  if (abs(xi) <= 1e-4) {
    sig.star <- alpha * sig  # when xi = 0, theta.xi = 1
    t.y <- (theta * exp(-(y - mu) / sig))^(1 / alpha)
    d2dls <- -t.y * (y - mu) / (sig.star)
    d2dls[these] <- d2dls[these] - 1 + (y[these] - mu[these]) /
      (sig.star[these])
  } else {
    # theta.xi <- theta^xi
    mu.star  <- mu + sig * ((theta.xi) - 1) / xi
    sig.star <- alpha * sig * (theta.xi)
    xi.star  <- alpha * xi
    y.star <- (y - mu) / sig
    t.y <- (1 + xi * y.star) / theta.xi
    d2dls <- -y.star * t.y^(-1 / xi.star - 1) / (alpha * theta.xi)

    d2dls[these] <- d2dls[these] - 1 + y.star[these] * (xi.star + 1) /
      (t.y[these] * alpha * theta.xi[these])
  }
  d2dls[is.na(y)] <- 0

  grad <- d1dls + d2dls
  return(grad)
}

logpost.mu.grad.test <- function(mu, Xb, tau, Qb, y, ls, xi) {
  sig <- exp(ls)
  d1dmu <- as.vector(-tau * Qb %*% (mu - Xb))

  mu.star  <- mu
  sig.star <- sig
  xi.star  <- xi
  t.y <- 1 + xi.star * (y - mu.star) / sig.star
  d2dmu <- (xi.star + 1) / (sig.star * t.y) - t.y^(-1 / xi.star - 1) / sig.star

  grad <- d1dmu + d2dmu
  return(grad)
}

logpost.mu.test <- function(mu, Xb, tau, Qb, y, ls, xi) {
  sig <- exp(ls)
  lp1 <- -0.5 * tau * quad.form(Qb, mu - Xb)

  mu.star  <- mu
  sig.star <- sig
  xi.star  <- xi
  t.y <- (1 + xi.star * (y - mu.star) / sig.star)^(-1 / xi.star)
  lp2 <- (xi.star + 1) * log(t.y) - t.y

  logpost <- lp1 + sum(lp2)

  return(logpost)
}

logpost.logsig.test <- function(mu, Xb, tau, Qb, y, ls, xi) {
  sig <- exp(ls)
  lp1 <- -0.5 * tau * quad.form(Qb, ls - Xb)

  mu.star  <- mu
  sig.star <- sig
  xi.star  <- xi
  lp2 <- dgev(x = y, loc = mu.star, scale = sig.star, shape = xi.star,
              log = TRUE)

  logpost <- lp1 + sum(lp2)

  return(logpost)
}

logpost.logsig.grad.test <- function(mu, Xb, tau, Qb, y, ls, xi) {
  sig <- exp(ls)
  d1dlogsig <- as.vector(-tau * Qb %*% (ls - Xb))

  mu.star  <- mu
  sig.star <- sig
  xi.star  <- xi
  y.star <- (y - mu.star) / sig.star
  t.y <- 1 + xi.star * y.star
  d2dlogsig <- -1 + y.star * ((xi.star + 1) / t.y - t.y^(-1 / xi.star - 1))

  grad <- d1dlogsig + d2dlogsig
  return(grad)
}

dPS.Rcpp <- function(A, alpha, bins) {
  if (is.null(dim(A))) {
    ns <- length(A)
    nt <- 1
    A <- matrix(A, ns, nt)  # turn it into a matrix
  } else {
    ns <- nrow(A)
    nt <- ncol(A)
  }

  mid.points <- bins$MidPoints
  bin.width  <- bins$BinWidth

  results <- dPSCPP(A=A, alpha=alpha, mid_points=mid.points,
                    bin_width=bin.width)
  return(results)
}

ld.old <- function(u, A, alpha){
  psi <- pi * u
  c <- (sin(alpha * psi) / sin(psi))^(1 / (1 - alpha))
  c <- c * sin((1 - alpha) * psi) / sin(alpha * psi)
  logd <- log(alpha) - log(1 - alpha) - (1 / (1 - alpha)) * log(A) +
    log(c) - c * (1/A^(alpha / (1 - alpha)))
  exp(logd)
}

dPS.old <- function(A, alpha, bins){
  l <- 0
  for(j in 1:bins$npts){
    l <- l + bins$BinWidth[j] * ld(bins$MidPoints[j], A, alpha)
  }
  l <- ifelse(A > 0, log(l), -Inf)
  return(l)
}


######################################################
##########    FUNCTIONS USED FOR PREDICTION  #########
######################################################

rGEV <- function(n, mu, sig, xi) {
  tau <- runif(n)
  x <- -1 / log(tau)
  x <- x^(xi) - 1
  x <- mu + sig * x / xi
  return(x)
}

proj.beta <- function(B, d12, d22, S11inv, tau, logrho) {
  #B<-n-vector of observed beta (minus the mean)
  ns <- nrow(d22)
  rho <- exp(logrho)

  S22 <- exp(-d22 / rho) / tau
  S12 <- exp(-d12 / rho) / tau
  S11inv <- S11inv * tau

  P2 <- S12 %*% S11inv
  P1 <- S22 - S12 %*% S11inv %*% t(S12)
  P1 <- t(chol(P1))

  Bnew <- P2 %*% B + P1 %*% rnorm(ns)
  return(Bnew)
}

######################################################
####  FUNCTION TO COMPUTE THE RANDOM EFFECTS  ########
######################################################
# this is theta^(1 / alpha)
make.theta <- function(FAC, logs, alpha) {
  #theta is nxnF
  #s is nFxnt
  #alpha in (0,1)
  if (length(logs) == 1) {
    xxx <- (FAC^(1 / alpha)) * exp(logs)
  }
  if (length(logs) > 1) {
    xxx <- (FAC^(1 / alpha)) %*% exp(logs)
  }
  return(xxx)
}

######################################################
########  2d gaussian kernel basis functions  ########
######################################################
stdKern <- function(w, single = FALSE) {
  if (single) { K <- w / sum(w) }
  if (!single) { K <- sweep(w, 1, rowSums(w), "/") }
  return(K)
}

make.kern <- function(d2, logrho) {
  rho2 <- exp(logrho)^2
  w <- exp(-0.5 * d2 / rho2)
  return(w)
}

getW <- function(dw2, rho) {
  w <- stdW(makeW(dw2 = dw2, rho = rho))
  return(w)
}

# get the kernel weighting
makeW <- function(dw2, rho) {
  w <- exp(-0.5 * dw2 / (rho^2))

  return(w)
}

# standardize the kernel weights
stdW <- function(x, single = FALSE) {
  if (single) {x <- x / sum(x)}
  if (!single) {x <- sweep(x, 1, rowSums(x), "/")}
  return(x)
}


add.basis.X <- function(X, B, time.interact = FALSE) {
  ns <- dim(X)[1]
  nt <- dim(X)[2]
  np <- dim(X)[3]
  if (time.interact) {
    nB <- dim(B)[2] * 2
  } else {
    nB <- dim(B)[2]
  }

  # create a new X array that will be big enough to hold the basis functions
  newX <- array(0, dim = c(ns, nt, np + nB))

  # copy over old X information
  newX[, , 1:np] <- X

  if (time.interact) {
    for (t in 1:nt) {
      B.interact <- B * X[, t, 2]
      newX[, t, (np + 1):(np + nB)] <- cbind(B, B.interact)
    }
  } else {
    for (t in 1:nt) {
      newX[, t, (np + 1):(np + nB)] <- B
    }
  }

  return(newX)
}

rep.basis.X <- function(X, newB, time.interact = FALSE) {
  ns <- dim(X)[1]
  nt <- dim(X)[2]
  np <- dim(X)[3]
  nB <- dim(newB)[2]
  if (time.interact) {
    start <- np - 2 * nB + 1
  } else {
    start <- np - nB + 1
  }
  end   <- np


  if (time.interact) {
    for (t in 1:nt) {
      newB.interact <- newB * X[, t, 2]
      X[, t, start:end] <- cbind(newB, newB.interact)
    }
  } else {
    for (t in 1:nt) {
      X[, t, start:end] <- newB
    }
  }

  return(X)
}

######################################################
########            OTHER FUNCTIONS        ###########
######################################################

get.level <- function(logs, cuts){
  if (is.matrix(logs)) {
    level <- getLevelCPP(A = logs, cuts = cuts)
  } else if (length(logs) > 1) {
    A <- as.matrix(logs, length(logs), 1)
    level <- as.vector(getLevelCPP(A = A, cuts = cuts))
  } else {
    level <- sum(logs > cuts) + 1
  }

  return(level)
}

get.level.mat <- function(logs, cuts) {
  level <- logs
  for (i in 1:nrow(logs)) { for (j in 1:ncol(logs)) {
    level[i, j] <- sum(logs[i, j] > cuts) + 1
  }}

  return(level)
}

logdet <- function(X){
  determinant(X)$modulus
}

trunc <- function(x, eps = 0.1){
  x <- ifelse(x < eps, eps, x)
  x <- ifelse(x > 1 - eps, 1 - eps, x)
  return(x)
}

getXBeta <- function(X, beta) {
  np <- length(beta)
  if (np != dim(X)[3]) {
    stop("X has the wrong dimensions")
  }
  XBeta <- 0
  for (p in 1:np) {
    XBeta <- XBeta + X[, , p] * beta[p]
  }

  return(XBeta)
}

# SSE for row of Y - EC
SSE.B <- function(B1, B2, B.star, Y, alpha, lambda = 1000){

  BB  <- B1^(1 / alpha)
  B2  <- B.star
  EC  <- sweepC2plus(X = B2, y = BB)
  EC  <- rowSumsC(EC^alpha)

  # penalty term is to make sure that the bases sum to 1
  sse <- sum((Y - EC)^2, na.rm = TRUE) + lambda * (sum(B1) - 1)^2

  return(sse)
}

SSE.B.grad <- function(B1, B.star, Y, alpha, lambda = 1000, exclude = 1){

  B1.star <- B1^(1 / alpha)
  B2   <- B.star

  BB <- sweepC2plus(X = B2, y = B1.star)
  EC0  <- rowSumsC(BB^alpha)

  EC1  <- BB^(alpha - 1)
  EC1  <- sweepC2times(EC1, B1.star / B1)
  EC1  <- sweepC1times(EC1, Y - EC0)

  grad <- -2 * colSums(EC1, na.rm = TRUE) +
    2 * lambda * (sum(B1) - 1)

  return(grad)
}


make.EC  <- function(B, alpha){
  Ba    <- B^(1 / alpha)
  EC    <- NULL
  for(j in 1:nrow(B)){
    BB <- sweep(Ba, 2, Ba[j, ], "+")
    EC <- cbind(EC, rowSums(BB^alpha))
  }

  return(EC)
}

SSE.rhoalpha <- function(rho, dw2, Y, alpha) {
  w <- getW(rho = rho, dw2 = dw2)
  w <- w^(1 / alpha)

  n <- ncol(Y)
  EC <- getECRhoAlphaC(w = w, alpha = alpha)
  diag(EC) <- NA

  if (any(is.nan(EC))) {
    return(Inf)
  }

  sse <- sum((Y - EC)^2, na.rm = TRUE)

  return(sse)
}

getGPSS <- function(Qb, param, Xb) {
  # get the sum of squares for a gaussian process with prec matrix Qb
  if (is.matrix(Xb)) {  # this takes less time than the slower version
    SS <- diag(quad.form(Qb, param - Xb))
  } else {
    SS <- quad.form(Qb, param - Xb)
  }
  return(SS)
}

getGPSS.slow <- function(Qb, param, Xb) {
  # get the sum of squares for a gaussian process with prec matrix Qb
  if (is.matrix(Xb)) {
    nt <- ncol(Xb)
    SS <- rep(0, nt)
    for (t in 1:nt) {
      SS[t] <- quad.form(Qb, param[, t] - Xb[, t])
    }
  } else {
    SS <- quad.form(Qb, param - Xb)
  }
  return(SS)
}

# Performs kernel smoothing of the extremal coefficient matrix.
Ksmooth <- function(ECmat, s = NULL, bw = NULL){

  n           <- nrow(ECmat)
  diag(ECmat) <- 0
  E1          <- ifelse(ECmat == 0, 0, 1)
  if (is.null(s)) {s <- 1:n}
  if (is.null(bw)) {bw <- 2 * min(dist(s))}

  d2       <- as.matrix(dist(s) / bw)^2
  W        <- exp(-d2)
  diag(W)  <- 0

  num      <- W %*% ECmat %*% W
  den      <- W %*% E1 %*% W

  ECsmooth <- num / den

  return(ECsmooth)
}

#### Plotting
theme_clean <- function(base_size = 12) {
  require(grid)
  theme_grey(base_size) %+replace%
    theme(
      axis.title      =   element_blank(),
      axis.text       =   element_blank(),
      panel.background    =   element_blank(),
      panel.grid      =   element_blank(),
      axis.ticks.length   =   unit(0,"cm"),
      panel.margin    =   unit(0.5,"lines"),
      plot.margin     =   unit(c(0.5,0.5,0.5,1),"lines"),
      complete = TRUE
    )
}

map.ga.ggplot <- function(Y, counties = NULL, main = "", fill.legend = "",
                          midpoint = NULL, limits = NULL) {
  require(ggplot2)
  require(maps)
  if (is.null(midpoint)) {
    midpoint <- 1.5
  }
  georgia <- map("county", "georgia", fill = TRUE, col = "transparent",
                 plot = FALSE)
  subregion <- sapply(strsplit(georgia$names, ","), function(x) x[2])
  county_map <- map_data(map = "county", region = "georgia")

  # a hack in case the data is in a different order than subregion
  if (is.null(counties)) {
    cat("To guarantee accurate maps, it is recommended to include",
        "a list of counties.")
    basis <- data.frame(Y, subregion)
  } else {
    basis <- data.frame(Y, subregion = counties)
  }
  extcoef_map <- merge(county_map, basis, all.x = TRUE)

  if (is.null(limits)) {
    limits <- c(min(Y), max(Y))
  }

  # using fill = Y because that's the column of extcoef_map with the actual data
  p <- ggplot(extcoef_map, aes(x = long, y = lat, group = group, fill = Y))
  p <- p + geom_polygon(colour = "grey", aes(fill = Y))
  p <- p + expand_limits(x = extcoef_map$long, y = extcoef_map$lat)
  p <- p + coord_map("polyconic")
  p <- p + ggtitle(main)  # make the title
  p <- p + labs(fill = fill.legend)
  p <- p + scale_fill_gradient2(low = "dodgerblue4", high = "firebrick4",
                                mid = "#ffffff", midpoint = midpoint,
                                limits = limits)
  p <- p + theme_clean()
  return(p)
}

map.sc.ggplot <- function(Y, main = "", fill.legend = "", midpoint = NULL) {
  require(ggplot2)
  require(maps)
  if (is.null(midpoint)) {
    midpoint <- 1.5
  }
  sc <- map("county", "south carolina", fill = TRUE, col = "transparent",
            plot = FALSE)
  subregion <- sapply(strsplit(sc$names, ","), function(x) x[2])
  county_map <- map_data(map = "county", region = "south carolina")

  basis <- data.frame(Y, subregion)
  extcoef_map <- merge(county_map, basis, all.x = TRUE)

  # using fill = Y because that's the column of extcoef_map with the actual data
  p <- ggplot(extcoef_map, aes(x = long, y = lat, group = group, fill = Y))
  p <- p + geom_polygon(colour = "grey", aes(fill = Y))
  p <- p + expand_limits(x = extcoef_map$long, y = extcoef_map$lat)
  p <- p + coord_map("polyconic")
  p <- p + labs(title = main, fill = fill.legend)
  p <- p + scale_fill_gradient2(low = "dodgerblue4", high = "firebrick4",
                                mid = "#ffffff", midpoint = midpoint)
  p <- p + theme_clean()
  return(p)
}

################################################################
# Arguments:
#   preds(iters, npreds): mcmc predictions for test site/day
#   probs(nprobs): sample quantiles for scoring
#   validate(npreds): validation data
#
# Returns:
#   score(nprobs): a single quantile score per quantile
################################################################
QuantScore <- function(preds, probs, validate) {

  #   nt <- ncol(validate)  # number of prediction days
  #   np <- nrow(validate)  # number of prediction sites
  npreds <- length(validate)
  nprobs <- length(probs)  # number of quantiles to find quantile score

  # we need to know the predicted quantiles for each site and day in the
  # validation set.
  # nprobs x npreds
  pred.quants <- apply(preds, 2, quantile, probs=probs, na.rm=T)

  scores.sites <- matrix(NA, nprobs, npreds)

  for (q in 1:nprobs) {
    diff <- pred.quants[q, ] - validate
    i <- diff >= 0  # diff >= 0 means qhat is larger
    scores.sites[q, ] <- 2 * (i - probs[q]) * diff
  }

  scores <- apply(scores.sites, 1, mean, na.rm=T)

  return(scores)
}

BrierScore <- function(preds, validate, thresh) {
  # iters <- nrow(post.prob)
  # np    <- ncol(post.prob)

  # scores <- rep(NA, iters)
  # for (i in 1:iters) {
  #   scores[i] <- mean((validate - post.prob[i, ])^2)
  # }
  probs <- rep(NA, length(validate))
  for (i in 1:length(validate)) {
    probs[i] <- mean(preds[, i] > thresh[i])
  }
  score <- mean(((validate >= thresh) - probs)^2)

  return(score)
}


get.pw.ec <- function(Y, nq = 100, qlim = c(0, 1), site.idx = 1,
                      verbose = FALSE, update = NULL) {
  # get the pairwise chi as an average over nq quantiles
  # between qlim[1] and qlim[2]
  # if qlim[2] == 1, then we'll set it to the max quantile for the two sites

  if (site.idx == 2) {  # each column represents a site
    Y <- t(Y)  # transform to rows
  }

  ns <- nrow(Y)
  nt <- ncol(Y)

  if (is.null(update)) {
    update <- floor(ns / 4)
  }

  ec <- matrix(0, ns, ns)
  eps <- .Machine$double.eps^0.5

  qlims <- matrix(0, nrow = (ns * ns - ns) / 2 + ns, ncol = 2)
  qlim.idx <- 1
  these <- ns
  for (i in 1:ns) {
    for (j in i:ns) {

      # only want to include years which have both observations
      these.ij   <- which(colSums(is.na(Y[c(i, j), ])) == 0)
      U <- Y[c(i, j), these.ij]
      for (k in 1:2) {
        U[k, ] <- rank(U[k, ]) / (length(these.ij) + 1)
      }
      colmax.ij  <- apply(U, 2, max)
      min.ij     <- max(min(U), qlim[1]) + eps
      max.ij     <- min(max(U), qlim[2]) - eps
      if (max.ij < min.ij) {
        cat("  i:", i, "j:", j, "\n")
      }
      quantiles  <- seq(min.ij, max.ij, length = nq)
      if (max(quantiles) == 1) {
        cat("  i:", i, "j:", j, "\n")
      }
      qhat <- rep(NA, length(quantiles))
      Q.ij <- rep(NA, length(quantiles))
      for (q in 1:length(quantiles)) {
        qhat[q] <- mean(U[1, ] < quantiles[q])
        Q.ij[q] <- mean(colmax.ij < quantiles[q])
      }

      # we're using qhat here because we want to make sure that the quantile by
      # which we divide is only as precise as the data will allow.
      ec[i, j] <- ec[j, i] <- mean(log(Q.ij) / log(qhat))
      ec[i, j] <- ec[j, i] <- min(ec[i, j], 2)  # keep it within the range
      ec[i, j] <- ec[j, i] <- max(ec[i, j], 1)  # keep it within the range

      qlims[qlim.idx, ] <- c(min.ij, max.ij)
      qlim.idx <- qlim.idx + 1
    }
    if (verbose & (i %% update == 0)) {
      cat("  Finished i =", i, "\n")
    }
  }

  return(list(ec = ec, qlims = qlims))
}


get.pw.ec.fmado <- function(Y, thresh = NULL, thresh.quant = FALSE,
                            qlim = c(0, 1), site.idx = 1) {

  if (site.idx == 2) {  # each column represents a site
    Y <- t(Y)  # transform to rows
  }

  ns <- nrow(Y)
  nt <- ncol(Y)

  # we want a way to handle POT and max-stable in the same funciton.
  # if the data are max-stable, we can just set the threshold at -Inf
  # for all sites.
  if (is.null(thresh)) {
    thresh <- rep(-Inf, ns)
  } else if (length(thresh) == 1) {
    cat("\t Using the same threshold for all sites \n")
    thresh <- rep(thresh, ns)
  } else if (length(thresh) != ns) {
    stop("If defined, thresh must be of length: 1, or length: ns.")
  }

  if (!thresh.quant) {  # find empirical cdf at threshold values
    for (i in 1:ns) {
      thresh[i] <- mean(Y[i, ] < thresh[i], na.rm = TRUE)
    }
  }

  # get values of empirical cdf for Y
  # need t(apply) because apply gives back nt x ns
  Y <- t(apply(Y, 1, rank, na.last = "keep")) / (rowSums(is.finite(Y)) + 1)

  # shift and scale to account for threshold
  Y <- (Y - thresh) / (1 - thresh)
  Y[Y <= 0] <- 0

  if (is.null(update)) {
    update <- floor(ns / 4)
  }

  # qlims <- matrix(0, nrow = (ns * ns - ns) / 2 + ns, ncol = 2)
  # qlim.idx <- 1
  fmado <- madogramCPP(data = Y)
  fmado[fmado >= 1 / 6] <- 1 / 6
  ec <- (1 + 2 * fmado) / (1 - 2 * fmado)

  return(list(ec = ec, fmadogram = fmado))
}

bspline.2d <- function(s, scale = TRUE, df.x, df.y) {
  ns <- nrow(s)

  if (scale) {
    s.scale <- min(diff(range(s[, 1])), diff(range(s[, 2])))
    s[, 1] <- (s[, 1] - min(s[, 1])) / s.scale
    s[, 2] <- (s[, 2] - min(s[, 2])) / s.scale
  }

  B.x <- bs(s[, 1], df = df.x, Boundary.knots = c(-0.1, 1.1))
  B.y <- bs(s[, 2], df = df.y, Boundary.knots = c(-0.1, 1.1))

  B <- matrix(NA, nrow = ns, ncol = df.x * df.y)
  for (i in 1:ncol(B.x)) {
    for (j in 1:ncol(B.y)) {
      B[, (i - 1) * df.y + j] <- B.x[, i] * B.y[, j]
    }
  }

  # if the basis has no weight for any sites, remove from group
  keep.bases <- which(colSums(B) > 0)
  B <- B[, keep.bases]

  return(B)
}

mrl.plot <- function (data, umin = min(data), umax = max(data) - 0.1,
                      conf = 0.95, nint = 100, xlab = NULL, ylab = NULL,
                      main = NULL) {
  x <- xu <- xl <- numeric(nint)
  u <- seq(umin, umax, length = nint)
  for (i in 1:nint) {
    data <- data[data > u[i]]
    x[i] <- mean(data - u[i])
    sdev <- sqrt(var(data))
    n <- length(data)
    xu[i] <- x[i] + (qnorm((1 + conf)/2) * sdev)/sqrt(n)
    xl[i] <- x[i] - (qnorm((1 + conf)/2) * sdev)/sqrt(n)
  }

  if (is.null(xlab)) {
    xlab <- "u"
  }
  if (is.null(ylab)) {
    ylab <- "Mean Excess"
  }
  if (is.null(main)) {
    main <- ""
  }
  plot(u, x, type = "l", xlab = xlab, ylab = ylab, main = main,
       ylim = c(min(xl[!is.na(xl)]), max(xu[!is.na(xu)])),
       cex.axis = 1.5, cex.main = 1.5, cex.lab = 1.5, lwd = 1.25)
  lines(u[!is.na(xl)], xl[!is.na(xl)], lty = 2)
  lines(u[!is.na(xu)], xu[!is.na(xu)], lty = 2)
}

#############################################################:
###           OTHER FUNCTIONS USED IN THE MCMC            ###:
#############################################################:
grad_loglike_betamu <- function(beta1, X1, y, theta, ls, xi, thresh,
                                alpha) {
  mu <- 0
  p.mu <- dim(X1)[3]
  for (j in 1:p.mu) {
    mu <- mu + X1[, , j] * beta1[j]
  }

  sigma    <- exp(ls)
  mu_star  <- mu + sigma * ((theta^xi) - 1) / xi
  sig_star <- alpha * sigma * (theta^xi)
  xi_star  <- alpha * xi

  tx_star  <- (1 + xi_star * (y - mu_star) / sig_star)

  this.grad <- -(tx_star^(-1 / xi_star - 1)) / sig_star +
    (y > thresh) * (xi_star + 1) / (sig_star * tx_star)
  this.grad[is.na(y)] <- 0  # for the missing data

  grad <- rep(0, p.mu)
  for (j in 1:p.mu) {
    grad[j] <- sum(this.grad * X1[, , j])
  }

  return(grad)
}

grad_logpost_betamu <- function(beta1, beta.mu, beta.sd, X1, y, theta, ls,
                                xi, thresh, alpha) {
  mu <- 0
  p.mu <- dim(X1)[3]
  for (j in 1:p.mu) {
    mu <- mu + X1[, , j] * beta1[j]
  }

  sig <- exp(ls)

  tx  <- (1 + xi * (y - mu) / sig)

  this.grad <- - (tx^(-1 / (alpha * xi) - 1)) * theta^(1 / alpha) /
    (alpha * sig) +
    (y > thresh) * (alpha * xi + 1) / (alpha * tx * sig)
  this.grad[is.na(y)] <- 0  # for the missing data

  grad <- rep(0, p.mu)
  for (j in 1:p.mu) {
    grad[j] <- sum(this.grad * X1[, , j]) - (beta1[j] - beta.mu) / beta.sd^2
  }

  return(grad)
}

hess_logpost_betamu <- function(beta1, beta.mu, beta.sd, X1, y, theta, ls,
                                xi, thresh, alpha) {
  mu <- 0
  p.mu <- dim(X1)[3]
  for (j in 1:p.mu) {
    mu <- mu + X1[, , j] * beta1[j]
  }

  sigma <- exp(ls)
  tx <- (1 + xi * (y - mu) / sigma)

  this.hess <- - (tx^(-1 / (alpha * xi) - 2) * (alpha * xi + 1) *
                    theta^(1 / alpha)) / (alpha^2 * sigma^2) +
    (y > thresh) * (alpha * xi + 1) * xi / (alpha * tx^2 * sigma^2)
  this.hess[is.na(y)] <- 0  # for missing data

  hess <- matrix(0, p.mu, p.mu)
  for (j in 1:p.mu) {
    for (i in j:p.mu) {
      hess[i, j] <- hess[j, i] <- sum(this.hess * X1[, , i] * X1[, , j])
    }
  }

  hess <- hess - diag(1 / beta.sd^2, p.mu)  # account for prior

  return(hess)
}

hess_loglike_betamu <- function(beta1, X1, y, theta, ls, xi, thresh,
                                alpha) {
  mu <- 0
  p.mu <- dim(X1)[3]
  for (j in 1:p.mu) {
    mu <- mu + X1[, , j] * beta1[j]
  }

  # print(mean(mu))
  sigma    <- exp(ls)
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
      hess[i, j] <- hess[j, i] <- sum(this.hess * X1[, , i] * X1[, , j])
    }
  }

  return(hess)
}

loglike_mu <- function(beta1, X1, y, theta, ls, xi, thresh, alpha) {
  mu <- 0
  p.mu <- dim(X1)[3]
  for (j in 1:p.mu) {
    mu <- mu + X1[, , j] * beta1[j]
  }

  ll <- loglike(y = y, theta = theta, mu = mu, ls = ls, xi = xi,
                thresh = thresh, alpha = alpha)

  return(sum(ll[!is.na(y)]))
}

logpost_mu <- function(beta1, beta.mu, beta.sd, X1, y, theta, ls, xi,
                       thresh, alpha) {
  mu <- 0
  p.mu <- dim(X1)[3]
  for (j in 1:p.mu) {
    mu <- mu + X1[, , j] * beta1[j]
  }

  ll <- loglike(y = y, theta = theta, mu = mu, ls = ls, xi = xi,
                thresh = thresh, alpha = alpha)

  lp <- sum(dnorm(beta1, beta.mu, beta.sd, log = TRUE))

  return(sum(ll[!is.na(y)]) + lp)
}

grad_loglike_betasig <- function(beta2, X2, y, theta, mu, xi, thresh,
                                 alpha) {
  ls <- 0
  p.sig <- dim(X2)[3]
  for (j in 1:p.sig) {
    ls <- ls + X2[, , j] * beta2[j]
  }

  sigma    <- exp(ls)
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
    grad[j] <- sum(this.grad * X2[, , j])
  }

  return(grad)
}

grad_logpost_betasig <- function(beta2, beta.mu, beta.sd, X2, y, theta, mu,
                                 xi, thresh, alpha) {
  ls <- 0
  p.sig <- dim(X2)[3]
  for (j in 1:p.sig) {
    ls <- ls + X2[, , j] * beta2[j]
  }

  sigma <- exp(ls)
  res <- y - mu

  tx <- (1 + xi * res / sigma)

  this.grad <- - theta^(1 / alpha) * res * tx^(-1 / (alpha * xi) - 1) /
    (alpha * sigma) +
    (y > thresh) * ((alpha * xi + 1) * res / (alpha * tx * sigma) - 1)

  this.grad[is.na(y)] <- 0

  grad <- rep(0, p.sig)
  for (j in 1:p.sig) {
    grad[j] <- sum(this.grad * X2[, , j]) - (beta2[j] - beta.mu) / beta.sd^2
  }

  return(grad)
}

hess_logpost_betasig <- function(beta2, beta.mu, beta.sd, X2, y, theta, mu,
                                 xi, thresh, alpha) {
  ls <- 0
  p.sig <- dim(X2)[3]
  for (j in 1:p.sig) {
    ls <- ls + X2[, , j] * beta2[j]
  }

  sigma <- exp(ls)
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
      hess[i, j] <- hess[j, i] <- sum(this.hess * X2[, , i] * X2[, , j])
    }
  }

  hess <- hess - diag(1 / beta.sd^2, p.sig)  # account for prior

  return(hess)
}

hess_loglike_betasig <- function(beta2, X2, y, theta, mu, xi, thresh,
                                 alpha) {
  ls <- 0
  p.sig <- dim(X2)[3]
  for (j in 1:p.sig) {
    ls <- ls + X2[, , j] * beta2[j]
  }

  sigma    <- exp(ls)
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
      hess[i, j] <- hess[j, i] <- sum(this.hess * X2[, , i] * X2[, , j])
    }
  }

  return(hess)
}

loglike_sig <- function(beta2, X2, y, theta, mu, xi, thresh, alpha) {
  ls <- 0
  p.sig <- dim(X2)[3]
  for (j in 1:p.sig) {
    ls <- ls + X2[, , j] * beta2[j]
  }

  ll <- loglike(y = y, theta = theta, mu = mu, ls = ls, xi = xi,
                thresh = thresh, alpha = alpha)

  return(sum(ll))
}

logpost_sig <- function(beta2, beta.mu, beta.sd, X2, y, theta, mu, xi,
                        thresh, alpha) {
  ls <- 0
  p.sig <- dim(X2)[3]
  for (j in 1:p.sig) {
    ls <- ls + X2[, , j] * beta2[j]
  }

  ll <- loglike(y = y, theta = theta, mu = mu, ls = ls, xi = xi,
                thresh = thresh, alpha = alpha)

  lp <- sum(dnorm(beta2, beta.mu, beta.sd, log = TRUE))
  return(sum(ll[!is.na(y)]) + lp)
}

# ###########  PS functions  ############
#
# get.level <- function(A, cuts){
#   lev <- A * 0 + 1
#   for (j in 1:length(cuts)) {
#     lev <- ifelse(A > cuts[j], j + 1, lev)
#   }
#   return(lev)
# }