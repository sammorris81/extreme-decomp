# Rcpp functions
if (!exists("dPSCPP")) {
  sourceCpp(file = "llps.cpp")
}

if (!exists("ifelsematCPP")) {
  sourceCpp(file = "ifelse.cpp")
}

if (!exists("madogramCPP")) {
  sourceCpp(file = "ec.cpp")
}

#############################################################
########   FUNCTIONS TO COMPUTE INITIAL VALUES    ###########
#############################################################

get.inits.mu <- function(y, logsig = 0, xi = 0.1){
  m <- median(y, na.rm = TRUE)
  mu <- m - exp(logsig) * (log(2)^(-xi) - 1) / xi
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
########           GEV FUNCTIONS           ###########
######################################################

# loglike <- function(y, mu, logsig, xi, theta, alpha) {
#   missing <- is.na(y)
#   theta <- theta^alpha
#   mu.star <- mu + exp(logsig) * (theta^xi - 1) / xi
#   sig.star <- alpha * exp(logsig) * (theta^xi)
#   xi.star <- alpha * xi
#   ttt <- (1 + xi.star * (y - mu.star) / sig.star)^(-1 / xi.star)
#   lll <- -log(sig.star) + (xi.star + 1) * log(ttt) - ttt
#   lll[missing] <- 0
#   return(lll)
# }

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
    B.interact <- B * X[, t, 2]
    newX[, t, (np + 1):(np + nB)] <- cbind(B, B.interact)
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
  sum(logs > cuts) + 1
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

rtnorm <- function(mn, sd = 0.25, fudge = 0){
  upper <- pnorm(1 - fudge, mn, sd)
  lower <- pnorm(fudge, mn, sd)
  if (is.matrix(mn)) {
    U <- matrix(runif(prod(dim(mn)), lower, upper), dim(mn)[1], dim(mn)[2])
  }
  if (!is.matrix(mn)) {
    U <- runif(length(mn), lower, upper)
  }
  return(qnorm(U, mn, sd))
}

dtnorm <- function(y, mn, sd = 0.25, fudge = 0){
  upper <- pnorm(1 - fudge, mn, sd)
  lower <- pnorm(fudge, mn, sd)
  l <- dnorm(y, mn, sd, log = TRUE) - log(upper - lower)
  return(l)
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

# update candidate standard deviation
mhUpdate <- function(acc, att, MH, nattempts = 50,
                     target.min = 0.3, target.max = 0.6,
                     lower = 0.8, higher = 1.2) {
  acc.rate     <- acc / att
  these.update <- att > nattempts
  these.low    <- (acc.rate < target.min) & these.update
  these.high   <- (acc.rate > target.max) & these.update

  MH[these.low]  <- MH[these.low] * lower
  MH[these.high] <- MH[these.high] * higher

  acc[these.update] <- 0
  att[these.update] <- 0

  results <- list(acc=acc, att=att, MH=MH)
  return(results)
}

dPS.Rcpp <- function(a, alpha, mid.points, bin.width) {
  if (is.null(dim(a))) {
    ns <- length(a)
    nt <- 1
    a <- matrix(a, ns, nt)  # turn it into a matrix
  } else {
    ns <- nrow(a)
    nt <- ncol(a)
  }

  results <- dPSCPP(a=a, alpha=alpha, mid_points=mid.points,
                    bin_width=bin.width)
  return(results)
}

#### Logposterior - needs to be done for each timepoint
logpost.mu <- function(mu, Xb, tau, Qb, y, logsig, xi) {
  sig <- exp(logsig)
  lp1 <- -0.5 * tau * quad.form(Qb, mu - Xb)
  # lp1 <- 0

  mu.star  <- mu
  sig.star <- sig
  xi.star  <- xi
  t.y <- (1 + xi.star * (y - mu.star) / sig.star)^(-1 / xi.star)
  lp2 <- (xi.star + 1) * log(t.y) - t.y
  # lp2 <- 0

  logpost <- lp1 + sum(lp2)

  return(logpost)
}

logpost.mu.grad <- function(mu, Xb, tau, Qb, y, logsig, xi) {
  sig <- exp(logsig)
  d1dmu <- as.vector(-tau * Qb %*% (mu - Xb))
  # d1dmu <- 0

  mu.star  <- mu
  sig.star <- sig
  xi.star  <- xi
  t.y <- 1 + xi.star * (y - mu.star) / sig.star
  d2dmu <- (xi.star + 1) / (sig.star * t.y) - t.y^(-1 / xi.star - 1) / sig.star
  # d2dmu <- 0

  grad <- d1dmu + d2dmu
  return(grad)
}

#### Logposterior - needs to be done for each timepoint
logpost.logsig <- function(mu, Xb, tau, Qb, y, logsig, xi) {
  sig <- exp(logsig)
  lp1 <- -0.5 * tau * quad.form(Qb, logsig - Xb)
  # lp1 <- 0

  mu.star  <- mu
  sig.star <- sig
  xi.star  <- xi
  # t.y <- (1 + xi.star * (y - mu.star) / sig.star)^(-1 / xi.star)
  # lp2 <- -log(sig.star) + (xi.star + 1) * log(t.y) - t.y
  lp2 <- dgev(x = y, loc = mu.star, scale = sig.star, shape = xi.star,
              log = TRUE)
  # lp2 <- 0

  logpost <- lp1 + sum(lp2)

  return(logpost)
}

logpost.logsig.grad <- function(mu, Xb, tau, Qb, y, logsig, xi) {
  sig <- exp(logsig)
  d1dlogsig <- as.vector(-tau * Qb %*% (logsig - Xb))
  # d1dlogsig <- 0

  mu.star  <- mu
  sig.star <- sig
  xi.star  <- xi
  y.star <- (y - mu.star) / sig.star
  t.y <- 1 + xi.star * y.star
  d2dlogsig <- -1 + y.star * ((xi.star + 1) / t.y - t.y^(-1 / xi.star - 1))

  grad <- d1dlogsig + d2dlogsig
  return(grad)
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


################################################################
# Arguments:
#   preds(iters, npreds): mcmc predictions for test site/day
#   probs(nprobs): sample quantiles for scoring
#   validate(npreds): validation data
#   thresh(npreds): threshold for the site at which the prediction
#                    is made
#
# Returns:
#   scores: list of scores (bs and qs)
#   scores$bs: brier score for exceeding the scoring threshold
#   scores$qs: quantile score for validation data that exceeds thresh
#              at the site
################################################################
Score <- function(preds, probs, validate, thresh) {
  these.qs <- validate > thresh

  # find the Brier scores for all sites

  # only get the quantile scores for these.qs
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
  fmado <- ifelse(fmado >= 1 / 6, 1 / 6, fmado)
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