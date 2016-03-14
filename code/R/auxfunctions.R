# Rcpp functions
if (!exists("dPSCPP")) {
  sourceCpp(file = "llps.cpp")
}

if (!exists("ifelsematCPP")) {
  sourceCpp(file = "ifelse.cpp")
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

loglike <- function(y, mu, logsig, xi, theta, alpha) {
  missing <- is.na(y)
  theta <- theta^alpha
  mu.star <- mu + exp(logsig) * (theta^xi - 1) / xi
  sig.star <- alpha * exp(logsig) * (theta^xi)
  xi.star <- alpha * xi
  ttt <- (1 + xi.star * (y - mu.star) / sig.star)^(-1 / xi.star)
  lll <- -log(sig.star) + (xi.star + 1) * log(ttt) - ttt
  lll[missing] <- 0
  return(lll)
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



######################################################
########            OTHER FUNCTIONS        ###########
######################################################

get.level <- function(logs, cuts){
  sum(logs > cuts) + 1
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

# update candidate standard deviation
mhUpdate <- function(acc, att, mh, nattempts = 50, lower = 0.8, higher = 1.2) {
  acc.rate     <- acc / att
  these.update <- att > nattempts
  these.low    <- (acc.rate < 0.25) & these.update
  these.high   <- (acc.rate > 0.50) & these.update
  
  mh[these.low]  <- mh[these.low] * lower
  mh[these.high] <- mh[these.high] * higher
  
  acc[these.update] <- 0
  att[these.update] <- 0
  
  results <- list(acc=acc, att=att, mh=mh)
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
      panel.margin    =   unit(0,"lines"),
      plot.margin     =   unit(c(0,0,0,0),"lines"),
      complete = TRUE
    )
}

map.ga.ggplot <- function(Y, main = "", fill.legend = "", midpoint = NULL, 
                          limits = NULL) {
  require(ggplot2)
  require(maps)
  if (is.null(midpoint)) {
    midpoint <- 1.5
  }
  georgia <- map("county", "georgia", fill = TRUE, col = "transparent",
                 plot = FALSE)
  subregion <- sapply(strsplit(georgia$names, ","), function(x) x[2])
  county_map <- map_data(map = "county", region = "georgia")
  
  basis <- data.frame(Y, subregion)
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