# Rcpp functions
if (!exists("dPSCPP")) {
  sourceCpp(file = "llps.cpp")
}

if (!exists("ifelsematCPP")) {
  sourceCpp(file = "ifelse.cpp")
}

################################################################################
# Common data transformations
################################################################################
transform <- list(
  logit = function(x, lower=0, upper=1) {
    x <- (x - lower) / (upper - lower)
    return(log(x / (1 - x)))
  },
  inv.logit = function(x, lower=0, upper=1) {
    p <- exp(x) / (1 + exp(x))
    p <- p * (upper - lower) + lower
    return(p)
  },
  probit = function(x, lower=0, upper=1) {
    x <- (x - lower) / (upper - lower)
    return(qnorm(x))
  },
  inv.probit = function(x, lower=0, upper=1) {
    p <- pnorm(x)
    p <- p * (upper - lower) + lower
    return(p)
  },
  log = function(x) log(x),
  exp = function(x) exp(x),
  copula = function(dens) {
    this.dens <- paste("p", dens, sep = "")
    function(x, ...) qnorm(do.call(this.dens, args = list(x, ...)))
  },
  inv.copula = function(dens) {
    this.dens <- paste("q", dens, sep = "")
    function(x, ...) do.call(this.dens, args = list(pnorm(x), ...))
  }
)


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

theme_clean <- function(base_size = 12) {
  require(grid)
  theme_grey(base_size) %+replace%
    theme(
      axis.title      =   element_blank(),
      axis.text       =   element_blank(),
      panel.background    =   element_blank(),
      panel.grid      =   element_blank(),
      axis.ticks.length   =   unit(0,"cm"),
      axis.ticks.margin   =   unit(0,"cm"),
      panel.margin    =   unit(0,"lines"),
      plot.margin     =   unit(c(0,0,0,0),"lines"),
      complete = TRUE
    )
}

map.ga.ggplot <- function(Y, main = "", fill.legend = "") {
  require(ggplot2)
  require(maps)
  georgia <- map("county", "georgia", fill = TRUE, col = "transparent",
                 plot = FALSE)
  subregion <- sapply(strsplit(georgia$names, ","), function(x) x[2])
  county_map <- map_data(map = "county", region = "georgia")
  
  basis <- data.frame(Y, subregion)
  extcoef_map <- merge(county_map, basis, all.x = TRUE)
  
  # using fill = Y because that's the column of extcoef_map with the actual data
  p <- ggplot(extcoef_map, aes(x = long, y = lat, group = group, fill = Y))
  p <- p + geom_polygon(colour = "grey", aes(fill = Y))
  p <- p + expand_limits(x = extcoef_map$long, y = extcoef_map$lat)
  p <- p + coord_map("polyconic")
  p <- p + labs(title = main, fill = fill.legend)
  p <- p + scale_fill_gradient2(low = "dodgerblue4", high = "firebrick4", 
                                mid = "#ffffff", midpoint = median(Y))
  p <- p + theme_clean()
  return(p)
}