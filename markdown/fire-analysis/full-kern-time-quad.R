rm(list=ls())

source(file = "./package_load.R", chdir = T)

# Number of bases: we can do 1 basis function, but it's not a particularly 
# interesting case because the basis function is 1 at all locations due to the 
# fact that they need to add up to 1 across all basis functions at each site.
# opting for 2, 5, 10, 15, and then looking a max stable method with fixed 
# alpha.
method <- "kern" # using kern for the results from abba
L      <- 15  # will be using this to get basis functions for covariates
results.file <- paste("./cv-results/", method, "-", L, ".RData", sep = "")

# load in the data
load(file = "../../code/analysis/fire/georgia_preprocess/fire_data.RData")
# transpose Y because preprocessed forest fire data is Y[t, i]
Y  <- t(Y)

# get the Georgia map and coordinates
# from georgia_preprocess in code/analysis/fire
load(file = "../../code/analysis/fire/georgia_preprocess/georgia_map.RData")
d <- rdist(cents)
diag(d) <- 0
n <- nrow(cents)

# standardize the locations
s <- cents
s[, 1] <- (s[, 1] - min(s[, 1])) / diff(range(s[, 1]))
s[, 2] <- (s[, 2] - min(s[, 2])) / diff(range(s[, 2]))

# find the extremal coefficients
ec.hat <- get.pw.ec(Y = Y, qlim = c(0.90, 1), verbose = TRUE, update = 50)$ec

################################################################################
## Estimate the rho and alpha
################################################################################

cat("Start basis function estimation \n")
# basis function estimates using only the training data
out       <- get.factors.EC(ec.hat, L = L, s = s)
B.est     <- out$est
ec.smooth <- out$EC.smooth

cat("Start estimation of rho and alpha \n")
# alpha and rho estimates using only the training data
out       <- get.rho.alpha(EC = ec.hat, s = s, knots = s)
rho       <- out$rho
ec.smooth <- out$EC.smooth
alpha     <- out$alpha
dw2       <- out$dw2
w         <- getW(rho = rho, dw2 = dw2)  # using w as basis functions in MCMC

################################################################################
#### Run the MCMC:
#### Use the basis functions with the MCMC
#### The response is the total acreage burned in a year 
####   Y[i, t] = acres burned in county i and year t 
####   X[i, t, p] = pth covariate for site i in year t
####     Using (1, time, sites, sites * time) where time = (t - nt / 2) / nt
################################################################################

ns <- nrow(Y)
nt <- ncol(Y)
np <- 3 + L * 3  # for a single year (int, t, t^2, B1...BL, t * (B1...BL),
                 #                    t^2 * (B1...BL))

## create covariate matrix for training
X <- array(1, dim = c(ns, nt, np))
for (i in 1:ns) {
  for (t in 1:nt) { 
    time <- (t - nt / 2) / nt
    X[i, t, 2:np] <- c(time, time^2, B.est[i, ], B.est[i, ] * time, 
                       B.est[i, ] * time^2) 
  }
}

################################################################################
## need spatially smoothed threshold - only using training data
################################################################################
thresh <- rep(0, ns)
neighbors <- 5
d <- rdist(s)
diag(d) <- 0

# take the 5 closest neighbors when finding the threshold
for (i in 1:ns) {
  these <- order(d[i, ])[2:(neighbors + 1)]  # the closest is always site i
  thresh[i] <- quantile(Y[these, ], probs = 0.95, na.rm = TRUE)
}
thresh <- matrix(thresh, ns, nt)

################################################################################
## run the MCMC
################################################################################
iters  <- 30000
burn   <- 20000
update <- 1000

# iters <- 200; burn <- 150; update <- 50  # for testing

cat("Start mcmc fit \n")
set.seed(6262)  # mcmc
# note: we use w as our basis functions because the MCMC should be identical
#       once the bandwidth and alpha terms are estimated.
fit <- ReShMCMC(y = Y, X = X, thresh = thresh, B = w, alpha = alpha, 
                iters = iters, burn = burn, update = update, iterplot = FALSE)

cat("Finished fit and predict \n")

save(B.est, thresh, alpha, rho, fit, file = results.file)
