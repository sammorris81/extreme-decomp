rm(list=ls())

library(splines)
library(maps)
library(maptools)
library(fields)
library(ggplot2)
library(gridExtra)
library(rapport)
library(Rcpp)
source(file = "../../code/analysis/fire/adj.R", chdir = TRUE)
source(file = "../../code/R/auxfunctions.R", chdir = TRUE)
source(file = "../../code/R/PCAX.R", chdir = TRUE)
source(file = "../../code/R/mcmc.R")
# load(file = "../code/analysis/fire/gaCntyFires.RData")
load(file = "../../code/analysis/fire/georgia_preprocess/fire_data.RData")

# get the Georgia map and coordinates
# from georgia_preprocess in code/analysis/fire
load(file = "../../code/analysis/fire/georgia_preprocess/georgia_map.RData")
d <- rdist(cents)
diag(d) <- 0
n <- nrow(cents)

# chi[i, j] is the chi statistic between sites i and j
load(file = "../../code/analysis/fire/georgia_preprocess/chi.RData")
chi.hat <- ifelse(chi <= 0, 0, chi)
ec.hat  <- 2 - chi.hat
image.plot(1:n, 1:n, chi.hat, main = "estimated chi")
image.plot(1:n, 1:n, ec.hat, main = "estimated EC")

# plot the estimated extremal coefficients for a few counties
these.counties <- sample(x = 159, size = 10, replace = FALSE)
subregion <- sapply(strsplit(georgia$names, ","), function(x) x[2])

p.1 <- map.ga.ggplot(Y = ec.hat[, 1], 
                     main = paste("Extremal Coefficients for", 
                                  subregion[these.counties[1]], "county"),
                     fill.legend = "EC")

p.2 <- map.ga.ggplot(Y = ec.hat[, 2], 
                     main = paste("Extremal Coefficients for", 
                                  subregion[these.counties[2]], "county"),
                     fill.legend = "EC")

################################################################################
#### Sort and smooth:
#### I don't know that we need to sort, but it makes the plot of the ECs look 
#### better
################################################################################
# sort the counties so it makes sense geographically

# sort the sites by longitude
sites.order.long <- order(cents[, 1])
cents.order.long <- cents[sites.order.long, ]
image.plot(1:n, 1:n, rdist(cents.order.long))

# # more complicated sorting - doesn't appear to improve plot at all
# # dade (#41 the northwest county) is the reference
# sites.first <- 41  # Dade
# sites.last  <- 20  # Camden
# sites.order.for <- order(d[sites.first, ])
# sites.order.rev <- order(d[sites.last, ])
# sites.order <- rep(NA, n)
# sites.order[1:71] <- sites.order.for[1:71]
# sites.order[(n-70):n] <- sites.order.rev[71:1]
# 
# idx <- 72
# for (i in 72:n) {
#   if (!(sites.order.for[i] %in% sites.order)) {
#     sites.order[idx] <- sites.order.for[i] 
#     idx <- idx + 1
#   }
# }
# length(unique(sites.order))
# cents.order <- cents[sites.order, ]
# image.plot(1:n, 1:n, rdist(cents.order))

# Do smoothing
L <- 12

out       <- get.factors.EC(ec.hat,L=L,s=cents)
B.est     <- out$est
alphahat  <- out$alpha
ec.smooth <- out$EC.smooth
ec.est    <- make.EC(B.est, alphahat)

print(out$pct)

# Plot some of the results
p.B1 <- map.ga.ggplot(Y = B.est[, 1], main = "Basis function 1",
                      fill.legend = "Basis value")
p.B2 <- map.ga.ggplot(Y = B.est[, 2], main = "Basis function 2",
                      fill.legend = "Basis value")
p.B3 <- map.ga.ggplot(Y = B.est[, 3], main = "Basis function 3",
                      fill.legend = "Basis value")
p.B4 <- map.ga.ggplot(Y = B.est[, 4], main = "Basis function 4",
                      fill.legend = "Basis value")

grid.arrange(p.B1, p.B2, p.B3, p.B4, ncol = 2, widths = c(1.5, 1.5), 
             top = "Basis functions 1 -- 4")

################################################################################
#### Run the MCMC:
#### Use the basis functions with the MCMC
#### The response is the total acreage burned in a year 
####   Y[i, t] = acres burned in county i and year t 
####   X[i, t, p] = pth covariate for site i in year t
####     Using (1, time, B, B * time) where time = (t - nt / 2) / nt
################################################################################

## transpose Y because preprocessed forest fire data is Y[t, i]
Y <- t(Y)
nt <- ncol(Y)
ns <- nrow(Y)
np <- ncol(out$est) * 2 + 2

## create covariate matrix
X <- array(1, dim = c(ns, nt, np))
for (i in 1:ns) {
  for (t in 1:nt) {
    time <- (t - nt / 2) / nt
    X[i, t, 2:np] <- c(time, out$est[i, ], out$est[i, ] * time) 
  }
}

## need spatially smoothed threshold
thresh <- rep(0, ns)

# first standardize the centroid locations
cents.std <- cents
cents.rng <- apply(cents.std, 2, range)
cents.std[, 1] <- (cents.std[, 1] - cents.rng[1, 1]) / 
  (cents.rng[2, 1] - cents.rng[1, 1])
cents.std[, 2] <- (cents.std[, 2] - cents.rng[1, 2]) / 
  (cents.rng[2, 2] - cents.rng[1, 2])

# get the distances between standardized quantiles
d.std <- rdist(cents.std)
diag(d.std) <- 0

# set the radius to q(0.05) of the distances and take q(0.95) across the sites
# within the radius across all the years
rad <- quantile(d.std[upper.tri(d)], probs = 0.05)
for (i in 1:ns) {
  these <- which(d.std[i, ] < rad)
  thresh[i] <- quantile(Y[these, ], probs = 0.95)
}

# plot smoothed quantile
p.3 <- map.ga.ggplot(Y = thresh, 
                     main = paste("Spatially smoothed q(0.95) of acreage burned, ",
                                  "bandwidth = ", round(rad, 3)),
                     fill.legend = "Acreage burned")
p.3

train <- seq(1, 159, 2)
test  <- seq(2, 158, 2)

fit.1 <- ReShMCMC(y = Y[train, ], X = X[train, , ], thresh = thresh[train], 
                B = out$est[train, , drop = FALSE], alpha = out$alpha, 
                iters = 300, burn = 100, update = 10, iterplot = TRUE)

fit.2 <- ReShMCMC(y = Y, X = X, thresh = thresh, 
                B = out$est, alpha = out$alpha, 
                iters = 300, burn = 100, update = 10, iterplot = TRUE)

fit.3 <- ReShMCMC(y = Y[train, ], X = X[train, , 1:12], thresh = thresh[train], 
                  B = out$est[train, 1:5, drop = FALSE], alpha = out$alpha, 
                  iters = 300, burn = 100, update = 10, iterplot = TRUE)

fit.4 <- ReShMCMC(y = Y, X = X[, , 1:12], thresh = thresh, 
                  B = out$est[, 1:5, drop = FALSE], alpha = out$alpha, 
                  iters = 300, burn = 100, update = 10, iterplot = TRUE)

y.pred <- pred.ReShMCMC(mcmcoutput = fit, X.pred = X[test, , ], 
                        B = out$est[test, , drop = FALSE], alpha = out$alpha, 
                        start = 1, end = 200, update = 10)

rm(list=ls())
library(splines)
library(maps)
library(maptools)
library(fields)
library(ggplot2)
# library(gridExtra)  # not available on hpc
# library(rapport)    # not available on hpc
library(Rcpp)
source(file = "../../../usefulR/usefulfunctions.R", chdir = TRUE)
source(file = "../../code/analysis/fire/adj.R", chdir = TRUE)
source(file = "../../code/R/auxfunctions.R", chdir = TRUE)
source(file = "../../code/R/PCAX.R", chdir = TRUE)
source(file = "../../code/R/mcmc.R")
# load(file = "../code/analysis/fire/gaCntyFires.RData")
load(file = "../../code/analysis/fire/georgia_preprocess/fire_data.RData")

if (Sys.info()["nodename"] == "sam-ubuntu") {
  setMKLthreads(1)
  do.upload <- TRUE
} else if (Sys.info()["sysname"] == "Darwin") {
  do.upload <- TRUE
} else {
  do.upload <- FALSE
}

# test for missing
# get the Georgia map and coordinates
# from georgia_preprocess in code/analysis/fire
load(file = "../../code/analysis/fire/georgia_preprocess/georgia_map.RData")
d <- rdist(cents)
diag(d) <- 0
n <- nrow(cents)

Y <- t(Y)

# set up the 5 fold cross validation
n.tot <- nrow(Y) * ncol(Y)
set.seed(28)  #cv
nfolds <- 5
cv.idx <- get.cv.test(n = n.tot, nfolds = nfolds)

Y.tst <- Y
Y.tst[cv.idx[[1]]] <- NA




fit.1 <- ReShMCMC(y = Y, X = X, thresh = thresh, 
                  B = out$est, alpha = out$alpha, 
                  iters = 300, burn = 100, update = 10, iterplot = TRUE)
