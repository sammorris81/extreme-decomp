rm(list=ls())
source(file = "./package_load.R", chdir = T)
library(gridExtra)
library(SpatialExtremes)

# setMKLthreads(5)
################################################################################
#### Load in the data ##########################################################
################################################################################
# Doing this first because want to do a little cleaning up of non necessary
# items in the RData file.
load(file = "precip.RData")
Y <- Yvec
# only keep s and Yvec
elev.old <- elev

# find the elevations
elev <- rep(0, nrow(s))
for (i in 1:nrow(s)) {
  row <- which(s1 == s[i, 1])
  col <- which(s2 == s[i, 2])
  elev[i] <- elev.old[row, col]
}
rm(Yvec, Ymat, s1, s2)

################################################################################
#### For computational restraints, only use s[, 1] > -90 #######################
################################################################################
# This is basically the cutoff used by Reich and Shaby (2012)
keep.these <- which(s[, 1] > -90)
s <- s[keep.these, ]
Y <- Y[keep.these, ]
elev <- elev[keep.these]

# scale the sites
s.scale        <- s
s.scale.factor <- min(diff(range(s[, 1])), diff(range(s[, 2])))
s.min          <- apply(s, 2, min)
s.scale[, 1]   <- (s[, 1] - s.min[1]) / s.scale.factor
s.scale[, 2]   <- (s[, 2] - s.min[2]) / s.scale.factor
save(Y, s, s.scale, elev, year, file = "precip_preprocess.RData")

################################################################################
#### Get ec.hat and smoothed ec matrix                        ##################
################################################################################
# compute the distances
d <- rdist(s.scale)
diag(d) <- 0
ns <- nrow(s.scale)
nt <- ncol(Y)

nknots <- 10
ec.hat <- get.ec.fmad(Y = Y, s = s.scale)
ec.smooth <- KsmoothCV(ec.hat, s.scale)$EC
alpha.hats <- EstimateAlpha(ec.hat = ec.hat, d = d, n0 = 50)
save(ec.hat, alpha.hats, ec.smooth, file = "all-extcoef.RData")

#### Try to precalculate the basis functions #########
#### Hoping to save a little time in the analysis ####
load("precip_preprocess.RData")
load("all-extcoef.RData")
# openblas.set.num.threads(4)
cents.grid     <- s.scale

# Empirical basis functions
cat("Starting estimation of empirical basis functions \n")

out <- get.factors.EC(EC.smooth = ec.smooth,
                      alpha.hat = alpha.hats, L = L, s = s.scale,
                      maxit = 2000, n_starts = 3)
B.ebf <- out$est
cat("  Finished ebf. \n", sep = "")

filename <- paste("basis_functions/ebf-", L, "-all.RData", sep = "")
save(B.ebf, file = filename)

# Gaussian kernel functions
set.seed(5687 + L)  # knots + L
cat("Starting estimation of Gaussian kernels \n")
knots <- cover.design(cents.grid, nd = L)$design

out <- get.rho(EC.smooth = ec.smooth, alpha.hat = alpha.hats,
               s = s.scale, knots = knots,
               init.rho = 0.3)
B.gsk <- getW(rho = out$rho, dw2 = out$dw2)
cat("  Finished gsk. \n", sep = "")

filename <- paste("basis_functions/gsk-", L, "-all.RData", sep = "")
save(B.gsk, knots, file = filename)

