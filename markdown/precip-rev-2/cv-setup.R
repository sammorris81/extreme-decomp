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
#### Preprocess locations and data and setup cross-validation ##################
################################################################################
# compute the distances
d <- rdist(s.scale)
diag(d) <- 0
ns <- nrow(s.scale)
nt <- ncol(Y)

# set up the 5 fold cross validation
n.tot <- nrow(Y) * ncol(Y)
set.seed(28)  # cv
nfolds <- 5  # picking 5 because data are max-stable and 64 years of data

# stratifying the selection so each training set location loses only 20% of
# observations
cv.idx <- get.cv.test.strat(data = Y, nfolds = nfolds, idx = 1)

# loop over the list for cross-validation and get the basis functions
# before getting started

ec.hat <- ec.smooth <- vector(mode = "list", length = nfolds)
alpha.hats <-rep(0, nfolds)
for (fold in 1:nfolds) {
  Y.tst <- Y
  Y.tst[cv.idx[[fold]]] <- NA

  # build ec matrix: ns x ns
  ec.hat[[fold]] <- get.ec.fmad(Y = Y.tst, s = s.scale)
  ec.smooth[[fold]] <- KsmoothCV(ec.hat[[fold]], s.scale)$EC
  alpha.hats[fold] <- EstimateAlpha(ec.hat = ec.hat[[fold]], d = d, n0 = 50)

  cat("finished fold:", fold, "\n")
}
save(cv.idx, ec.hat, alpha.hats, ec.smooth, file = "cv-extcoef.RData")

#### Try to precalculate the basis functions #########
#### Hoping to save a little time in the analysis ####
load("precip_preprocess.RData")
load("cv-extcoef.RData")
nfolds <- length(cv.idx)
# openblas.set.num.threads(4)
cents.grid     <- s.scale

nknots <- c(5, 10, 15, 20, 25, 30, 35, 40)

for (L in nknots) {
  # Empirical basis functions
  cat("Starting estimation of empirical basis functions \n")
  B.ebf <- B.gsk <- vector(mode = "list", length = nfolds)

  for (fold in 1:nfolds) {
    out <- get.factors.EC(EC.smooth = ec.smooth[[fold]], 
                          alpha.hat = alpha.hats[fold], L = L, s = s.scale,
                          maxit = 2000, n_starts = 3)
    B.ebf[[fold]] <- out$est
    cat("  Finished fold ", fold, " of ", nfolds, " for ebf. \n", sep = "")
  }

  filename <- paste("basis_functions/ebf-", L, ".RData", sep = "")
  save(B.ebf, file = filename)

  # Gaussian kernel functions
  set.seed(5687 + L)  # knots + L
  cat("Starting estimation of Gaussian kernels \n")
  knots <- cover.design(cents.grid, nd = L)$design
  B.gsk <- vector(mode = "list", length = nfolds)
  for (fold in 1:nfolds) {
    out <- get.rho(EC.smooth = ec.smooth[[fold]], alpha.hat = alpha.hats[fold],
                   s = s.scale, knots = knots,
                   init.rho = 0.3)
    B.gsk[[fold]] <- getW(rho = out$rho, dw2 = out$dw2)
    cat("  Finished fold ", fold, " of ", nfolds, " for gsk. \n", sep = "")
  }

  filename <- paste("basis_functions/gsk-", L, ".RData", sep = "")
  save(B.gsk, knots, file = filename)

  cat("Finished L = ", L, ".\n", sep = "")
}