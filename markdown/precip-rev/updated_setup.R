# This file takes the revised basis functions and folds from Brian and
# recomputes the Gaussian kernel functions with the new folds.

rm(list = ls())
source("./package_load.R")
load("./precip_preprocess.RData")
ns   <- nrow(Y)
nt   <- ncol(Y)

# Get folds
nknots <- c(5, 10, 15, 20, 25, 30, 35, 40)
nfolds <- 5  # picking 5 because data are max-stable and 64 years of data

set.seed(0820)
# stratifying the selection so each training set location loses only 20% of
# observations
cv.idx <- get.cv.test.strat(data = Y, nfolds = nfolds, idx = 1)
save(cv.idx, file = "cv-extcoef.RData")

#### Try to precalculate the basis functions #########
#### Hoping to save a little time in the analysis ####
load("precip_preprocess.RData")
load("cv-extcoef.RData")
nfolds <- length(cv.idx)
# openblas.set.num.threads(4)
s.scale        <- s
s.scale.factor <- min(diff(range(s[, 1])), diff(range(s[, 2])))
s.min          <- apply(s, 2, min)
s.scale[, 1]   <- (s[, 1] - s.min[1]) / s.scale.factor
s.scale[, 2]   <- (s[, 2] - s.min[2]) / s.scale.factor
cents.grid     <- s.scale

# Load preprocessed precip data
load("./precip_preprocess.RData")
for (L in nknots) {
  # Storage
  ec.smooth <- B.ebf <- B.gsk <- vector(mode = "list", length = nfolds)
  alphas.ebf <- alphas.gsk <- rep(0, nfolds)
  set.seed(5687 + L)  # knots + L
  knots <- cover.design(s, nd = L)$design
  
  for (fold in seq_len(nfolds)) {
    Y.tst <- Y
    Y.tst[cv.idx[[fold]]] <- NA
    ec.hat <- get.chi(Y.tst)
    
    # Get the basis function using the EBF method
    out.ebf <- get.factors.EC(EChat = ec.hat, L = L, s = s.scale, 
                              verbose = TRUE)
    B.ebf[[f]] <- out.ebf$est
    ec.smooth[[f]] <- out.ebf$EC.smooth
    alphas.ebf[f] <- out.ebf$alpha
    
    # Get the basis functions using the GSK method
    out   <- get.rho.alpha(EC = ec.hat, s = s.scale, knots = knots,
                           init.rho = 0.3)
    B.gsk[[f]] <- getW(rho = out$rho, dw2 = out$dw2)
    alphas.gsk[f]  <- out$alpha

    cat("  Finished fold ", f, " of ", nfolds, " for gsk. \n", sep = "")
  }

  filename <- paste("ebf-", L, ".RData", sep = "")
  save(B.ebf, ec.smooth, alphas.ebf, file = filename)
  
  filename <- paste("gsk-", L, ".RData", sep = "")
  save(B.gsk, alphas.gsk, knots, file = filename)

  cat("Finished L = ", L, ".\n", sep = "")
}