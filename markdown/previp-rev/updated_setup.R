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
fold <- matrix(sample(1:nfolds, ns * nt, replace = TRUE), ns, nt)

# Load preprocessed precip data
load("./precip_preprocess.RData")

# Run GSK smoother
ec.hat <- vector(mode = "list", length = nfolds)

#### NEED get.chi FUNCTION FROM BRIAN

for (fold in 1:nfolds) {
  Y.tst <- Y
  Y.tst[fold == f] <- NA

  # # build ec matrix: ns x ns
  # ec <- get.pw.ec.fmado(Y = Y.tst)
  # ec.hat[[fold]] <- ec$ec
  this.ec <- fmadogram(data = t(Y.tst), coord = s, which = "ext")
  this.ec[this.ec >= 2] <- 2
  ec <- matrix(1, nrow(Y), nrow(Y))
  ec[lower.tri(ec)] <- this.ec[, 3]
  ec[upper.tri(ec)] <- t(ec)[upper.tri(ec)]
  ec.hat[[fold]] <- ec

  cat("finished fold:", fold, "\n")
}

for (L in nknots) {
  # Gaussian kernel functions
  set.seed(5687 + L)  # knots + L
  cat("Starting estimation of Gaussian kernels \n")
  knots <- cover.design(s, nd = L)$design
  B.gsk <- vector(mode = "list", length = nfolds)
  for (fold in 1:nfolds) {
    Y.tst <- Y
    Y.tst[fold == f] <- NA
    ec.hat <- get.chi(Y.tst)
    out   <- get.rho.alpha(EC = ec.hat, s = s, knots = knots,
                           init.rho = 3)
    B.gsk[[fold]] <- getW(rho = out$rho, dw2 = out$dw2)
    alphas[fold]  <- out$alpha

    cat("  Finished fold ", fold, " of ", nfolds, " for gsk. \n", sep = "")
  }

  filename <- paste("gsk-", L, ".RData", sep = "")
  save(B.gsk, alphas, knots, file = filename)

  cat("Finished L = ", L, ".\n", sep = "")
}