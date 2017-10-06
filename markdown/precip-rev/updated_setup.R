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

# Save the cv.folds for use when fitting the model.
cv.idx <- vector(mode = "list", length = nfolds)
for (f in seq_len(nfolds)) {
  cv.idx[[f]] <- fold == f
}
save(cv.idx, file = "cv-extcoef.RData")

# Load preprocessed precip data
load("./precip_preprocess.RData")
for (L in nknots) {
  # Brian did the basis functions for EBF separately
  ebf.basis.file <- paste0("./from-bjr/basis_L", L, ".RData")
  load(ebf.basis.file)
  ec.smooth <- B.ebf <- vector(mode = "list", length = nfolds)
  alphas <- rep(0, nfolds)
  for (f in seq_len(nfolds)) {
    B.ebf[[f]] <- OUTPUT[[f]]$est
    ec.smooth[[f]] <- OUTPUT[[f]]$EC.smooth
    alphas[f] <- OUTPUT[[f]]$alpha
  }

  filename <- paste("ebf-", L, ".RData", sep = "")
  save(B.ebf, ec.smooth, alphas, file = filename)

  # Gaussian kernel functions
  set.seed(5687 + L)  # knots + L
  cat("Starting estimation of Gaussian kernels \n")
  knots <- cover.design(s, nd = L)$design
  B.gsk <- vector(mode = "list", length = nfolds)
  alphas <- rep(0, nfolds)
  for (f in 1:nfolds) {
    Y.tst <- Y
    Y.tst[fold == f] <- NA
    ec.hat <- get.chi(Y.tst)
    out   <- get.rho.alpha(EC = ec.hat, s = s, knots = knots,
                           init.rho = 5)
    B.gsk[[f]] <- getW(rho = out$rho, dw2 = out$dw2)
    alphas[f]  <- out$alpha

    cat("  Finished fold ", f, " of ", nfolds, " for gsk. \n", sep = "")
  }

  filename <- paste("gsk-", L, ".RData", sep = "")
  save(B.gsk, alphas, knots, file = filename)

  cat("Finished L = ", L, ".\n", sep = "")
}