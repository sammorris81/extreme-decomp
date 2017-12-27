rm(list=ls())
source(file = "./package_load.R", chdir = T)
library(gridExtra)
library(SpatialExtremes)

#### Try to precalculate the basis functions #########
#### Hoping to save a little time in the analysis ####
load("precip_preprocess.RData")
load("cv-extcoef.RData")
nfolds <- length(cv.idx)
# openblas.set.num.threads(4)
cents.grid     <- s.scale

L <- 3

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
  print(out$rho)
  cat("  Finished fold ", fold, " of ", nfolds, " for gsk. \n", sep = "")
}

filename <- paste("basis_functions/gsk-", L, ".RData", sep = "")
save(B.gsk, knots, file = filename)

cat("Finished L = ", L, ".\n", sep = "")
