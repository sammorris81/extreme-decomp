rm(list=ls())
source(file = "./package_load.R", chdir = T)
library(gridExtra)
library(SpatialExtremes)
library(parallel)
options(warn = 0)

L <- 2

#### Try to precalculate the basis functions #########
#### Hoping to save a little time in the analysis ####
load("precip_preprocess.RData")
load("cv-extcoef.RData")
nfolds <- length(cv.idx)

# Empirical basis functions
cat("Starting estimation of empirical basis functions \n")
get.factors.EC.mclapply <- function(x, ec.smooth, alpha.hat, L, s) {
  this.ec.smooth <- ec.smooth[[x]]
  this.alpha.hat <- alpha.hat[x]
  out <- get.factors.EC(EC.smooth = this.ec.smooth,
                        alpha.hat = this.alpha.hat,
                        L = L,
                        s = s)
  return(out)
}

out <- mclapply(1:nfolds,
                get.factors.EC.mclapply,
                ec.smooth = ec.smooth,
                alpha.hat = alpha.hats,
                L = L,
                s = s.scale,
                mc.cores = 4)

ebf.pct <- B.ebf <- vector(mode = "list", length = nfolds)

for (i in 1:nfolds) {
  B.ebf[[i]] <- out[[i]]$est
  ebf.pct[[i]] <- out[[i]]$pct
}

filename <- paste("basis_functions/ebf-", L, ".RData", sep = "")
save(B.ebf, file = filename)
