rm(list=ls())
library(splines)
library(maps)
library(maptools)
library(fields)
library(evd)
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

if (Sys.info()["nodename"] == "cwl-mth-sam-001" | 
    Sys.info()["nodename"] == "cwl-mth-sam-002") {
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
Y <- t(Y)  # Y is originally nt x ns

# set up the 5 fold cross validation
n.tot <- nrow(Y) * ncol(Y)
set.seed(28)  #cv
nfolds <- 5
cv.idx <- get.cv.test(n = n.tot, nfolds = nfolds)

# loop over the list for cross-validation and get the basis functions
# before getting started, find the upper quantile limit

fold <- 1; i <- 1; j <- 2
# check qlims - recording the min and max for all qlims on each fold
qlim.min.range <- matrix(0, nrow = nfolds, ncol = 2)
qlim.max.range <- matrix(0, nrow = nfolds, ncol = 2)
for (fold in 1:nfolds) {
  Y.tst <- Y
  Y.tst[cv.idx[[fold]]] <- NA
  
  # build ec matrix: ns x ns
  ec <- get.pw.ec(Y = Y.tst, qlim = c(0.95, 1), verbose = TRUE, update = 50)
  qlim.min.range[fold, ] <- range(ec$qlims[, 1])
  qlim.max.range[fold, ] <- range(ec$qlims[, 2])
  
  # run smoother
  
  cat("finished fold:", fold, "\n")
}



