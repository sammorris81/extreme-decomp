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
Y <- t(Y)  # when loading in, Y is nt x ns

# set up the 5 fold cross validation
n.tot <- nrow(Y) * ncol(Y)
set.seed(28)  #cv
nfolds <- 5
cv.idx <- get.cv.test(n = n.tot, nfolds = nfolds)

# loop over the list for cross-validation and get the basis functions
# before getting started, find the upper quantile limit

for (i in 1:nfolds) {
  Y.tst <- Y
  Y.tst[cv.idx[[1]]] <- NA
  
  # build ec matrix
  chi <- matrix(NA, nrow(Y), nrow(Y))
  for (i in 1:nrow(Y)) {
    for (j in i:nrow(Y)) {
      chi.ij <- chiplot(Y[c(i, j), ], which = 1, ask = FALSE)
      use.id <- which(chi.ij > 0.95)
      chi.ij <- mean(chiplot(Y[, c(i, j)], which = 1, ask = FALSE)$chi[95:100, 2])
      chi[i, j] <- chi[j, i] <- chi.ij
      if (j %% 50 == 0) {
        print(paste("j:", j))
      }
    }
    if (i %% 10 == 0) {
      print(paste("i:", i))
    }
  }
  
  # run smoother
}


Y <- t(Y)
temp <- chiplot(Y[c(1, 2), ], nq = 1000, which = 1, ask = FALSE)
temp <- chiplot(Y[, c(1, 2)], nq = 1000, which = 1, ask = FALSE)


data <- na.omit(Y[, c(1, 2)])
data <- cbind(rank(data[, 1])/(n + 1), rank(data[, 2])/(n + 1))
rowmax <- apply(data, 1, max)
rowmin <- apply(data, 1, min)
eps <- .Machine$double.eps^0.5
qlim2 <- c(min(rowmax) + eps, max(rowmin) - eps)
