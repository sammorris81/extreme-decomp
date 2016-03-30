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

################################################################################
#### Load in the data ##########################################################
################################################################################
# Doing this first because want to do a little cleaning up of non necessary 
# items in the RData file.                                                  
load(file = "precip.RData")
Y <- Yvec
# only keep s and Yvec
rm(Yvec, Ymat, s1, s2)

################################################################################
#### Load functions and other misc setup #######################################
################################################################################
source(file = "../../../usefulR/usefulfunctions.R", chdir = TRUE)
source(file = "../../code/R/auxfunctions.R", chdir = TRUE)
source(file = "../../code/R/PCAX.R", chdir = TRUE)
source(file = "../../code/R/mcmc.R")

# we don't need to upload the beowulf if it's running on beowulf
if (Sys.info()["nodename"] == "cwl-mth-sam-001" | 
    Sys.info()["nodename"] == "cwl-mth-sam-002") {
  setMKLthreads(1)
  do.upload <- TRUE
} else if (Sys.info()["sysname"] == "Darwin") {
  do.upload <- TRUE
} else {
  do.upload <- FALSE
}

################################################################################
#### Preprocess locations and data and setup cross-validation ##################
################################################################################
# get distance matrix
d <- rdist(s)
diag(d) <- 0
ns <- nrow(s)
nt <- ncol(Y)

# set up the 5 fold cross validation
n.tot <- nrow(Y) * ncol(Y)
set.seed(28)  # cv
nfolds <- 5  # picking 5 because data are max-stable and 64 years of data

# stratifying the selection so each training set location loses only 20% of 
# observations
cv.idx <- get.cv.test.strat(data = Y, nfolds = nfolds, idx = 1)

# loop over the list for cross-validation and get the basis functions
# before getting started, find the upper quantile limit

# check qlims - recording the min and max for all qlims on each fold
qlim.min.range <- matrix(0, nrow = nfolds, ncol = 2)
qlim.max.range <- matrix(0, nrow = nfolds, ncol = 2)

ec.hat <- vector(mode = "list", length = nfolds)
for (fold in 1:nfolds) {
  Y.tst <- Y
  Y.tst[cv.idx[[fold]]] <- NA
  
  # build ec matrix: ns x ns
  ec <- get.pw.ec(Y = Y.tst, qlim = c(0, 1), verbose = TRUE, update = 50)
  ec.hat[[fold]] <- ec$ec
  qlim.min.range[fold, ] <- range(ec$qlims[, 1])
  qlim.max.range[fold, ] <- range(ec$qlims[, 2])
  
  cat("finished fold:", fold, "\n")
}

save(cv.idx, ec.hat, file = "cv-extcoef.RData")

library(ggplot2)
library(gridExtra)
# just want to see if it looks as weird when we run all the data
ec <- get.pw.ec(Y = Y, qlim = c(0.90, 1), verbose = TRUE, update = 50)$ec
p.1 <- map.ga.ggplot(Y = ec[, 4], 
                     main = paste("Extremal Coefficients full data"),
                     fill.legend = "EC")

p.2 <- map.ga.ggplot(Y = ec.hat[[1]][, 4], 
                     main = paste("Extremal Coefficients cross validation"),
                     fill.legend = "EC")

p.3 <- map.ga.ggplot(Y = ec.hat[[2]][, 4], 
                     main = paste("Extremal Coefficients cross validation"),
                     fill.legend = "EC")

p.4 <- map.ga.ggplot(Y = ec.hat[[4]][, 4], 
                     main = paste("Extremal Coefficients cross validation"),
                     fill.legend = "EC")

# grid.arrange(p.1, p.2, ncol = 2, widths = c(1.5, 1.5), 
#              top = "EC comparison for CV")

grid.arrange(p.1, p.2, p.3, p.4, ncol = 2, widths = c(1.5, 1.5), 
             top = "EC comparison for CV")

p.1 <- map.ga.ggplot(Y = ec[, 10], 
                     main = paste("Extremal Coefficients full data"),
                     fill.legend = "EC")

p.2 <- map.ga.ggplot(Y = ec.hat[[1]][, 10], 
                     main = paste("Extremal Coefficients cross validation"),
                     fill.legend = "EC")

p.3 <- map.ga.ggplot(Y = ec.hat[[2]][, 10], 
                     main = paste("Extremal Coefficients cross validation"),
                     fill.legend = "EC")

p.4 <- map.ga.ggplot(Y = ec.hat[[4]][, 10], 
                     main = paste("Extremal Coefficients cross validation"),
                     fill.legend = "EC")

grid.arrange(p.1, p.2, p.3, p.4, ncol = 2, widths = c(1.5, 1.5), 
             top = "EC comparison for CV")

d <- rdist(cents)
diag(d) <- 0
n <- nrow(cents)

# standardize the locations
s <- cents
s[, 1] <- (s[, 1] - min(s[, 1])) / diff(range(s[, 1]))
s[, 2] <- (s[, 2] - min(s[, 2])) / diff(range(s[, 2]))

# get a reasonable bandwidth for the kernel smoother 
d <- rdist(s)
diag(d) <- 0
ksmooth.bw <- quantile(d[upper.tri(d)], probs = 0.05)

cat("Start basis function estimation \n")
# basis function estimates using only the training data
L <- 10
out       <- get.factors.EC(ec, L = L, s = s, bw = ksmooth.bw)
B.est     <- out$est
ec.smooth <- out$EC.smooth

plots <- vector(mode = "list", length = L)
for (i in 1:L) {
  title <- paste("basis function", i)
  plots[[i]] <- map.ga.ggplot(Y = B.est[, i], main = title,
                              midpoint = median(B.est[, i]))
}
multiplot(plotlist = plots, cols = 4)
