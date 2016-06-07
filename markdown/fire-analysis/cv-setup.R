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
ns <- nrow(cents)
Y <- t(Y)  # Y is originally nt x ns

# set up the 10 fold cross validation
set.seed(28)  #cv
nfolds <- 10
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
  ec <- get.pw.ec(Y = Y.tst, qlim = c(0.90, 1), verbose = TRUE, update = 50)
  ec.hat[[fold]] <- ec$ec
  qlim.min.range[fold, ] <- range(ec$qlims[, 1])
  qlim.max.range[fold, ] <- range(ec$qlims[, 2])

  cat("finished fold:", fold, "\n")
}

# standardize the locations
s <- cents
s.scale <- min(diff(range(s[, 1])), diff(range(s[, 2])))
s.min   <- apply(s, 2, min)
s[, 1] <- (s[, 1] - s.min[1]) / s.scale
s[, 2] <- (s[, 2] - s.min[2]) / s.scale

thresh80 <- thresh90 <- thresh95 <- thresh99 <- rep(0, ns)
neighbors <- 5
d <- rdist(s)
diag(d) <- 0

save(cv.idx, ec.hat, file = "cv-extcoef.RData")

#### Try to precalculate the basis functions #########
#### Hoping to save a little time in the analysis ####
load(file = "../../code/analysis/fire/georgia_preprocess/fire_data.RData")
load(file = "../../code/analysis/fire/georgia_preprocess/georgia_map.RData")
load("cv-extcoef.RData")
nfolds <- length(cv.idx)
d <- rdist(cents)
diag(d) <- 0
n <- nrow(cents)

# get candidate knot grid for Gaussian kernel functions
grid.x <- seq(min(cents[, 1]), max(cents[, 1]), length = 100)
grid.y <- seq(min(cents[, 2]), max(cents[, 2]), length = 100)
cents.grid <- as.matrix(expand.grid(grid.x, grid.y))
inGA <- map.where("state", x = cents.grid[, 1], y = cents.grid[, 2])
cents.grid <- cents.grid[inGA == "georgia", ]
cents.grid <- cents.grid[rowSums(is.na(cents.grid)) == 0, ]

# standardize the locations
s <- cents
s.scale <- min(diff(range(s[, 1])), diff(range(s[, 2])))
s.min   <- apply(s, 2, min)
s[, 1] <- (s[, 1] - s.min[1]) / s.scale
s[, 2] <- (s[, 2] - s.min[2]) / s.scale
cents.grid[, 1] <- (cents.grid[, 1] - s.min[1]) / s.scale
cents.grid[, 2] <- (cents.grid[, 2] - s.min[2]) / s.scale

nknots <- c(5, 10, 15, 20, 25, 30, 35, 40)

for (L in nknots[7:8]) {
  # Empirical basis functions
  cat("Starting estimation of empirical basis functions \n")
  alphas <- rep(0, nfolds)
  ec.smooth <- B.ebf <- vector(mode = "list", length = nfolds)
  for (fold in 1:nfolds) {
    out               <- get.factors.EC(ec.hat[[fold]], L = L, s = s)
    B.ebf[[fold]]     <- out$est
    ec.smooth[[fold]] <- out$EC.smooth
    alphas[fold]      <- out$alpha

    cat("  Finished fold ", fold, " of ", nfolds, " for ebf. \n", sep = "")
  }

  filename <- paste("ebf-", L, ".RData", sep = "")
  save(B.ebf, ec.smooth, alphas, file = filename)

  # Gaussian kernel functions
  set.seed(5687 + L)  # knots + L
  cat("Starting estimation of Gaussian kernels \n")
  knots <- cover.design(cents.grid, nd = L)$design
  B.gsk <- vector(mode = "list", length = nfolds)
  for (fold in 1:nfolds) {
    out   <- get.rho.alpha(EC = ec.hat[[fold]], s = s, knots = knots,
                           init.rho = 0.25)
    B.gsk[[fold]] <- getW(rho = out$rho, dw2 = out$dw2)
    alphas[fold]  <- out$alpha

    cat("  Finished fold ", fold, " of ", nfolds, " for gsk. \n", sep = "")
  }

  filename <- paste("gsk-", L, ".RData", sep = "")
  save(B.gsk, alphas, knots, file = filename)

  cat("Finished L = ", L, ".\n", sep = "")
}

#### looks like L = 35 is after things settle down.
# get pairwise extremal coefficients
# build ec matrix: ns x ns
ec <- get.pw.ec.fmado(Y = Y)
ec.hat <- ec$ec
L <- 35

# Empirical basis functions
cat("Starting estimation of empirical basis functions \n")
out       <- get.factors.EC(ec.hat, L = L, s = s)
B.ebf     <- out$est
ec.smooth <- out$EC.smooth
alpha     <- out$alpha

filename <- paste("ebf-", L, "-all.RData", sep = "")
save(B.ebf, ec.smooth, alpha, file = filename)

# Gaussian kernel functions
cat("Starting estimation of Gaussian kernels \n")
set.seed(5687)
knots <- cover.design(cents.grid, nd = L)$design
out   <- get.rho.alpha(EC = ec.hat, s = s, knots = knots)
B.gsk <- getW(rho = out$rho, dw2 = out$dw2)

filename <- paste("gsk-", L, "-all.RData", sep = "")
save(B.gsk, alpha, knots, file = filename)

#### looks like L = 35 is after things settle down.
# get pairwise extremal coefficients
# build ec matrix: ns x ns
ec <- get.pw.ec.fmado(Y = Y)
ec.hat <- ec$ec
L <- 25

# Empirical basis functions
cat("Starting estimation of empirical basis functions \n")
out       <- get.factors.EC(ec.hat, L = L, s = s)
B.ebf     <- out$est
ec.smooth <- out$EC.smooth
alpha     <- out$alpha

filename <- paste("ebf-", L, "-all.RData", sep = "")
save(B.ebf, ec.smooth, alpha, file = filename)

# Gaussian kernel functions
cat("Starting estimation of Gaussian kernels \n")
set.seed(5687)
knots <- cover.design(cents.grid, nd = L)$design
out   <- get.rho.alpha(EC = ec.hat, s = s, knots = knots)
B.gsk <- getW(rho = out$rho, dw2 = out$dw2)

filename <- paste("gsk-", L, "-all.RData", sep = "")
save(B.gsk, alpha, knots, file = filename)

# plot cumsum against basis function
L <- 25
file <- paste("ebf-", L, "-all.RData", sep = "")
load(file)
v <- colSums(B.ebf) / ns
plot(1:L, cumsum(v), ylim = c(0, 1),
     main = paste("Fire analysis (", L, " knots)", sep = ""),
     ylab = "Cumulative relative contribution",
     xlab = "Knot")
plotname <- paste("plots/firev-", L, ".pdf", sep = "")
dev.print(device = pdf, file = plotname,
          width = 6, height = 6)

L <- 35
file <- paste("ebf-", L, "-all.RData", sep = "")
load(file)
v <- colSums(B.ebf) / ns
plot(1:L, cumsum(v), ylim = c(0, 1),
     main = paste("Fire analysis (", L, " knots)", sep = ""),
     ylab = "Cumulative relative contribution",
     xlab = "Knot")
plotname <- paste("plots/firev-", L, ".pdf", sep = "")
dev.print(device = pdf, file = plotname,
          width = 6, height = 6)

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
