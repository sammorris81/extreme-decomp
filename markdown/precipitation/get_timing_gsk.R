rm(list = ls())
source(file = "./package_load.R", chdir = T)

Ls <- seq(5, 40, by = 5)
timing <- matrix(NA, 16, 5)
rownames(timing) <- paste(c(rep("cur", 8), rep("fut", 8)), rep(Ls, 2), sep = "-")
for (L in Ls) {
  for (cv in 1:5) {

    # Number of bases: 5, 10, 15, 20
    process <- "gsk"      # ebf: empirical basis functions, gsk: gaussian kernels
    margin  <- "gsk"      # ebf: empirical basis functions, gsk: gaussian kernels
    time    <- "current"  # current or future
    # L       <- 10         # number of knots to use for the basis functions
    # cv      <- 1          # which cross-validation set to use

    loc.fun <- scale.fun <- ~ time + elev # + B1 + B2 + B3 + B4 + B5 + 0

    # fit the model and get predictions
    source(file = "./fitmodel_timing.R")

    # Number of bases: 5, 10, 15, 20
    process <- "gsk"      # ebf: empirical basis functions, gsk: gaussian kernels
    margin  <- "gsk"      # ebf: empirical basis functions, gsk: gaussian kernels
    time    <- "future"   # current or future
    # L       <- 10         # number of knots to use for the basis functions
    # cv      <- 1          # which cross-validation set to use

    loc.fun <- scale.fun <- ~ time + elev # + B1 + B2 + B3 + B4 + B5 + 0

    # fit the model and get predictions
    source(file = "./fitmodel_timing.R")
  }}
