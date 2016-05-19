Ls <- c(40, 35, 30, 25, 20, 15, 10, 5)

for (L in Ls) {
  rm(list=ls())

  source(file = "./package_load.R", chdir = T)
  # Number of bases: 5, 10, 15, 20
  process <- "gsk"      # ebf: empirical basis functions, gsk: gaussian kernels
  margin  <- "gsk"      # ebf: empirical basis functions, gsk: gaussian kernels
  time    <- "current"  # current or future
  # L       <- 40         # number of knots to use for the basis functions
  cv      <- 4          # which cross-validation set to use

  loc.fun <- scale.fun <- ~ time + elev # + B1 + B2 + B3 + B4 + B5 + 0

  # fit the model and get predictions
  source(file = "./fitmodel.R")

  rm(list=ls())

  source(file = "./package_load.R", chdir = T)
  # Number of bases: 5, 10, 15, 20
  process <- "gsk"      # ebf: empirical basis functions, gsk: gaussian kernels
  margin  <- "gsk"      # ebf: empirical basis functions, gsk: gaussian kernels
  time    <- "future"   # current or future
  # L       <- 40         # number of knots to use for the basis functions
  cv      <- 4          # which cross-validation set to use

  loc.fun <- scale.fun <- ~ time + elev # + B1 + B2 + B3 + B4 + B5 + 0

  # fit the model and get predictions
  source(file = "./fitmodel.R")
}