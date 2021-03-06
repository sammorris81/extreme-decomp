rm(list=ls())

source(file = "./package_load.R", chdir = T)
# Number of bases: 5, 10, 15, 20
process <- "ebf"      # ebf: empirical basis functions, gsk: gaussian kernels
margin  <- "gsk"      # ebf: empirical basis functions, gsk: gaussian kernels
time    <- "current"  # current or future
L       <- 35         # number of knots to use for the basis functions

loc.fun <- scale.fun <- ~ time + elev # + B1 + B2 + B3 + B4 + B5 + 0

# fit the model and get predictions
source(file = "./fitmodel_nocv.R")
