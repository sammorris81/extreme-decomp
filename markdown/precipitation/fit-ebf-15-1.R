rm(list=ls())

source(file = "./package_load.R", chdir = T)
# Number of bases: 5, 10, 15, 20
process <- "ebf"      # ebf: empirical basis functions, gsk: gaussian kernels
margin  <- "gsk"      # ebf: empirical basis functions, gsk: gaussian kernels
time    <- "current"  # current or future
L       <- 15          # number of knots to use for the basis functions
cv      <- 1          # which cross-validation set to use

# fit the model and get predictions
source(file = "./fitmodel.R")