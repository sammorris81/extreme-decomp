rm(list=ls())

source(file = "./package_load.R", chdir = T)
# Number of bases: 5, 10, 15, 20
process <- "ebf"      # ebf: empirical basis functions, gsk: gaussian kernels
margin  <- "gsk"      # ebf: empirical basis functions, gsk: gaussian kernels
time    <- "future"   # current or future
L       <- 10         # number of knots to use for the basis functions

# fit the model and get predictions
source(file = "./fitmodel_nocv.R")
