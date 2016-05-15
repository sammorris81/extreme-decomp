rm(list=ls())

source(file = "./package_load.R", chdir = T)
# Number of bases: 5, 10, 15, 20
process <- "gsk" # ebf: empirical basis functions, gsk: gaussian kernels
margin  <- "gsk" # ebf: empirical basis functions, gsk: gaussian kernels
L       <- 35    # number of knots to use for the basis functions
results.file <- paste("./cv-results/", process, "-", margin, "-", L,
                      "-all.RData", sep = "")
table.file   <- paste("./cv-tables/", process, "-", margin, "-", L,
                      "-all.txt", sep = "")

# fit the model and get predictions
source(file = "./fitmodel_nocv.R")
