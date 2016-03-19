rm(list=ls())

source(file = "./package_load.R", chdir = T)
# Number of bases: 5, 10, 15, 20
process <- "ebf" # ebf: empirical basis functions, gsk: gaussian kernels
margin  <- "gsk" # ebf: empirical basis functions, gsk: gaussian kernels
L       <- 10    # number of knots to use for the basis functions
cv      <- 5     # which cross-validation set to use
results.file <- paste("./cv-results/", process, "-", margin, "-", L, "-", cv,
                      ".RData", sep = "")
table.file   <- paste("./cv-tables/", process, "-", margin, "-", L, "-", cv,
                      ".txt", sep = "")

# fit the model and get predictions
source(file = "./fitmodel.R")