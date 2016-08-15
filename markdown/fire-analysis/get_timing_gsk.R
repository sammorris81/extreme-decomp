rm(list = ls())
source(file = "./package_load.R", chdir = T)

Ls <- seq(5, 40, by = 5)
timing <- matrix(NA, 8, 5)
rownames(timing) <- Ls
for (L in Ls) {
  for (cv in 1:5) {
    # Number of bases: 5, 10, 15, 20
    process <- "gsk" # ebf: empirical basis functions, gsk: gaussian kernels
    margin  <- "gsk" # ebf: empirical basis functions, gsk: gaussian kernels
    # L       <- 5     # number of knots to use for the basis functions
    # cv      <- 1     # which cross-validation set to use
    # results.file <- paste("./cv-results/", process, "-", margin, "-", L, "-", cv,
    #                       ".RData", sep = "")
    table.file   <- paste("./cv-tables/", process, "-timing.txt", sep = "")

    # fit the model and get predictions
    source(file = "./fitmodel_timing.R")
  }
}
