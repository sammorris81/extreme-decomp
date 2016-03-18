rm(list=ls())

source(file = "./package_load.R", chdir = T)
# Number of bases: 4, 9, 16, 25 (using to match the 2d B-splines)
process <- "ebf" # ebf: empirical basis functions, gsk: gaussian kernels
margin  <- "2bs" # ebf: empirical basis functions, 2bs: B Splines in 2d
L       <- 9     # number of knots for ebf or df^2 for 2bs
# if margin == 2bs and L = 4, using 2nd order spatial trend
cv      <- 2     # which cross-validation set to use
results.file <- paste("./cv-results/", process, "-", margin, "-", L, "-", cv, 
                      ".RData", sep = "")
table.file   <- paste("./cv-tables/", process, "-", margin, "-", L, "-", cv, 
                      ".txt", sep = "")

# fit the model and get predictions
source(file = "./fitmodel.R"