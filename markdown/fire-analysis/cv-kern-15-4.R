rm(list=ls())

source(file = "./package_load.R", chdir = T)

# Number of bases: we can do 1 basis function, but it's not a particularly 
# interesting case because the basis function is 1 at all locations due to the 
# fact that they need to add up to 1 across all basis functions at each site.
# opting for 2, 5, 10, 15, and then looking a max stable method with fixed 
# alpha.
method <- "kern" # using kern for the results from abba
L      <- 15  # will be using this to get basis functions for covariates
cv     <- 4   # which cross-validation set to use
results.file <- paste("./cv-results/", method, "-", L, "-", cv, ".RData", sep = "")
table.file   <- paste("./cv-tables/", method, "-", L, "-", cv, ".txt", sep = "")


# fit the model and get predictions
source(file = "./fitmodel_kern.R")


rm(list=ls())

source(file = "./package_load.R", chdir = T)

# Number of bases: we can do 1 basis function, but it's not a particularly 
# interesting case because the basis function is 1 at all locations due to the 
# fact that they need to add up to 1 across all basis functions at each site.
# opting for 2, 5, 10, 15, and then looking a max stable method with fixed 
# alpha.
method <- "kern" # using kern for the results from abba
L      <- 15  # will be using this to get basis functions for covariates
cv     <- 9   # which cross-validation set to use
results.file <- paste("./cv-results/", method, "-", L, "-", cv, ".RData", sep = "")
table.file   <- paste("./cv-tables/", method, "-", L, "-", cv, ".txt", sep = "")


# fit the model and get predictions
source(file = "./fitmodel_kern.R")