rm(list=ls())

source(file = "./package_load.R", chdir = T)

# Number of bases: In this method we are combaring the performance of the abba
# function to our proposed method which uses basis functions. So, here
# there aren't a set of basis functions.
L <- "ms" # using ms for the results from abba
cv <- 4   # which cross-validation set to use
results.file <- paste("./cv-results/", L, "-", cv, ".RData", sep = "")
table.file   <- paste("./cv-tables/", L, "-", cv, ".txt", sep = "")

# fit the model and get predictions
source(file = "./fitmodel_ms.R")


rm(list=ls())

source(file = "./package_load.R", chdir = T)

# Number of bases: we can do 1 basis function, but it's not a particularly 
# interesting case because the basis function is 1 at all locations due to the 
# fact that they need to add up to 1 across all basis functions at each site.
# opting for 2, 5, 10, 15, and then looking a max stable method with fixed 
# alpha.
L <- "ms" # number of basis functions
cv <- 9   # which cross-validation set to use
results.file <- paste("./cv-results/", L, "-", cv, ".RData", sep = "")
table.file   <- paste("./cv-tables/", L, "-", cv, ".txt", sep = "")

# fit the model and get predictions
source(file = "./fitmodel.R")