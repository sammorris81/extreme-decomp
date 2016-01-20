rm(list=ls())

library(splines)
library(maps)
library(maptools)
library(fields)
library(ggplot2)
# library(gridExtra)  # not available on hpc
# library(rapport)    # not available on hpc
library(Rcpp)
source(file = "../../../usefulR/usefulfunctions.R", chdir = TRUE)
source(file = "../../code/analysis/fire/adj.R", chdir = TRUE)
source(file = "../../code/R/auxfunctions.R", chdir = TRUE)
source(file = "../../code/R/PCAX.R", chdir = TRUE)
source(file = "../../code/R/mcmc.R")
# load(file = "../code/analysis/fire/gaCntyFires.RData")
load(file = "../../code/analysis/fire/georgia_preprocess/fire_data.RData")

# 

# Number of bases: we can do 1 basis function, but it's not a particularly 
# interesting case because the basis function is 1 at all locations due to the 
# fact that they need to add up to 1 across all basis functions at each site.
# opting for 2, 5, 10, 15, and then looking a max stable method with fixed 
# alpha.
L <- 2
cv <- 1

source(file = "./fitmodel.R")

cat("Finished fitting model \n")