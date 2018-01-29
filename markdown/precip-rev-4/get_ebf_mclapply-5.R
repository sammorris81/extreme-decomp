rm(list=ls())
source(file = "./package_load.R", chdir = T)
library(gridExtra)
library(SpatialExtremes)
library(parallel)
options(warn = 0)

L <- 5
mc.cores <- 4

source("./get_ebf_mclapply.R")
