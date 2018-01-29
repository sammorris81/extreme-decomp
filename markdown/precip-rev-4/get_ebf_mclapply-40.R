rm(list=ls())
source(file = "./package_load.R", chdir = T)
library(gridExtra)
library(SpatialExtremes)
library(parallel)
options(warn = 0)

L <- 40
mc.cores <- 5

source("./get_ebf_mclapply.R")
