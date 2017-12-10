library(splines)
library(maps)
library(maptools)
library(fields)
library(ggplot2)
library(emulator)
# library(gridExtra)  # not available on hpc
# library(rapport)    # not available on hpc
library(Rcpp)
library(mvtnorm)
library(extRemes)
library(compiler)
enableJIT(3)
source(file = "../../../usefulR/usefulfunctions.R", chdir = TRUE)
source(file = "../../code/R/auxfunctions.R", chdir = TRUE)
source(file = "../../code/R/get_ebf.R", chdir = TRUE)
source(file = "../../code/R/get_rho.R", chdir = TRUE)
source(file = "../../code/R/mcmc.R", chdir = TRUE)
source(file = "../../code/R/updatemodel.R", chdir = TRUE)
source(file = "./MappingFunctions.R")

if (Sys.info()["nodename"] == "cwl-mth-sam-001") {
  openblas.set.num.threads(1)
  do.upload <- TRUE
} else if (Sys.info()["nodename"] == "sammorris-macbookpro.roam.corp.google.com") {
  do.upload <- FALSE
} else if (Sys.info()["sysname"] == "Darwin") {
  do.upload <- TRUE
} else if (Sys.info()["sysname"] == "Windows") {
  do.upload <- FALSE
  setMKLthreads(1)
} else {
  do.upload <- FALSE
  # set number of threads to use
  openblas.set.num.threads(1)
}

options(warn = 2)
