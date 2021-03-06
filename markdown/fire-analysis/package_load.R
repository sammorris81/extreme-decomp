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
source(file = "../../code/analysis/fire/adj.R", chdir = TRUE)
source(file = "../../code/R/auxfunctions.R", chdir = TRUE)
source(file = "../../code/R/PCAX.R", chdir = TRUE)
source(file = "../../code/R/mcmc.R", chdir = TRUE)
source(file = "../../code/R/updatemodel.R", chdir = TRUE)

if (Sys.info()["nodename"] == "cwl-mth-sam-001") {
  openblas.set.num.threads(1)
  do.upload <- TRUE
} else if (Sys.info()["sysname"] == "Darwin") {
  do.upload <- TRUE
} else {
  do.upload <- FALSE
}