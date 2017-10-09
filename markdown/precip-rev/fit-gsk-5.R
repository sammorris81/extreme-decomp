L <- 5                # number of knots to use for the basis functions
cvs <- seq_len(5)     # cross-validation sets
process <- "gsk"      # ebf: empirical basis functions, gsk: gaussian kernels
margin  <- "gsk"      # ebf: empirical basis functions, gsk: gaussian kernels
for (cv in cvs) {
  # don't lose which setting we're running
  rm(list=setdiff(ls(), c("cv", "cvs", "L", "process", "margin")))
  gc()
  source(file = "./package_load.R", chdir = T)
  time    <- "current"  # current or future
  source(file = "./fitmodel.R")
  
  # don't lose which setting we're running
  rm(list=setdiff(ls(), c("cv", "cvs", "L", "process", "margin")))
  gc()
  source(file = "./package_load.R", chdir = T)
  # fit the model and get predictions
  source(file = "./fitmodel.R")
}
