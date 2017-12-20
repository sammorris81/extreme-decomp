L <- 15               # number of knots to use for the basis functions
cv <- 4
cvs <- seq_len(5)     # cross-validation sets
process <- "gsk"      # ebf: empirical basis functions, gsk: gaussian kernels
margin  <- "gsk"      # ebf: empirical basis functions, gsk: gaussian kernels
times <- c("current", "future")
# for (cv in cvs) {
  for (time in times) {
    # don't lose which setting we're running
    rm(list=setdiff(ls(), 
                    c("cv", "cvs", "L", "process", "margin", "time", "times")))
    gc()
    source(file = "./package_load.R", chdir = T)
    source(file = "./fitmodel.R")
  }
# }
