rm(list = ls())

load(file = "../../code/analysis/fire/georgia_preprocess/fire_data.RData")
source(file = "package_load.R")
load(file = "cv-extcoef.RData")
nfolds <- length(cv.idx)
procs <- c("basis", "kern")  # what process determines the spatial structure
nprocs <- length(procs)
bases   <- c(2, 5, 10, 15, 20)
nbases  <- length(bases)
probs.for.qs <- c(0.95, 0.96, 0.97, 0.98, 0.99, 0.995)  # always check fitmodel
nprobs <- length(probs.for.qs)

files <- list.files(path = "cv-results/")

# load in the data
load(file = "../../code/analysis/fire/georgia_preprocess/fire_data.RData")
Y <- t(Y)

for (i in 1:length(files)) {
  split      <- unlist(strsplit(unlist(strsplit(files[i], "-")), "[.]"))
  table.file <- paste("cv-tables/", split[1], "-", split[2], "-", split[3], 
                      ".txt", sep = "")
  
  results.file <- paste("cv-results/", files[i], sep = "")
  load(file = results.file)
  # files are named by the number of basis functions which skips numbers
  fold    <- as.numeric(split[3])
  Y.tst   <- Y[sort(cv.idx[[fold]])]
  qs.update <- QuantScore(preds = fit$y.pred, probs = probs.for.qs, 
                        validate = Y.tst)
  results[1:length(probs.for.qs)] <- qs.update
  write.table(results, file = table.file)
  
  upload.pre <- "samorris@hpc.stat.ncsu.edu:~/repos-git/extreme-decomp/markdown/"
  upload.pre <- paste(upload.pre, "fire-analysis/cv-tables/", sep = "")
  upload.cmd <- paste("scp ", table.file, " ", upload.pre, sep = "")
  
  system(upload.cmd)
  save(B.est, alpha, fit, cv.idx, results, file = results.file)
}

