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
Y <- Y.all <- t(Y)

# get the Georgia map and coordinates
# from georgia_preprocess in code/analysis/fire
load(file = "../../code/analysis/fire/georgia_preprocess/georgia_map.RData")
d <- rdist(cents)
diag(d) <- 0
n <- nrow(cents)

# standardize the locations
s <- cents
s[, 1] <- (s[, 1] - min(s[, 1])) / diff(range(s[, 1]))
s[, 2] <- (s[, 2] - min(s[, 2])) / diff(range(s[, 2]))

for (i in 1:length(files)) {
  Y <- Y.all
  split      <- unlist(strsplit(unlist(strsplit(files[i], "-")), "[.]"))
  table.file <- paste("cv-tables/", split[1], "-", split[2], "-", split[3], 
                      ".txt", sep = "")
  
  results.file <- paste("cv-results/", files[i], sep = "")
  load(file = results.file)
  # files are named by the number of basis functions which skips numbers
  fold    <- as.numeric(split[3])
  
  # separate out the test and training data
  Y.tst   <- Y[sort(cv.idx[[fold]])]  # make sure that the order matches preds
  Y[cv.idx[[fold]]] <- NA  # remove the testing data
  
  thresh <- rep(0, nrow(Y))
  neighbors <- 5
  d <- rdist(s)
  diag(d) <- 0
  
  # take the 5 closest neighbors when finding the threshold
  for (j in 1:nrow(Y)) {
    these <- order(d[j, ])[2:(neighbors + 1)]  # the closest is always site i
    thresh[j] <- quantile(Y[these, ], probs = 0.95, na.rm = TRUE)
  }
  thresh <- matrix(thresh, nrow(Y), ncol(Y))
  thresh.tst <- thresh[sort(cv.idx[[fold]])]
  
  qs.update <- QuantScore(preds = fit$y.pred, probs = probs.for.qs, 
                          validate = Y.tst)
  bs.update <- BrierScore(preds = fit$y.pred, validate = Y.tst, 
                          thresh = thresh.tst)
  
  results.old <- results
  results <- c(qs.update, bs.update, tail(results.old, 2))
  names(results) <- c(probs.for.qs, "bs", "timing", "system")
  write.table(results, file = table.file)
  
  upload.pre <- "samorris@hpc.stat.ncsu.edu:~/repos-git/extreme-decomp/markdown/"
  upload.pre <- paste(upload.pre, "fire-analysis/cv-tables/", sep = "")
  upload.cmd <- paste("scp ", table.file, " ", upload.pre, sep = "")
  
  system(upload.cmd)
  save(B.est, alpha, fit, cv.idx, results, file = results.file)
  print(paste("Dataset", i, "of", length(files), "finished"))
}



