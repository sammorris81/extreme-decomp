# Already set in other file:
# process:= what kind of spatial process (ebf, gsk)
# margin := how to construct marginal basis functions
# cv     := which cross-validation testing set to use
# L      := the number of basis functions to use
options(warn = 2)
library(compiler)
library(Rcpp)
enableJIT(3)
load(file = "../../code/analysis/fire/georgia_preprocess/fire_data.RData")
load(file = "./cv-extcoef.RData")
source(file = "../../code/R/auxfunctions.R", chdir = TRUE)

files.res <- list.files(path = "cv-results/")
files.tab <- list.files(path = "cv-tables/")
nfiles <- length(files.res)
for (i in 1:nfiles) {
  split     <- unlist(strsplit(unlist(strsplit(files.res[i], "-")), "[.]"))
  cv <- as.numeric(split[4])
  Y.all <- Y <- t(Y)
  Y.tst <- Y[cv.idx[[cv]]]  # save the testing data to validate
  thresh90.tst <- thresh90[cv.idx[[cv]]]
  thresh95.tst <- thresh95[cv.idx[[cv]]]
  thresh99.tst <- thresh99[cv.idx[[cv]]]
  filename <- paste("cv-results/", files.res[i], sep = "")
  table.file <- paste("cv-tables/", files.tab[i], sep = "")
  load(filename)
  system <- results[10]

  # calculate the scores
  probs.for.qs <- c(0.95, 0.96, 0.97, 0.98, 0.99, 0.995)
  qs.results <- QuantScore(preds = fit$y.pred[5001:10000, ],
                           probs = probs.for.qs, validate = Y.tst)
  bs.results95 <- BrierScore(preds = fit$y.pred[5001:10000, ],
                             validate = Y.tst, thresh = thresh95.tst)
  bs.results99 <- BrierScore(preds = fit$y.pred[5001:10000, ],
                             validate = Y.tst, thresh = thresh99.tst)
  results <- c(qs.results, bs.results95, bs.results99, fit$timing)
  # results <- c(results, Sys.info()["nodename"])
  results <- c(results, system)
  names(results) <- c(probs.for.qs, "bs-95", "bs-99", "timing", "system")

  write.table(results, file = table.file)
  print(paste("Finished with ", i, sep = ""))
}