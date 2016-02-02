rm(list = ls())

load(file = "../../code/analysis/fire/georgia_preprocess/fire_data.RData")
load(file = "cv-extcoef.RData")
nfolds <- length(cv.idx)
basis  <- c(2, 5, 10, 15, "kern")
nbases <- length(basis)
probs.for.qs <- c(0.95, 0.96, 0.97, 0.98, 0.99, 0.995)  # always check fitmodel
nprobs <- length(probs.for.qs)

files <- list.files(path = "cv-tables/")
qs.results <- vector(mode = "list", length = nbases)  # each element is a matrix

for (b in 1:nbases) {
  qs.results[[b]]  <- matrix(NA, nfolds, nprobs)
  colnames(qs.results[[b]])  <- probs.for.qs
  rownames(qs.results[[b]])  <- paste("fold:", 1:nfolds)
}

# timing is a data.frame that contains time, hostname, basis, and fold
timing <- data.frame(timing = double(), hostname = factor(), 
                     basis = factor(), fold = factor())
for (i in 1:length(files)) {
  split     <- unlist(strsplit(unlist(strsplit(files[i], "-")), "[.]"))
  # files are named by the number of basis functions which skips numbers
  setting   <- which(basis == split[1])  
  fold      <- as.numeric(split[2])
  table.set <- read.table(paste("cv-tables/", files[i], sep = ""), 
                          stringsAsFactors = FALSE)
  
  # first extract the timing information from the end of the vector
  timing.tail <- tail(table.set$x, 2)
  timing.row <- data.frame(timing = as.numeric(timing.tail[1]), 
                           host = timing.tail[2],
                           basis = split[1], fold = as.factor(fold))
  timing <- rbind(timing, timing.row)
  qs.results[[setting]][fold, ]  <- as.numeric(table.set$x[1:nprobs])
}

# combine lists into a single matrix that averages qs over all folds for 
# each entry in probs.for.qs
qs.results.mn <- qs.results.se <- matrix(NA, nbases, nprobs)
for (b in 1:nbases) {
  qs.results.mn[b, ]  <- apply(qs.results[[b]][1:5,], 2, mean, na.rm = TRUE)
  qs.results.se[b, ]  <- apply(qs.results[[b]][1:5,], 2, sd, na.rm = TRUE) / sqrt(5)
}

colnames(qs.results.mn) <- colnames(qs.results.se) <- probs.for.qs
rownames(qs.results.mn) <- rownames(qs.results.se) <- basis
round(qs.results.mn, 3)
round(qs.results.se, 3)

for (setting in 1:nsettings) {
  print(length(finished.sets[[setting]]))
}

# how many have finished
colSums(!is.na(bs.results[[1]]))
colSums(!is.na(bs.results[[2]]))
colSums(!is.na(bs.results[[3]]))
colSums(!is.na(bs.results[[4]]))
colSums(!is.na(bs.results[[5]]))
colSums(!is.na(bs.results[[6]]))

round(bs.results.combined[, 1] / bs.results.combined[, 3], 4)
round(bs.results.combined[, 2] / bs.results.combined[, 3], 4)
round(auc.results.combined[, 1], 4)
round(auc.results.combined[, 2], 4)
round(auc.results.combined[, 3], 4)

# Check for differences
# First do Friedman test (one-way repeated measures)
#   friedman.test(y ~ trt | block, data)
# Then follow up with the Wilcoxon, Nemenyi, McDonald-Thompson test
# pWNMT(x, b, trt, method, n.mc)
#     x: list of values
#     b: vector of blocks (only needed if x is a vector)
#     trt: vector of treatments
#     method: "Exact", "Monte Carlo" or "Asymptotic"

library(NSM3)
set.seed(6727)  #npar
groups <- as.factor(rep(1:nmethods, each=50))
dataset <- as.factor(rep(1:50, times=nmethods))
results.friedman <- matrix(0, nsettings, 2)
colnames(results.friedman) <- c("bs", "auc")

for (setting in 1:nsettings) {
  scores <- as.vector(bs.results[[setting]])
  combine <- data.frame(scores, groups, dataset)
  results.friedman[setting, 1] <- friedman.test(scores ~ groups | dataset,
                                                data=combine)$p.value
  
  scores <- as.vector(auc.results[[setting]])
  combine <- data.frame(scores, groups, dataset)
  results.friedman[setting, 2] <- friedman.test(scores ~ groups | dataset,
                                                data=combine)$p.value
}

# posthoc is  Wilcoxon, Nemenyi, McDonald-Thompson test
bs.results.wnmt  <- matrix(0, choose(nmethods, 2), nsettings)
auc.results.wnmt <- matrix(0, choose(nmethods, 2), nsettings)
for (setting in 1:nsettings) {
  scores <- as.vector(bs.results[[setting]])
  combine <- data.frame(scores, groups, dataset)
  bs.results.wnmt[, setting] <- pWNMT(x=combine$scores, b=combine$dataset,
                                      trt=combine$groups, n.mc=20000)$p.val
  
  scores <- as.vector(auc.results[[setting]])
  combine <- data.frame(scores, groups, dataset)
  auc.results.wnmt[, setting] <- pWNMT(x=combine$scores, b=combine$dataset,
                                       trt=combine$groups, n.mc=20000)$p.val
  
  print(paste("setting:", setting))
}

save(bs.results, bs.results.combined, bs.results.wnmt,
     auc.results, auc.results.combined, auc.results.wnmt,
     results.friedman,
     file = "results.RData")

# # look at a few iteration plots
# set <- 1
# setting <- 1
# dataset <- paste("sim-results/", setting, "-", set, ".RData", sep = "")
# load(dataset)
# 
# par(mfrow = c(4, 5))
# for (i in 1:4) {
#   plot(log(fit.gev$a[, i, ]), type = "l", 
#        main = paste("log(a[", i, "])", sep = ""))
# }
# plot(fit.gev$beta, type = "l", main = bquote(beta[0]))
# 
# for (i in 8:11) {
#   plot(log(fit.gev$a[, i, ]), type = "l", 
#        main = paste("log(a[", i, "])", sep = ""))
# }
# plot(fit.gev$alpha, type = "l", main = bquote(alpha))
# 
# for (i in 1:4) {
#   plot(fit.gev$b[, i, ], type = "l", main = paste("b[", i, "]", sep = ""))
# }
# plot(fit.gev$rho, type = "l", main = bquote(rho))
# 
# for (i in 6:10) {
#   plot(fit.gev$b[, i, ], type = "l", main = paste("b[", i, "]", sep = ""))
# }