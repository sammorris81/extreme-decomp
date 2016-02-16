rm(list = ls())

load(file = "../../code/analysis/fire/georgia_preprocess/fire_data.RData")
load(file = "cv-extcoef.RData")
nfolds <- length(cv.idx)
procs <- c("basis", "kern")  # what process determines the spatial structure
nprocs <- length(procs)
bases   <- c(2, 5, 10, 15, 20)
nbases  <- length(bases)
probs.for.qs <- c(0.95, 0.96, 0.97, 0.98, 0.99, 0.995)  # always check fitmodel
nprobs <- length(probs.for.qs)

files <- list.files(path = "cv-tables/")
qs.results <- vector(mode = "list", length = nbases)  # each element is a matrix
bs.results <- matrix(NA, nbases * nprocs, nfolds)

for (b in 1:(nbases * nprocs)) {
  qs.results[[b]]  <- matrix(NA, nfolds, nprobs)
  colnames(qs.results[[b]])  <- probs.for.qs
  rownames(qs.results[[b]])  <- paste("fold:", 1:nfolds)
}

# timing is a data.frame that contains time, hostname, basis, and fold
timing <- data.frame(timing = double(), hostname = factor(), 
                     proc = factor(), basis = factor(), fold = factor())
for (i in 1:(length(files) - 1)) {  # last file is timing.txt
  split     <- unlist(strsplit(unlist(strsplit(files[i], "-")), "[.]"))
  # files are named by the number of basis functions which skips numbers
  proc.idx  <- which(procs == split[1])
  basis.idx <- which(bases == split[2])  
  idx       <- (proc.idx - 1) * nbases + basis.idx
  fold      <- as.numeric(split[3])
  table.set <- read.table(paste("cv-tables/", files[i], sep = ""), 
                          stringsAsFactors = FALSE)
  
  # first extract the timing information from the end of the vector
  timing.tail <- tail(table.set$x, 2)
  timing.row <- data.frame(timing = as.numeric(timing.tail[1]), 
                           host = timing.tail[2],
                           proc = split[1], basis = split[2], 
                           fold = as.factor(fold))
  timing <- rbind(timing, timing.row)
  qs.results[[idx]][fold, ]  <- as.numeric(table.set$x[1:nprobs])
  bs.results[idx, fold] <- as.numeric(table.set$x[nprobs + 1])
}

# combine lists into a single matrix that averages qs over all folds for 
# each entry in probs.for.qs
# CHECK to make sure you're only including the folds that you want
qs.results.mn <- qs.results.se <- matrix(NA, nbases * 2, nprobs)
for (p in 1:nprocs) {
  for (b in 1:nbases) {
    this.row <- (p - 1) * nbases + b
    this.qs <- qs.results[[this.row]][1:5, ]
    qs.results.mn[this.row, ]  <- apply(this.qs, 2, mean, 
                                        na.rm = TRUE)
    qs.results.se[this.row, ]  <- apply(this.qs, 2, sd, 
                                        na.rm = TRUE) / sqrt(10)
  }
}

colnames(qs.results.mn) <- colnames(qs.results.se) <- probs.for.qs
these.rownames <- paste(rep(procs, each = nbases), 
                        rep(bases, times = nprocs))
rownames(qs.results.mn) <- rownames(qs.results.se) <- these.rownames
rownames(bs.results) <- these.rownames

bs.results.mn <- apply(bs.results, 1, mean)
bs.results.se <- apply(bs.results, 1, sd) / sqrt(10)

round(qs.results.mn, 3)
round(qs.results.se, 3)

matplot(x = probs.for.qs, y = t(qs.results.mn),
        col = c(rep("dodgerblue4", nbases), 
                rep("firebrick4", nbases)), 
        bg = c(rep("dodgerblue1", nbases),
               rep("firebrick1", nbases)),
        type = "b", pch = c(rep(21, nbases), rep(22, nbases)),
        lty = c(rep(1, nbases), rep(2, nbases)),
        cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex = 1.5, lwd=1.5,
        main = "Average quantile score at selected quantiles",
        ylab = "Mean quantile score", xlab = "Quantile")
legend("topright", legend = c("Basis functions", "Kernel"), 
       lty = c(1, 2), pch = c(21, 22), 
       col = c("dodgerblue4", "firebrick4"), 
       pt.bg = c("dodgerblue1", "firebrick1"),
       title = "Spatial Process", 
       lwd = 1.5, cex=1.5)

# look at timing - try to keep like computers together (in minutes)
timing <- read.table(file = "./cv-tables/timing.txt")
timing.mn <- apply(timing, 1, mean)

dev.new(width = 4, height = 3)
par(mfrow=c(1, 1), mar=c(5.1, 5.1, 4.1, 2.1))
ylim = range(timing.mn)
ylim[2] <- ylim[2] + 100
plot(bases, timing.mn[1:5], type = "b", 
     ylim = ylim, 
     xaxt = "n", xlab = "Number of basis functions", 
     ylab = "Timing", 
     # main = "Timing comparison of 100 iterations: Basis functions vs kernel", 
     cex.lab=1.5, cex.axis=1.5, cex.main=1.5, lwd=1.5)
axis(1, at = bases, labels = TRUE)
lines(bases, timing.mn[6:10], lty = 2, type = "b", lwd = 1.5)
legend("topleft", legend = c("Basis functions", "Kernel"), 
       lty = c(1, 2), lwd = 1.5, pch = 1,
       title = "Spatial Process", 
       cex=1.5)
dev.print(device = pdf, file = "plots/timing.pdf")
dev.off()


rm(list=ls())
library(fields)
# look at results
# load in the data
load(file = "../../code/analysis/fire/georgia_preprocess/fire_data.RData")
Y <- Y.all <- t(Y)
ns <- nrow(Y)
nt <- ncol(Y)

fold <- 1
load(paste("./cv-results/basis-10-", fold, ".RData", sep = ""))
fit.basis <- fit
load(paste("./cv-results/kern-10-", fold, ".RData", sep = ""))
fit.kern <- fit

Y.tst <- Y[sort(cv.idx[[fold]])]
Y[cv.idx[[fold]]] <- NA

################################################################################
## need spatially smoothed threshold - only using training data
################################################################################
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

thresh <- rep(0, ns)
neighbors <- 5
d <- rdist(s)
diag(d) <- 0

# take the 5 closest neighbors when finding the threshold
for (i in 1:ns) {
  these <- order(d[i, ])[2:(neighbors + 1)]  # the closest is always site i
  thresh[i] <- quantile(Y[these, ], probs = 0.95, na.rm = TRUE)
}
thresh <- matrix(thresh, ns, nt)
thresh.tst <- thresh[sort(cv.idx[[fold]])]

this.Y <- Y.tst[1]
this.thresh <- thresh.tst[1]
quantile(fit.basis$y.pred[, 1], probs = c(0.025, 0.975))
quantile(fit.kern$y.pred[, 1], probs = c(0.025, 0.975))

this.Y <- Y[sort(cv.idx[[1]])[2]]
quantile(fit.basis$y.pred[, 2], probs = c(0.025, 0.975))
quantile(fit.kern$y.pred[, 2], probs = c(0.025, 0.975))

this.Y <- Y[sort(cv.idx[[1]])[711]]
quantile(fit.basis$y.pred[, 711], probs = c(0.025, 0.975))
quantile(fit.kern$y.pred[, 711], probs = c(0.025, 0.975))

dim(fit.basis$beta2)

plot(fit.basis$beta1[, 1], type = "l")
plot(fit.kern$beta1[, 1], type = "l")
plot(fit.basis$beta1[, 6], type = "l")
plot(fit.kern$beta1[, 6], type = "l")
plot(fit.basis$beta2[, 1], type = "l")
plot(fit.kern$beta2[, 1], type = "l")
plot(fit.basis$beta2[, 6], type = "l")
plot(fit.kern$beta2[, 6], type = "l")

hist(fit.basis$y.pred[, 1])
hist(fit.kern$y.pred[, 1])

