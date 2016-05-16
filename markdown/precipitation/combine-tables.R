rm(list = ls())

load(file = "cv-extcoef.RData")
nfolds <- length(cv.idx)
procs  <- c("ebf", "gsk")  # what process determines the spatial structure
nprocs <- length(procs)
times  <- c("current", "future")  # is this current or future data
ntimes <- length(times)
margs  <- "gsk"
nmargs <- length(margs)
bases  <- c(5, 10, 15, 20, 25, 30, 35, 40)
nbases <- length(bases)
probs.for.qs <- c(0.95, 0.96, 0.97, 0.98, 0.99, 0.995)  # always check fitmodel
probs.for.bs <- c(0.95, 0.99)
nprobs.qs <- length(probs.for.qs)
nprobs.bs <- length(probs.for.bs)

files <- list.files(path = "cv-tables/")
# each element of these lists is a matrix - including an extra for gsk-gsk-all
qs.results <- vector(mode = "list", length = nbases * nprocs * nmargs)
bs.results <- vector(mode = "list", length = nbases * nprocs * nmargs)

for (b in 1:(ntimes * nbases * nmargs * nprocs + 1)) {
  qs.results[[b]] <- matrix(NA, nfolds, nprobs.qs)
  bs.results[[b]] <- matrix(NA, nfolds, nprobs.bs)
  colnames(qs.results[[b]]) <- probs.for.qs
  colnames(bs.results[[b]]) <- probs.for.bs
  rownames(qs.results[[b]]) <- paste("fold:", 1:nfolds)
  rownames(bs.results[[b]]) <- paste("fold:", 1:nfolds)
}

# timing is a data.frame that contains time, hostname, basis, and fold
timing <- data.frame(timing = double(), hostname = factor(),
                     proc = factor(), time = factor(),
                     basis = factor(), fold = factor())
rownames <- rep(0, ntimes * nbases * nmargs * nprocs)
for (i in 1:(length(files))) {  # last file is timing.txt
  split     <- unlist(strsplit(unlist(strsplit(files[i], "-")), "[.]"))
  # files are named by the number of basis functions which skips numbers
  proc.idx  <- which(procs == split[1])
  time.idx  <- which(times == split[2])
  basis.idx <- which(bases == split[3])

  # idx: 1 - 5: ebf spatial, ebf marginal
  # idx: 6 - 10: ebf spatial, gsk marginal
  # idx: 11 - 15: gsk spatial, ebf marginal
  # idx: 16 - 20: gsk spatial, gsk marginal
  # idx: 21: gsk spatial, gsk marginal - knots at all counties

  if (as.numeric(split[3]) == 159) {
    idx <- length(rownames)
  } else {
    idx       <- (proc.idx - 1) * (nbases * ntimes) +
      (time.idx - 1) * nbases + basis.idx
  }

  if (time.idx == 1) {
    rownames[idx] <- paste(split[1], "-cur-", split[3], sep = "")
  } else {
    rownames[idx] <- paste(split[1], "-fut-", split[3], sep = "")
  }

  fold      <- as.numeric(split[4])
  table.set <- read.table(paste("cv-tables/", files[i], sep = ""),
                          stringsAsFactors = FALSE)

  # first extract the timing information from the end of the vector
  timing.tail <- tail(table.set$x, 2)
  timing.row <- data.frame(timing = as.numeric(timing.tail[1]),
                           host = timing.tail[2],
                           proc = split[1], time = split[2],
                           basis = split[3], fold = as.factor(fold))
  timing <- rbind(timing, timing.row)
  qs.results[[idx]][fold, ] <- as.numeric(table.set$x[1:nprobs.qs])
  bs.start <- nprobs.qs + 1
  bs.end   <- bs.start + 1
  bs.results[[idx]][fold, ] <- as.numeric(table.set$x[bs.start:bs.end])
}

# combine lists into a single matrix that averages qs over all folds for
# each entry in probs.for.qs
# CHECK to make sure you're only including the folds that you want
qs.results.mn <- qs.results.se <- matrix(NA, nbases * ntimes * 2, nprobs.qs)
bs.results.mn <- bs.results.se <- matrix(NA, nbases * ntimes * 2, nprobs.bs)
# idx: 1 - 5: ebf spatial, ebf marginal
# idx: 6 - 10: ebf spatial, gsk marginal
# idx: 11 - 15: gsk spatial, ebf marginal
# idx: 16 - 20: gsk spatial, gsk marginal
# idx: 21: gsk spatial, gsk marginal, knots at all locations
for (p in 1:nprocs) {
  for (t in 1:ntimes) {
    for (b in 1:nbases) {
      this.row <- (p - 1) * (nbases * ntimes) +
        (t - 1) * nbases + b
      this.qs <- qs.results[[this.row]]
      this.bs <- bs.results[[this.row]]
      qs.results.mn[this.row, ] <- apply(this.qs, 2, mean,
                                         na.rm = TRUE)
      bs.results.mn[this.row, ] <- apply(this.bs, 2, mean,
                                         na.rm = TRUE)
      qs.results.se[this.row, ] <- apply(this.qs, 2, sd,
                                         na.rm = TRUE) / sqrt(5)
      bs.results.se[this.row, ] <- apply(this.bs, 2, sd,
                                         na.rm = TRUE) / sqrt(5)
    }
  }
}

rownames(bs.results.mn) <- rownames
rownames(qs.results.mn) <- rownames

these.ebf.cur <- 1:8
these.ebf.fut <- 9:16
these.gsk.cur <- 17:24
these.gsk.fut <- 25:32

# Brier scores
quartz(width = 8, height = 8)
par(mfrow = c(2, 2))
plot(seq_along(these.ebf.cur), bs.results.mn[these.ebf.cur, 1], type = "l",
     main = "Current: Brier score for q(0.95)",
     ylab = "Brier score", xlab = "Knots",
     ylim = range(bs.results.mn[c(these.ebf.cur, these.gsk.cur), 1]),
     xaxt = "n")
lines(seq_along(these.gsk.cur), bs.results.mn[these.gsk.cur, 1], lty = 2)
axis(1, at = 1:8, labels = seq(5, 40, by = 5))
legend("topright", lty = c(1, 2), legend = c("EBF   ", "GSK   "))

plot(seq_along(these.ebf.fut), bs.results.mn[these.ebf.fut, 1], type = "l",
     main = "Future: Brier score for q(0.95)",
     ylab = "Brier score", xlab = "Knots",
     ylim = range(bs.results.mn[c(these.ebf.fut, these.gsk.fut), 1]),
     xaxt = "n")
lines(seq_along(these.gsk.fut), bs.results.mn[these.gsk.fut, 1], lty = 2)
axis(1, at = 1:8, labels = seq(5, 40, by = 5))
legend("topright", lty = c(1, 2), legend = c("EBF   ", "GSK   "))

plot(seq_along(these.ebf.cur), bs.results.mn[these.ebf.cur, 2], type = "l",
     main = "Current: Brier score for q(0.99)",
     ylab = "Brier score", xlab = "Knots",
     ylim = range(bs.results.mn[c(these.ebf.cur, these.gsk.cur), 2]),
     xaxt = "n")
lines(seq_along(these.gsk.cur), bs.results.mn[these.gsk.cur, 2], lty = 2)
axis(1, at = 1:8, labels = seq(5, 40, by = 5))
legend("topright", lty = c(1, 2), legend = c("EBF   ", "GSK   "))

plot(seq_along(these.ebf.fut), bs.results.mn[these.ebf.fut, 2], type = "l",
     main = "Future: Brier score for q(0.99)",
     ylab = "Brier score", xlab = "Knots",
     ylim = range(bs.results.mn[c(these.ebf.fut, these.gsk.fut), 2]),
     xaxt = "n")
lines(seq_along(these.gsk.fut), bs.results.mn[these.gsk.fut, 2], lty = 2)
axis(1, at = 1:8, labels = seq(5, 40, by = 5))
legend("topright", lty = c(1, 2), legend = c("EBF   ", "GSK   "))
dev.print(device = pdf, file = "./plots/precip-bs.pdf")
dev.off()

#### Quantile scores
quartz(width = 8, height = 8)
par(mfrow = c(2, 2))
plot(seq_along(these.ebf.cur), qs.results.mn[these.ebf.cur, 1], type = "l",
     main = "Current: Quantile score for q(0.95)",
     ylab = "Quantile score", xlab = "Knots",
     ylim = range(qs.results.mn[c(these.ebf.cur, these.gsk.cur), 1]),
     xaxt = "n")
lines(seq_along(these.gsk.cur), qs.results.mn[these.gsk.cur, 1], lty = 2)
axis(1, at = 1:8, labels = seq(5, 40, by = 5))
legend("topright", lty = c(1, 2), legend = c("EBF   ", "GSK   "))

plot(seq_along(these.ebf.fut), qs.results.mn[these.ebf.fut, 1], type = "l",
     main = "Future: Quantile score for q(0.95)",
     ylab = "Quantile score", xlab = "Knots",
     ylim = range(qs.results.mn[c(these.ebf.fut, these.gsk.fut), 1]),
     xaxt = "n")
lines(seq_along(these.gsk.fut), qs.results.mn[these.gsk.fut, 1], lty = 2)
axis(1, at = 1:8, labels = seq(5, 40, by = 5))
legend("topright", lty = c(1, 2), legend = c("EBF   ", "GSK   "))

plot(seq_along(these.ebf.cur), qs.results.mn[these.ebf.cur, 5], type = "l",
     main = "Current: Quantile score for q(0.99)",
     ylab = "Quantile score", xlab = "Knots",
     ylim = range(qs.results.mn[c(these.ebf.cur, these.gsk.cur), 5]),
     xaxt = "n")
lines(seq_along(these.gsk.cur), qs.results.mn[these.gsk.cur, 5], lty = 2)
axis(1, at = 1:8, labels = seq(5, 40, by = 5))
legend("topright", lty = c(1, 2), legend = c("EBF   ", "GSK   "))

plot(seq_along(these.ebf.fut), qs.results.mn[these.ebf.fut, 5], type = "l",
     main = "Future: Quantile score for q(0.99)",
     ylab = "Quantile score", xlab = "Knots",
     ylim = range(qs.results.mn[c(these.ebf.fut, these.gsk.fut), 5]),
     xaxt = "n")
lines(seq_along(these.gsk.fut), qs.results.mn[these.gsk.fut, 5], lty = 2)
axis(1, at = 1:8, labels = seq(5, 40, by = 5))
legend("topright", lty = c(1, 2), legend = c("EBF   ", "GSK   "))
dev.print(device = pdf, file = "./plots/precip-qs.pdf")
dev.off()

# posterior of time coefficient
load("./cv-results/ebf-current-35-all.RData")
beta1 <- fit$beta1[, 2]
beta2 <- fit$beta2[, 2]
load("./cv-results/gsk-current-35-all.RData")
beta1 <- cbind(beta1, fit$beta1[, 2])
beta2 <- cbind(beta2, fit$beta2[, 2])
load("./cv-results/ebf-future-35-all.RData")
beta1 <- cbind(beta1, fit$beta1[, 2])
beta2 <- cbind(beta2, fit$beta2[, 2])
load("./cv-results/gsk-future-35-all.RData")
beta1 <- cbind(beta1, fit$beta1[, 2])
beta2 <- cbind(beta2, fit$beta2[, 2])
quartz(width = 12, height = 6)
par(mfrow = c(1, 2))
boxplot(beta1, xlab = "", xaxt = "n")
axis(1, at = 1:4,
     labels = c("Current EBF", "Current GSK", "Future EBF", "Future GSK"))
boxplot(beta2, xlab = "", xaxt = "n")
axis(1, at = 1:4,
     labels = c("Current EBF", "Current GSK", "Future EBF", "Future GSK"))
dev.print(device = pdf, file = "./plots/precip-post-time.pdf")
dev.off()

# add in 21st row
this.row <- 21
this.qs <- qs.results[[this.row]]
this.bs <- bs.results[[this.row]]
qs.results.mn[this.row, ] <- apply(this.qs, 2, mean,
                                   na.rm = TRUE)
bs.results.mn[this.row, ] <- apply(this.bs, 2, mean,
                                   na.rm = TRUE)
qs.results.se[this.row, ] <- apply(this.qs, 2, sd,
                                   na.rm = TRUE) / sqrt(10)
bs.results.se[this.row, ] <- apply(this.bs, 2, sd,
                                   na.rm = TRUE) / sqrt(10)

process.col <- c(rep("ebf", nbases * 2), rep("gsk", nbases * 2 + 1))
margin.col  <- c(rep("ebf", nbases), rep("gsk", nbases),
                 rep("ebf", nbases), rep("gsk", nbases + 1))
bases.col   <- c(rep(bases, 4), 159)
results.df <- data.frame(process = as.factor(process.col),
                         margin = as.factor(margin.col),
                         bases = as.integer(bases.col),
                         bs95 = as.double(100 * bs.results.mn[, 1]),
                         bs99 = as.double(100 * bs.results.mn[, 2]),
                         qs95 = as.double(qs.results.mn[, 1]),
                         qs99 = as.double(qs.results.mn[, 5]))

quartz(width = 16, height = 8)
par(mfrow = c(1, 2))
# look at plots by process/margin
ylim <- range(results.df$bs95)
ylim[2] <- ylim[2] + 0.4
these <- results.df$process == "ebf" & results.df$margin == "ebf"
plot(results.df$bases[these], results.df$bs95[these], type = "l",
     # main = "Brier score (x 100) for exceeding q(0.95)",
     ylim = ylim, ylab = "Brier score (x 100)",
     xlab = "Number of knots/basis functions",
     col = "dodgerblue1", lty = 1, cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.5)
these <- results.df$process == "ebf" & results.df$margin == "gsk"
lines(results.df$bases[these], results.df$bs95[these],
      col = "dodgerblue1", lty = 2)
these <- results.df$process == "gsk" & results.df$margin == "ebf"
lines(results.df$bases[these], results.df$bs95[these],
      col = "firebrick1", lty = 1)
these <- results.df$process == "gsk" & results.df$margin == "gsk" & results.df$bases != 159
lines(results.df$bases[these], results.df$bs95[these],
      col = "firebrick1", lty = 2)
curve(results.df$bs95[21] + x * 0, from = 5, to = 25, lty = 1, col = "black",
      add = TRUE)

legend("topright", legend = c("Process: ebf Margin: ebf",
                              "Process: ebf Margin: gsk",
                              "Process: gsk Margin: ebf",
                              "Process: gsk Margin: gsk",
                              "Process: gsk Margin: gsk, knots at all sites"),
       lty = c(1, 2, 1, 2),
       col = c("dodgerblue1", "dodgerblue1", "firebrick1", "firebrick1", "black"),
       cex = 1.5)


# look at plots by process/margin
ylim <- range(results.df$bs99)
ylim[2] <- ylim[2] + 0.15
these <- results.df$process == "ebf" & results.df$margin == "ebf"
plot(results.df$bases[these], results.df$bs99[these], type = "l",
     # main = "Brier score (x 100) for exceeding q(0.99)",
     ylim = ylim, ylab = "Brier score (x 100)",
     xlab = "Number of knots/basis functions",
     col = "dodgerblue1", lty = 1, cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.5)
these <- results.df$process == "ebf" & results.df$margin == "gsk"
lines(results.df$bases[these], results.df$bs99[these],
      col = "dodgerblue1", lty = 2)
these <- results.df$process == "gsk" & results.df$margin == "ebf"
lines(results.df$bases[these], results.df$bs99[these],
      col = "firebrick1", lty = 1)
these <- results.df$process == "gsk" & results.df$margin == "gsk" & results.df$bases != 159
lines(results.df$bases[these], results.df$bs99[these],
      col = "firebrick1", lty = 2)
curve(results.df$bs99[21] + x * 0, from = 5, to = 25, lty = 1, col = "black",
      add = TRUE)

legend("topright", legend = c("Process: ebf Margin: ebf",
                              "Process: ebf Margin: gsk",
                              "Process: gsk Margin: ebf",
                              "Process: gsk Margin: gsk",
                              "Process: gsk Margin: gsk, knots at all sites"),
       lty = c(1, 2, 1, 2),
       col = c("dodgerblue1", "dodgerblue1", "firebrick1", "firebrick1", "black"),
       cex = 1.5)
dev.print(device = pdf, file = "plots/bs-mean-fire.pdf")
dev.off()

quartz(width = 16, height = 8)
par(mfrow = c(1, 2))
ylim <- range(results.df$qs95[results.df$bases != 25])
ylim[2] <- ylim[2] + 10
# look at plots by process/margin
these <- results.df$process == "ebf" & results.df$margin == "ebf" & results.df$bases != 25
plot(results.df$bases[these], results.df$qs95[these], type = "l",
     # main = "Quantile score for q(0.95)",
     ylim = ylim, ylab = "Quantile score",
     xlim = c(5, 25), xlab = "Number of knots/basis functions",
     col = "dodgerblue1", lty = 1, cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.5)
these <- results.df$process == "ebf" & results.df$margin == "gsk"
lines(results.df$bases[these], results.df$qs95[these],
      col = "dodgerblue1", lty = 2)
these <- results.df$process == "gsk" & results.df$margin == "ebf" & results.df$bases != 25
lines(results.df$bases[these], results.df$qs95[these],
      col = "firebrick1", lty = 1)
these <- results.df$process == "gsk" & results.df$margin == "gsk" & results.df$bases != 159
lines(results.df$bases[these], results.df$qs95[these],
      col = "firebrick1", lty = 2)
curve(results.df$qs95[21] + x * 0, from = 5, to = 25, lty = 1, col = "black",
      add = TRUE)
legend("topleft", legend = c("Process: ebf Margin: ebf",
                             "Process: ebf Margin: gsk",
                             "Process: gsk Margin: ebf",
                             "Process: gsk Margin: gsk",
                             "Process: gsk Margin: gsk, knots at all sites"),
       lty = c(1, 2, 1, 2),
       col = c("dodgerblue1", "dodgerblue1", "firebrick1", "firebrick1", "black"),
       cex = 1.5)


# look at plots by process/margin
ylim <- range(results.df$qs99[results.df$bases != 25])
ylim[2] <- ylim[2] + 10
these <- results.df$process == "ebf" & results.df$margin == "ebf" & results.df$bases != 25
plot(results.df$bases[these], results.df$qs99[these], type = "l",
     # main = "Quantile score for q(0.99)",
     ylim = ylim, ylab = "Quantile score",
     xlim = c(5, 25), xlab = "Number of knots/basis functions",
     col = "dodgerblue1", lty = 1, cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.5)
these <- results.df$process == "ebf" & results.df$margin == "gsk"
lines(results.df$bases[these], results.df$qs99[these],
      col = "dodgerblue1", lty = 2)
these <- results.df$process == "gsk" & results.df$margin == "ebf" & results.df$bases != 25
lines(results.df$bases[these], results.df$qs99[these],
      col = "firebrick1", lty = 1)
these <- results.df$process == "gsk" & results.df$margin == "gsk" & results.df$bases != 159
lines(results.df$bases[these], results.df$qs99[these],
      col = "firebrick1", lty = 2)
curve(results.df$qs99[21] + x * 0, from = 5, to = 25, lty = 1, col = "black",
      add = TRUE)
legend("topleft", legend = c("Process: ebf Margin: ebf",
                             "Process: ebf Margin: gsk",
                             "Process: gsk Margin: ebf",
                             "Process: gsk Margin: gsk",
                             "Process: gsk Margin: gsk, knots at all sites"),
       lty = c(1, 2, 1, 2),
       col = c("dodgerblue1", "dodgerblue1", "firebrick1", "firebrick1", "black"),
       cex = 1.5)
dev.print(device = pdf, file = "plots/qs-mean-fire.pdf")
dev.off()



t.test(qs.results[[4]][, 1], qs.results[[9]][, 1], paired = TRUE)
t.test(qs.results[[4]][, 2], qs.results[[9]][, 2], paired = TRUE)
t.test(qs.results[[4]][, 3], qs.results[[9]][, 3], paired = TRUE)
t.test(qs.results[[4]][, 4], qs.results[[9]][, 4], paired = TRUE)
t.test(qs.results[[4]][, 5], qs.results[[9]][, 5], paired = TRUE)

t.test(bs.results[[4]][, 2], bs.results[[9]][, 2], paired = TRUE)
t.test(bs.results[[5]][, 2], bs.results[[10]][, 2], paired = TRUE)


colnames(qs.results.mn) <- colnames(qs.results.se) <- probs.for.qs
these.rownames <- paste(rep(procs, each = nbases),
                        rep(bases, times = nprocs))
rownames(qs.results.mn) <- rownames(qs.results.se) <- these.rownames
rownames(bs.results) <- these.rownames

bs.results.mn <- apply(bs.results, 1, mean)
bs.results.se <- apply(bs.results, 1, sd) / sqrt(10)

round(qs.results.mn, 3)
round(qs.results.se, 3)

dev.new(width = 4, height = 3)
par(mfrow=c(1, 1), mar=c(5.1, 5.1, 2.1, 2.1))
matplot(x = probs.for.qs, y = t(qs.results.mn),
        col = c(rep("dodgerblue4", nbases),
                rep("firebrick4", nbases)),
        bg = c(rep("dodgerblue1", nbases),
               rep("firebrick1", nbases)),
        type = "b", pch = c(rep(21, nbases), rep(22, nbases)),
        lty = c(rep(1, nbases), rep(2, nbases)),
        cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex = 1.5, lwd=1.5,
        # main = "Average quantile score at selected quantiles",
        ylab = "Mean quantile score", xlab = "Quantile")
legend("topright", legend = c("ECB", "GKB"),
       lty = c(1, 2), pch = c(21, 22),
       col = c("dodgerblue4", "firebrick4"),
       pt.bg = c("dodgerblue1", "firebrick1"),
       title = "Spatial Process",
       lwd = 1.5, cex=1.5)
dev.print(device = pdf, file = "../../LaTeX/plots/qs-mean.pdf")
dev.off()

# look at timing - try to keep like computers together (in minutes)
timing <- read.table(file = "./cv-tables/timing.txt")
timing.mn <- apply(timing, 1, mean)

dev.new(width = 4, height = 3)
par(mfrow=c(1, 1), mar=c(5.1, 5.1, 2.1, 2.1))
ylim = range(timing.mn)
ylim[2] <- ylim[2] + 100
plot(bases, timing.mn[1:5], type = "b",
     ylim = ylim,
     xaxt = "n", xlab = "Number of basis functions", ylab = "Timing (seconds)",
     col = "dodgerblue4", bg = "dodgerblue1", pch = 21,
     # main = "Timing comparison of 100 iterations: Basis functions vs kernel",
     cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex = 1.5, lwd=1.5)
axis(1, at = bases, labels = TRUE)
lines(bases, timing.mn[6:10], lty = 2, type = "b", lwd = 1.5,
      col = "firebrick4", bg = "firebrick1", pch = 22, cex = 1.5)
legend("topleft", legend = c("ECB", "GKF"),
       lty = c(1, 2), lwd = 1.5, pch = c(21, 22),
       col = c("dodgerblue4", "firebrick4"),
       pt.bg = c("dodgerblue1", "firebrick1"),
       title = "Spatial Process",
       cex=1.5)
dev.print(device = pdf, file = "../../LaTeX/plots/timing.pdf")
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

