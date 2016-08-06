rm(list = ls())

load(file = "../../code/analysis/fire/georgia_preprocess/fire_data.RData")
Y <- t(Y)
load(file = "cv-extcoef.RData")
nfolds <- length(cv.idx)
procs <- c("ebf", "gsk")  # what process determines the spatial structure
nprocs <- length(procs)
# margs <- c("ebf", "gsk")  # basis functions for the marginal distributions
margs <- "gsk"
nmargs <- length(margs)
bases   <- c(5, 10, 15, 20, 25, 30, 35, 40)
nbases  <- length(bases)
probs.for.qs <- c(0.95, 0.96, 0.97, 0.98, 0.99, 0.995)  # always check fitmodel
probs.for.bs <- c(0.95, 0.99)
nprobs.qs <- length(probs.for.qs)
nprobs.bs <- length(probs.for.bs)

files <- list.files(path = "cv-tables/")
files <- files[-c(41, 122)]
# each element of these lists is a matrix - including an extra for gsk-gsk-all
qs.results <- vector(mode = "list", length = nbases * nprocs * nmargs + 1)
bs.results <- vector(mode = "list", length = nbases * nprocs * nmargs + 1)

for (b in 1:(nbases * nmargs * nprocs + 1)) {
  qs.results[[b]] <- matrix(NA, nfolds, nprobs.qs)
  bs.results[[b]] <- matrix(NA, nfolds, nprobs.bs)
  colnames(qs.results[[b]]) <- probs.for.qs
  colnames(bs.results[[b]]) <- probs.for.bs
  rownames(qs.results[[b]]) <- paste("fold:", 1:nfolds)
  rownames(bs.results[[b]]) <- paste("fold:", 1:nfolds)
}

# timing is a data.frame that contains time, hostname, basis, and fold
timing <- data.frame(timing = double(), hostname = factor(),
                     proc = factor(), margin = factor(),
                     basis = factor(), fold = factor())
rownames <- rep(0, nbases * nmargs * nprocs + 1)
for (i in 1:(length(files))) {  # last file is timing.txt
  split     <- unlist(strsplit(unlist(strsplit(files[i], "-")), "[.]"))
  # files are named by the number of basis functions which skips numbers
  proc.idx  <- which(procs == split[1])
  margin.idx <- which(margs == split[2])
  basis.idx <- which(bases == split[3])

  if (as.numeric(split[3]) == 159) {
    idx <- length(rownames)
  } else {
    idx       <- (proc.idx - 1) * (nbases * nmargs) +
                 (margin.idx - 1) * nbases + basis.idx
  }

  rownames[idx] <- paste(split[1], "-", split[3], sep = "")

  fold      <- as.numeric(split[4])
  table.set <- read.table(paste("cv-tables/", files[i], sep = ""),
                          stringsAsFactors = FALSE)

  # first extract the timing information from the end of the vector
  timing.tail <- tail(table.set$x, 2)
  timing.row <- data.frame(timing = as.numeric(timing.tail[1]),
                           host = timing.tail[2],
                           proc = split[1], margin = split[2],
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
qs.results.mn <- qs.results.se <- matrix(NA, nbases * nmargs * 2 + 1, nprobs.qs)
bs.results.mn <- bs.results.se <- matrix(NA, nbases * nmargs * 2 + 1, nprobs.bs)
# idx: 1 - 5: ebf spatial, ebf marginal
# idx: 6 - 10: ebf spatial, gsk marginal
# idx: 11 - 15: gsk spatial, ebf marginal
# idx: 16 - 20: gsk spatial, gsk marginal
done.sets <- 1:10
for (p in 1:nprocs) {
  for (m in 1:nmargs) {
    for (b in 1:nbases) {
      this.row <- (p - 1) * (nbases * nmargs) +
                  (m - 1) * nbases + b
      this.qs <- qs.results[[this.row]]
      this.bs <- bs.results[[this.row]]
      qs.results.mn[this.row, ] <- apply(this.qs[done.sets, ], 2, mean,
                                         na.rm = TRUE)
      bs.results.mn[this.row, ] <- apply(this.bs[done.sets, ], 2, mean,
                                         na.rm = TRUE)
      qs.results.se[this.row, ] <- apply(this.qs[done.sets, ], 2, sd,
                                         na.rm = TRUE) / sqrt(10)
      bs.results.se[this.row, ] <- apply(this.bs[done.sets, ], 2, sd,
                                         na.rm = TRUE) / sqrt(10)
    }
  }
}

rownames(bs.results.mn) <- rownames
rownames(qs.results.mn) <- rownames

round(bs.results.mn * 100, 3)
round(qs.results.mn[, c(1, 5)], 3)

# Brier scores
these.ebf <- 1:8
these.gsk <- 9:16

quartz(width = 8, height = 8)
ylim <- range(bs.results.mn[c(these.ebf, these.gsk), ] * 100)
ylim[2] <- ylim[2] + 1.1  # give space for legend
plot(seq_along(these.ebf), bs.results.mn[these.ebf, 1] * 100, type = "b",
     # main = "Brier score (x 100)",
     ylab = "Brier score (x 100)", xlab = "Knots",
     ylim = ylim, xaxt = "n", bg = "dodgerblue1", pch = 21,
     cex = 1.5, cex.lab = 1.5, cex.axis = 1.5, lwd = 1.25)
lines(seq_along(these.gsk), bs.results.mn[these.gsk, 1] * 100, lty = 1,
      type = "b", bg = "firebrick1", pch = 21, cex = 1.5, lwd = 1.25)
lines(seq_along(these.ebf), bs.results.mn[these.ebf, 2] * 100, lty = 3,
      type = "b", bg = "dodgerblue1", pch = 21, cex = 1.5, lwd = 1.25)
lines(seq_along(these.gsk), bs.results.mn[these.gsk, 2] * 100, lty = 3,
      type = "b", bg = "firebrick1", pch = 21, cex = 1.5, lwd = 1.25)
axis(1, at = 1:8, labels = seq(5, 40, by = 5), cex.lab = 1.5, cex.axis = 1.5)
legend("topright",
       pt.bg = c("dodgerblue1", "firebrick1", "dodgerblue1", "firebrick1"),
       pch = c(21, 21, 21, 21), pt.lwd = 1.25, pt.cex = 1.5,
       cex = 1.5, lwd = 1.25, lty = c(1, 1, 3, 3),
       legend = c("q(0.95): EBF", "q(0.95): GSK",
                  "q(0.99): EBF", "q(0.99): GSK"))
dev.print(device = pdf, width = 8, height = 8, file = "plots/fire-bs.pdf")
dev.off()

quartz(width = 8, height = 8)
ylim <- range(qs.results.mn[c(these.ebf, these.gsk), c(1, 5)])
ylim[2] <- ylim[2] + 15  # give space for legend
plot(seq_along(these.ebf), qs.results.mn[these.ebf, 1], type = "b",
     # main = "Brier score (x 100)",
     ylab = "Quantile Score", xlab = "Knots",
     ylim = ylim, xaxt = "n", bg = "dodgerblue1", pch = 21,
     cex = 1.5, cex.lab = 1.5, cex.axis = 1.5, lwd = 1.25)
lines(seq_along(these.gsk), qs.results.mn[these.gsk, 1], lty = 1,
      type = "b", bg = "firebrick1", pch = 21, cex = 1.5, lwd = 1.25)
lines(seq_along(these.ebf), qs.results.mn[these.ebf, 5], lty = 3,
      type = "b", bg = "dodgerblue1", pch = 21, cex = 1.5, lwd = 1.25)
lines(seq_along(these.gsk), qs.results.mn[these.gsk, 5], lty = 3,
      type = "b", bg = "firebrick1", pch = 21, cex = 1.5, lwd = 1.25)
axis(1, at = 1:8, labels = seq(5, 40, by = 5), cex.lab = 1.5, cex.axis = 1.5)
legend("topright",
       pt.bg = c("dodgerblue1", "firebrick1", "dodgerblue1", "firebrick1"),
       pch = c(21, 21, 21, 21), pt.lwd = 1.25, pt.cex = 1.5,
       cex = 1.5, lwd = 1.25, lty = c(1, 1, 3, 3),
       legend = c("q(0.95): EBF", "q(0.95): GSK",
                  "q(0.99): EBF", "q(0.99): GSK"))
dev.print(device = pdf, width = 8, height = 8, file = "plots/fire-qs.pdf")
dev.off()

round(bs.results.mn * 100, 3)
round(qs.results.mn[, c(1, 5)], 3)

#### Panel for beta_time EBF
rm(list = ls())
source(file = "./package_load.R", chdir = T)
library(fields)
load("./cv-results/ebf-gsk-35-all.RData")
load("./ebf-35-all.RData")
load("./gsk-35-all.RData")
#### spatial setup ####
# get Georgia coordinates from georgia_preprocess in code/analysis/fire
load(file = "../../code/analysis/fire/georgia_preprocess/fire_data.RData")
Y <- t(Y)
load(file = "../../code/analysis/fire/georgia_preprocess/georgia_map.RData")

# get candidate knot grid for Gaussian kernel functions
grid.x <- seq(min(cents[, 1]), max(cents[, 1]), length = 100)
grid.y <- seq(min(cents[, 2]), max(cents[, 2]), length = 100)
cents.grid <- as.matrix(expand.grid(grid.x, grid.y))
inGA <- map.where("state", x = cents.grid[, 1], y = cents.grid[, 2])
cents.grid <- cents.grid[inGA == "georgia", ]
cents.grid <- cents.grid[rowSums(is.na(cents.grid)) == 0, ]

# standardize the locations
s <- cents
s.scale <- min(diff(range(s[, 1])), diff(range(s[, 2])))
s.min   <- apply(s, 2, min)
s[, 1] <- (s[, 1] - s.min[1]) / s.scale
s[, 2] <- (s[, 2] - s.min[2]) / s.scale
cents.grid[, 1] <- (cents.grid[, 1] - s.min[1]) / s.scale
cents.grid[, 2] <- (cents.grid[, 2] - s.min[2]) / s.scale

dw2 <- rdist(s, knots)^2
dw2[dw2 < 1e-4] <- 0
ns <- nrow(Y)
niters <- length(fit$bw)
L <- nrow(knots)
beta.mu.time <- beta.sig.time <- matrix(0, nrow = niters, ncol = ns)
for (i in 1:niters) {
  B <- makeW(dw2 = dw2, rho = fit$bw[i])
  beta.mu.time[i, ] <- fit$beta1[i, 2] + B %*% fit$beta1[i, (L + 3):(2 * L + 2)]
  beta.sig.time[i, ] <- fit$beta2[i, 2] + B %*% fit$beta2[i, (L + 3):(2 * L + 2)]
}

beta.mu.mean <- apply(beta.mu.time, 2, mean)
beta.mu.sd   <- apply(beta.mu.time, 2, sd)
beta.mu.ppos <- apply(beta.mu.time > 0, 2, mean)
beta.sig.mean <- apply(beta.sig.time, 2, mean)
beta.sig.sd   <- apply(beta.sig.time, 2, sd)
beta.sig.ppos <- apply(beta.sig.time > 0, 2, mean)

title <- bquote(paste("EBF: Posterior mean of ",
                      beta[paste(mu, ", time")]))
legend.title <- bquote(paste("E(", beta[paste(mu, ", time")], ")"))
p1 <- map.ga.ggplot(Y = beta.mu.mean, midpoint = 0,
              main = title, fill.legend = legend.title)

title <- bquote(paste("EBF: Posterior SD of ",
                      beta[paste(mu, ", time")]))
legend.title <- bquote(paste("SD(", beta[paste(mu, ", time")], ")"))
p2 <- map.ga.ggplot(Y = beta.mu.sd, midpoint = median(beta.mu.sd),
              main = title, fill.legend = legend.title)

title <- bquote(paste("EBF: Posterior P(",
                      beta[paste(mu, ", time")], "> 0)"))
legend.title <- bquote(paste("P(", beta[paste(mu, ", time")], "> 0)"))
p3 <- map.ga.ggplot(Y = beta.mu.ppos, midpoint = 0.5, limits = c(0, 1),
              main = title, fill.legend = legend.title)

title <- bquote(paste("EBF: Posterior mean of ",
                      beta[paste(sigma, ", time")]))
legend.title <- bquote(paste("E(", beta[paste(sigma, ", time")], ")"))
p4 <- map.ga.ggplot(Y = beta.sig.mean, midpoint = 0,
              main = title, fill.legend = legend.title)

title <- bquote(paste("EBF: Posterior SD of ",
                      beta[paste(sigma, ", time")]))
legend.title <- bquote(paste("SD(", beta[paste(sigma, ", time")], ")"))
p5 <- map.ga.ggplot(Y = beta.sig.sd, midpoint = median(beta.sig.sd),
              main = title, fill.legend = legend.title)

title <- bquote(paste("EBF: Posterior P(",
                      beta[paste(sigma, ", time")], "> 0)"))
legend.title <- bquote(paste("P(", beta[paste(sigma, ", time")], "> 0)"))
p6 <- map.ga.ggplot(Y = beta.sig.ppos, midpoint = 0.5, limits = c(0, 1),
              main = title, fill.legend = legend.title)
multiplot(p1, p2, p3, p4, p5, p6, cols = 2)

rm(list = ls())
source(file = "./package_load.R", chdir = T)
library(fields)
load("./cv-results/gsk-gsk-35-all.RData")
load("./ebf-35-all.RData")
load("./gsk-35-all.RData")
#### spatial setup ####
# get Georgia coordinates from georgia_preprocess in code/analysis/fire
load(file = "../../code/analysis/fire/georgia_preprocess/fire_data.RData")
Y <- t(Y)
load(file = "../../code/analysis/fire/georgia_preprocess/georgia_map.RData")

# get candidate knot grid for Gaussian kernel functions
grid.x <- seq(min(cents[, 1]), max(cents[, 1]), length = 100)
grid.y <- seq(min(cents[, 2]), max(cents[, 2]), length = 100)
cents.grid <- as.matrix(expand.grid(grid.x, grid.y))
inGA <- map.where("state", x = cents.grid[, 1], y = cents.grid[, 2])
cents.grid <- cents.grid[inGA == "georgia", ]
cents.grid <- cents.grid[rowSums(is.na(cents.grid)) == 0, ]

# standardize the locations
s <- cents
s.scale <- min(diff(range(s[, 1])), diff(range(s[, 2])))
s.min   <- apply(s, 2, min)
s[, 1] <- (s[, 1] - s.min[1]) / s.scale
s[, 2] <- (s[, 2] - s.min[2]) / s.scale
cents.grid[, 1] <- (cents.grid[, 1] - s.min[1]) / s.scale
cents.grid[, 2] <- (cents.grid[, 2] - s.min[2]) / s.scale

dw2 <- rdist(s, knots)^2
dw2[dw2 < 1e-4] <- 0
ns <- nrow(Y)
niters <- length(fit$bw)
L <- nrow(knots)
beta.mu.time <- beta.sig.time <- matrix(0, nrow = niters, ncol = ns)
for (i in 1:niters) {
  B <- makeW(dw2 = dw2, rho = fit$bw[i])
  beta.mu.time[i, ] <- fit$beta1[i, 2] + B %*% fit$beta1[i, (L + 3):(2 * L + 2)]
  beta.sig.time[i, ] <- fit$beta2[i, 2] + B %*% fit$beta2[i, (L + 3):(2 * L + 2)]
}

beta.mu.mean <- apply(beta.mu.time, 2, mean)
beta.mu.sd   <- apply(beta.mu.time, 2, sd)
beta.mu.ppos <- apply(beta.mu.time > 0, 2, mean)
beta.sig.mean <- apply(beta.sig.time, 2, mean)
beta.sig.sd   <- apply(beta.sig.time, 2, sd)
beta.sig.ppos <- apply(beta.sig.time > 0, 2, mean)

title <- bquote(paste("GSK: Posterior mean of ",
                      beta[paste(mu, ", time")]))
legend.title <- bquote(paste("E(", beta[paste(mu, ", time")], ")"))
p1 <- map.ga.ggplot(Y = beta.mu.mean, midpoint =0,
              main = title, fill.legend = legend.title)

title <- bquote(paste("GSK: Posterior SD of ",
                      beta[paste(mu, ", time")]))
legend.title <- bquote(paste("SD(", beta[paste(mu, ", time")], ")"))
p2 <- map.ga.ggplot(Y = beta.mu.sd, midpoint = median(beta.mu.sd),
              main = title, fill.legend = legend.title)

title <- bquote(paste("GSK: Posterior P(",
                      beta[paste(mu, ", time")], "> 0)"))
legend.title <- bquote(paste("P(", beta[paste(mu, ", time")], "> 0)"))
p3 <- map.ga.ggplot(Y = beta.mu.ppos, midpoint = 0.5, limits = c(0, 1),
              main = title, fill.legend = legend.title)

title <- bquote(paste("GSK: Posterior mean of ",
                      beta[paste(sigma, ", time")]))
legend.title <- bquote(paste("E(", beta[paste(sigma, ", time")], ")"))
p4 <- map.ga.ggplot(Y = beta.sig.mean, midpoint = 0,
              main = title, fill.legend = legend.title)

title <- bquote(paste("GSK: Posterior SD of ",
                      beta[paste(sigma, ", time")]))
legend.title <- bquote(paste("SD(", beta[paste(sigma, ", time")], ")"))
p5 <- map.ga.ggplot(Y = beta.sig.sd, midpoint = median(beta.sig.sd),
              main = title, fill.legend = legend.title)

title <- bquote(paste("GSK: Posterior P(",
                      beta[paste(sigma, ", time")], "> 0)"))
legend.title <- bquote(paste("P(", beta[paste(sigma, ", time")], "> 0)"))
p6 <- map.ga.ggplot(Y = beta.sig.ppos, midpoint = 0.5, limits = c(0, 1),
              main = title, fill.legend = legend.title)
multiplot(p1, p2, p3, p4, p5, p6, cols = 2)

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

process.col <- c(rep("ebf", nbases), rep("gsk", nbases))
# margin.col  <- c(rep("ebf", nbases), rep("gsk", nbases),
#                  rep("ebf", nbases), rep("gsk", nbases))
bases.col   <- rep(bases, 2)
results.df <- data.frame(process = as.factor(process.col),
                         # margin = as.factor(margin.col),
                         bases = as.integer(bases.col),
                         bs95 = as.double(100 * bs.results.mn[1:16, 1]),
                         bs99 = as.double(100 * bs.results.mn[1:16, 2]),
                         qs95 = as.double(qs.results.mn[1:16, 1]),
                         qs99 = as.double(qs.results.mn[1:16, 5]))

quartz(width = 16, height = 8)
par(mfrow = c(1, 2))
# look at plots by process/margin
ylim <- range(results.df$bs95)
ylim[2] <- ylim[2] + 0.02
these <- results.df$process == "ebf"
plot(results.df$bases[these], results.df$bs95[these], type = "l",
     main = "Brier score (x 100) for exceeding q(0.95)",
     ylim = ylim, ylab = "Brier score (x 100)",
     xlab = "Number of knots/basis functions",
     lty = 1, cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.5)
these <- results.df$process == "gsk"
lines(results.df$bases[these], results.df$bs95[these], lty = 2)

legend("topright", legend = c("EBF", "GSK"),
       lty = c(1, 2),
       cex = 1.5)


# look at plots by process/margin
ylim <- range(results.df$bs99)
ylim[2] <- ylim[2] + 0.02
these <- results.df$process == "ebf"
plot(results.df$bases[these], results.df$bs99[these], type = "l",
     main = "Brier score (x 100) for exceeding q(0.99)",
     ylim = ylim, ylab = "Brier score (x 100)",
     xlab = "Number of knots/basis functions",
     lty = 1, cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.5)
these <- results.df$process == "gsk"
lines(results.df$bases[these], results.df$bs99[these], lty = 2)

legend("topright", legend = c("EBF", "GSK"),
       lty = c(1, 2),
       cex = 1.5)
dev.print(device = pdf, file = "plots/fire-bs.pdf")
dev.off()

quartz(width = 16, height = 8)
par(mfrow = c(1, 2))
ylim <- range(results.df$qs95)
ylim[2] <- ylim[2] + 2
# look at plots by process/margin
these <- results.df$process == "ebf"
plot(results.df$bases[these], results.df$qs95[these], type = "l",
     main = "Quantile score for q(0.95)",
     ylim = ylim, ylab = "Quantile score",
     xlim = c(5, 40), xlab = "Number of knots/basis functions",
     lty = 1, cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.5)
these <- results.df$process == "gsk"
lines(results.df$bases[these], results.df$qs95[these], lty = 2)

legend("topright", legend = c("EBF", "GSK"),
       lty = c(1, 2), cex = 1.5)

# look at plots by process/margin
ylim <- range(results.df$qs99)
ylim[2] <- ylim[2] + 2
these <- results.df$process == "ebf"
plot(results.df$bases[these], results.df$qs99[these], type = "l",
     main = "Quantile score for q(0.99)",
     ylim = ylim, ylab = "Quantile score",
     xlim = c(5, 40), xlab = "Number of knots/basis functions",
     lty = 1, cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.5)
these <- results.df$process == "gsk"
lines(results.df$bases[these], results.df$qs99[these], lty = 2)

legend("topright", legend = c("EBF", "GSK"),
       lty = c(1, 2),
       cex = 1.5)
dev.print(device = pdf, file = "plots/fire-qs.pdf")
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

