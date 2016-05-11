rm(list=ls())
source(file = "./package_load.R", chdir = T)
library(compiler)
enableJIT(3)

setMKLthreads(5)
################################################################################
#### Load in the data ##########################################################
################################################################################
# Doing this first because want to do a little cleaning up of non necessary
# items in the RData file.
load(file = "precip.RData")
Y <- Yvec
# only keep s and Yvec
elev.old <- elev

# find the elevations
elev <- rep(0, nrow(s))
for (i in 1:nrow(s)) {
  row <- which(s1 == s[i, 1])
  col <- which(s2 == s[i, 2])
  elev[i] <- elev.old[row, col]
}
rm(Yvec, Ymat, s1, s2)

################################################################################
#### For computational restraints, only use s[, 1] > -90 #######################
################################################################################
# This is basically the cutoff used by Reich and Shaby (2012)
keep.these <- which(s[, 1] > -90)
s <- s[keep.these, ]
Y <- Y[keep.these, ]
elev <- elev[keep.these]
save(Y, s, elev, year, file = "precip_preprocess.RData")

################################################################################
#### Preprocess locations and data and setup cross-validation ##################
################################################################################
# get distance matrix
d <- rdist(s)
diag(d) <- 0
ns <- nrow(s)
nt <- ncol(Y)

# set up the 5 fold cross validation
n.tot <- nrow(Y) * ncol(Y)
set.seed(28)  # cv
nfolds <- 5  # picking 5 because data are max-stable and 64 years of data

# stratifying the selection so each training set location loses only 20% of
# observations
cv.idx <- get.cv.test.strat(data = Y, nfolds = nfolds, idx = 1)

# loop over the list for cross-validation and get the basis functions
# before getting started

ec.hat <- vector(mode = "list", length = nfolds)
for (fold in 1:nfolds) {
  Y.tst <- Y
  Y.tst[cv.idx[[fold]]] <- NA

  # build ec matrix: ns x ns
  ec <- get.pw.ec.fmado(Y = Y.tst)
  ec.hat[[fold]] <- ec$ec

  cat("finished fold:", fold, "\n")
}

save(cv.idx, ec.hat, file = "cv-extcoef.RData")

#### Try to precalculate the basis functions #########
#### Hoping to save a little time in the analysis ####
s.scale        <- s
s.scale.factor <- min(diff(range(s[, 1])), diff(range(s[, 2])))
s.min          <- apply(s, 2, min)
s.scale[, 1]   <- (s[, 1] - s.min[1]) / s.scale.factor
s.scale[, 2]   <- (s[, 2] - s.min[2]) / s.scale.factor
cents.grid     <- s.scale

nknots <- c(5, 10, 15, 20, 25, 30, 35, 40)

for (L in nknots[7:8]) {
  # Empirical basis functions
  cat("Starting estimation of empirical basis functions \n")
  alphas <- rep(0, nfolds)
  ec.smooth <- B.ebf <- vector(mode = "list", length = nfolds)
  for (fold in 1:nfolds) {
    out               <- get.factors.EC(ec.hat[[fold]], L = L, s = s.scale)
    B.ebf[[fold]]     <- out$est
    ec.smooth[[fold]] <- out$EC.smooth
    alphas[fold]      <- out$alpha

    cat("  Finished fold ", fold, " of ", nfolds, " for ebf. \n", sep = "")
  }

  filename <- paste("ebf-", L, ".RData", sep = "")
  save(B.ebf, ec.smooth, alphas, file = filename)

  # Gaussian kernel functions
  cat("Starting estimation of Gaussian kernels \n")
  knots <- cover.design(cents.grid, nd = L)$design
  B.gsk <- vector(mode = "list", length = nfolds)
  for (fold in 1:nfolds) {
    out   <- get.rho.alpha(EC = ec.hat[[fold]], s = s.scale, knots = knots)
    B.gsk[[fold]] <- getW(rho = out$rho, dw2 = out$dw2)

    cat("  Finished fold ", fold, " of ", nfolds, " for gsk. \n", sep = "")
  }

  filename <- paste("gsk-", L, ".RData", sep = "")
  save(B.gsk, alphas, knots, file = filename)

  cat("Finished L = ", L, ".\n", sep = "")
}

par(mfrow = c(1, 2))
quilt.plot(s.scale[, 1], s.scale[, 2], B.ebf[[1]][, 4],
           nx = length(unique(s.scale[, 1])), ny = length(unique(s.scale[, 2])))
quilt.plot(s.scale[, 1], s.scale[, 2], B.ebf1[[1]][, 4],
           nx = length(unique(s.scale[, 1])), ny = length(unique(s.scale[, 2])))

quilt.plot(s.scale[, 1], s.scale[, 2], B.gsk[[1]][, 4],
           nx = length(unique(s.scale[, 1])), ny = length(unique(s.scale[, 2])))
quilt.plot(s.scale[, 1], s.scale[, 2], B.gsk1[[1]][, 3],
           nx = length(unique(s.scale[, 1])), ny = length(unique(s.scale[, 2])))

#### looks like L = 35 is after things settle down.
# get pairwise extremal coefficients
# build ec matrix: ns x ns
ec <- get.pw.ec.fmado(Y = Y)
ec.hat <- ec$ec

# Empirical basis functions
cat("Starting estimation of empirical basis functions \n")
out       <- get.factors.EC(ec.hat, L = L, s = s.scale)
B.ebf     <- out$est
ec.smooth <- out$EC.smooth
alpha     <- out$alpha

filename <- paste("ebf-", L, "-all.RData", sep = "")
save(B.ebf, ec.smooth, alpha, file = filename)

# Gaussian kernel functions
cat("Starting estimation of Gaussian kernels \n")
knots <- cover.design(cents.grid, nd = L)$design
out   <- get.rho.alpha(EC = ec.hat, s = s.scale, knots = knots)
B.gsk <- getW(rho = out$rho, dw2 = out$dw2)

filename <- paste("gsk-", L, "-all.RData", sep = "")
save(B.gsk, alpha, knots, file = filename)


#### plot some of the basis functions ####
nx <- length(unique(s[, 1]))
ny <- length(unique(s[, 2]))
par(mfrow = c(3, 2))
quilt.plot(x = s[, 1], y = s[, 2], z = B.ebf[[1]][, 1], nx = nx, ny = ny)
quilt.plot(x = s[, 1], y = s[, 2], z = B.ebf[[1]][, 2], nx = nx, ny = ny)
quilt.plot(x = s[, 1], y = s[, 2], z = B.ebf[[1]][, 3], nx = nx, ny = ny)
quilt.plot(x = s[, 1], y = s[, 2], z = B.ebf[[1]][, 4], nx = nx, ny = ny)
quilt.plot(x = s[, 1], y = s[, 2], z = B.ebf[[1]][, 5], nx = nx, ny = ny)
quilt.plot(x = s[, 1], y = s[, 2], z = B.ebf[[1]][, 6], nx = nx, ny = ny)

#### plot some of the time series ####
library(colorspace)

set.seed(7568)  # plot
these <- sample(2622, 50)
color <- rainbow_hcl(n = 4)  # SE, SW, NE, NW

s.mid <- apply(s, 2, median)
s.these <- s[these, ]
colors <- rep(NA, length(these))
colors[s.these[, 1] >= s.mid[1] & s.these[, 2] < s.mid[2]] <- color[1]    # SE
colors[s.these[, 1] < s.mid[1] & s.these[, 2] < s.mid[2]] <- color[2]     # SW
colors[s.these[, 1] >=  s.mid[1] & s.these[, 2] >= s.mid[2]] <- color[3]  # NE
colors[s.these[, 1] <  s.mid[1] & s.these[, 2] >= s.mid[2]] <- color[4]   # NE

current <- 1:32
future  <- 33:64
quartz(width = 16, height = 8)
par(mfrow = c(1, 2))
for (i in 1:length(these)) {
  if (i == 1) {
    plot(Y[these[i], current], type = "l", ylim = range(Y[, current]),
         # main = "Yearly max daily precipitation (1969 - 2000)",
         ylab = "Max precipitation", xaxt = "n", xlab = "Year",
         col = colors[i], lwd = 1.5,
         cex.main = 1.5, cex.axis = 1.5, cex.lab = 1.5)
    axis(1, at = c(current), labels = year[current], cex.axis = 1.5)
  } else {
    lines(Y[these[i], 1:32], col = colors[i], lwd = 1.5)
  }
}

legend("bottomright", col = color, lty = 1, cex = 1.5, lwd = 1.5,
       legend = c("Southeast", "Southwest", "Northeast", "Northwest"))


for (i in 1:length(these)) {
  if (i == 1) {
    plot(Y[these[i], future], type = "l", ylim = range(Y[, future]),
         # main = "Yearly max daily precipitation (2039 - 2070)",
         ylab = "Max precipitation", xaxt = "n", xlab = "Year",
         col = colors[i], lwd = 1.5,
         cex.main = 1.5, cex.axis = 1.5, cex.lab = 1.5)
    axis(1, at = c(current), labels = year[current], cex.axis = 1.5)
  } else {
    lines(Y[these[i], future], col = colors[i], lwd = 1.5)
  }
}

legend("bottomright", col = color, lty = 1, cex = 1.5, lwd = 1.5,
       legend = c("Southeast", "Southwest", "Northeast", "Northwest"))
dev.print(device = pdf, file = "plots/precip-ts.pdf")

dev.off()

library(ggplot2)
library(gridExtra)
# just want to see if it looks as weird when we run all the data
ec <- get.pw.ec(Y = Y, qlim = c(0.90, 1), verbose = TRUE, update = 50)$ec
p.1 <- map.ga.ggplot(Y = ec[, 4],
                     main = paste("Extremal Coefficients full data"),
                     fill.legend = "EC")

p.2 <- map.ga.ggplot(Y = ec.hat[[1]][, 4],
                     main = paste("Extremal Coefficients cross validation"),
                     fill.legend = "EC")

p.3 <- map.ga.ggplot(Y = ec.hat[[2]][, 4],
                     main = paste("Extremal Coefficients cross validation"),
                     fill.legend = "EC")

p.4 <- map.ga.ggplot(Y = ec.hat[[4]][, 4],
                     main = paste("Extremal Coefficients cross validation"),
                     fill.legend = "EC")

# grid.arrange(p.1, p.2, ncol = 2, widths = c(1.5, 1.5),
#              top = "EC comparison for CV")

grid.arrange(p.1, p.2, p.3, p.4, ncol = 2, widths = c(1.5, 1.5),
             top = "EC comparison for CV")

p.1 <- map.ga.ggplot(Y = ec[, 10],
                     main = paste("Extremal Coefficients full data"),
                     fill.legend = "EC")

p.2 <- map.ga.ggplot(Y = ec.hat[[1]][, 10],
                     main = paste("Extremal Coefficients cross validation"),
                     fill.legend = "EC")

p.3 <- map.ga.ggplot(Y = ec.hat[[2]][, 10],
                     main = paste("Extremal Coefficients cross validation"),
                     fill.legend = "EC")

p.4 <- map.ga.ggplot(Y = ec.hat[[4]][, 10],
                     main = paste("Extremal Coefficients cross validation"),
                     fill.legend = "EC")

grid.arrange(p.1, p.2, p.3, p.4, ncol = 2, widths = c(1.5, 1.5),
             top = "EC comparison for CV")

d <- rdist(cents)
diag(d) <- 0
n <- nrow(cents)

# standardize the locations
s <- cents
s[, 1] <- (s[, 1] - min(s[, 1])) / diff(range(s[, 1]))
s[, 2] <- (s[, 2] - min(s[, 2])) / diff(range(s[, 2]))

# get a reasonable bandwidth for the kernel smoother
d <- rdist(s)
diag(d) <- 0
ksmooth.bw <- quantile(d[upper.tri(d)], probs = 0.05)

cat("Start basis function estimation \n")
# basis function estimates using only the training data
L <- 10
out       <- get.factors.EC(ec, L = L, s = s, bw = ksmooth.bw)
B.est     <- out$est
ec.smooth <- out$EC.smooth

plots <- vector(mode = "list", length = L)
for (i in 1:L) {
  title <- paste("basis function", i)
  plots[[i]] <- map.ga.ggplot(Y = B.est[, i], main = title,
                              midpoint = median(B.est[, i]))
}
multiplot(plotlist = plots, cols = 4)
