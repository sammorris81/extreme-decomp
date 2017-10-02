rm(list=ls())
source(file = "./package_load.R", chdir = T)
library(gridExtra)
library(SpatialExtremes)

# setMKLthreads(5)
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
# temp <- fmadogram(data = t(Y[1:5, ]), coord = s[1:5, ], which = "ext")
# ec.temp <- matrix(1, nrow(Y[1:5, ]), nrow(Y[1:5, ]))
# ec.temp[lower.tri(ec.temp)] <- temp[, 3]
# ec.temp[upper.tri(ec.temp)] <- t(ec.temp)[upper.tri(ec.temp)]
# ec.hat[[1]] <-

for (fold in 1:nfolds) {
  Y.tst <- Y
  Y.tst[cv.idx[[fold]]] <- NA

  # # build ec matrix: ns x ns
  # ec <- get.pw.ec.fmado(Y = Y.tst)
  # ec.hat[[fold]] <- ec$ec
  this.ec <- fmadogram(data = t(Y.tst), coord = s, which = "ext")
  this.ec[this.ec >= 2] <- 2
  ec <- matrix(1, nrow(Y), nrow(Y))
  ec[lower.tri(ec)] <- this.ec[, 3]
  ec[upper.tri(ec)] <- t(ec)[upper.tri(ec)]
  ec.hat[[fold]] <- ec

  cat("finished fold:", fold, "\n")
}

save(cv.idx, ec.hat, file = "cv-extcoef.RData")

#### Try to precalculate the basis functions #########
#### Hoping to save a little time in the analysis ####
load("precip_preprocess.RData")
load("cv-extcoef.RData")
nfolds <- length(cv.idx)
# openblas.set.num.threads(4)
s.scale        <- s
s.scale.factor <- min(diff(range(s[, 1])), diff(range(s[, 2])))
s.min          <- apply(s, 2, min)
s.scale[, 1]   <- (s[, 1] - s.min[1]) / s.scale.factor
s.scale[, 2]   <- (s[, 2] - s.min[2]) / s.scale.factor
cents.grid     <- s.scale

nknots <- c(5, 10, 15, 20, 25, 30, 35, 40)

for (L in nknots) {
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
  set.seed(5687 + L)  # knots + L
  cat("Starting estimation of Gaussian kernels \n")
  knots <- cover.design(cents.grid, nd = L)$design
  B.gsk <- vector(mode = "list", length = nfolds)
  for (fold in 1:nfolds) {
    out   <- get.rho.alpha(EC = ec.hat[[fold]], s = s.scale, knots = knots,
                           init.rho = 0.3)
    B.gsk[[fold]] <- getW(rho = out$rho, dw2 = out$dw2)
    alphas[fold]  <- out$alpha

    cat("  Finished fold ", fold, " of ", nfolds, " for gsk. \n", sep = "")
  }

  filename <- paste("gsk-", L, ".RData", sep = "")
  save(B.gsk, alphas, knots, file = filename)

  cat("Finished L = ", L, ".\n", sep = "")
}

# par(mfrow = c(1, 2))
# quilt.plot(s.scale[, 1], s.scale[, 2], B.ebf[[1]][, 4],
#            nx = length(unique(s.scale[, 1])), ny = length(unique(s.scale[, 2])))
# quilt.plot(s.scale[, 1], s.scale[, 2], B.ebf1[[1]][, 4],
#            nx = length(unique(s.scale[, 1])), ny = length(unique(s.scale[, 2])))
#
# quilt.plot(s.scale[, 1], s.scale[, 2], B.gsk[[1]][, 4],
#            nx = length(unique(s.scale[, 1])), ny = length(unique(s.scale[, 2])))
# quilt.plot(s.scale[, 1], s.scale[, 2], B.gsk1[[1]][, 3],
#            nx = length(unique(s.scale[, 1])), ny = length(unique(s.scale[, 2])))
#
# #### looks like L = 35 is after things settle down.
# # get pairwise extremal coefficients
# # build ec matrix: ns x ns
# ec <- get.pw.ec.fmado(Y = Y)
# ec.hat <- ec$ec
# L <- 35
#
# # Empirical basis functions
# cat("Starting estimation of empirical basis functions \n")
# out       <- get.factors.EC(ec.hat, L = L, s = s.scale)
# B.ebf     <- out$est
# ec.smooth <- out$EC.smooth
# alpha     <- out$alpha
#
# filename <- paste("ebf-", L, "-all.RData", sep = "")
# save(B.ebf, ec.smooth, alpha, file = filename)
#
# # Gaussian kernel functions
# cat("Starting estimation of Gaussian kernels \n")
# set.seed(5687)
# knots <- cover.design(cents.grid, nd = L)$design
# out   <- get.rho.alpha(EC = ec.hat, s = s.scale, knots = knots)
# B.gsk <- getW(rho = out$rho, dw2 = out$dw2)
#
# filename <- paste("gsk-", L, "-all.RData", sep = "")
# save(B.gsk, alpha, knots, file = filename)

#### also do L = 25 is after things settle down.
# get pairwise extremal coefficients
# build ec matrix: ns x ns
this.ec <- fmadogram(data = t(Y), coord = s, which = "ext")
this.ec[this.ec >= 2] <- 2
ec <- matrix(1, nrow(Y), nrow(Y))
ec[lower.tri(ec)] <- this.ec[, 3]
ec[upper.tri(ec)] <- t(ec)[upper.tri(ec)]
ec.hat <- ec
# ec <- get.pw.ec.fmado(Y = Y)
# ec.hat <- ec$ec
L <- 25

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
set.seed(5687)
knots <- cover.design(cents.grid, nd = L)$design
out   <- get.rho.alpha(EC = ec.hat, s = s.scale, knots = knots)
B.gsk <- getW(rho = out$rho, dw2 = out$dw2)

filename <- paste("gsk-", L, "-all.RData", sep = "")
save(B.gsk, alpha, knots, file = filename)

# plot cumsum against basis function
L <- 25
file <- paste("ebf-", L, "-all.RData", sep = "")
load(file)
v <- colSums(B.ebf) / ns
plot(1:L, cumsum(v), ylim = c(0, 1),
     main = paste("Precipitation analysis (", L, " knots)", sep = ""),
     ylab = "Cumulative relative contribution",
     xlab = "Knot")
plotname <- paste("plots/precipv-", L, ".pdf", sep = "")
dev.print(device = pdf, file = plotname,
          width = 6, height = 6)

# L <- 35
# file <- paste("ebf-", L, "-all.RData", sep = "")
# load(file)
# v <- colSums(B.ebf) / ns
# plot(1:L, cumsum(v), ylim = c(0, 1),
#      main = paste("Precipitation analysis (", L, " knots)", sep = ""),
#      ylab = "Cumulative relative contribution",
#      xlab = "Knot")
# plotname <- paste("plots/precipv-", L, ".pdf", sep = "")
# dev.print(device = pdf, file = plotname,
#           width = 6, height = 6)


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
these <- sample(697, 50)
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
ylim <- range(Y[these, ])
ylim[1] <- 65
par(mfrow = c(1, 2))
for (i in 1:length(these)) {
  if (i == 1) {
    plot(Y[these[i], current], type = "l", ylim = ylim,
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
    plot(Y[these[i], future], type = "l", ylim = ylim,
         # main = "Yearly max daily precipitation (2039 - 2070)",
         ylab = "Max precipitation", xaxt = "n", xlab = "Year",
         col = colors[i], lwd = 1.5,
         cex.main = 1.5, cex.axis = 1.5, cex.lab = 1.5)
    axis(1, at = c(current), labels = year[future], cex.axis = 1.5)
  } else {
    lines(Y[these[i], future], col = colors[i], lwd = 1.5)
  }
}

legend("bottomright", col = color, lty = 1, cex = 1.5, lwd = 1.5,
       legend = c("Southeast", "Southwest", "Northeast", "Northwest"))
dev.print(device = pdf, file = "plots/precip-ts.pdf")

dev.off()

#### Generate basis function maps ####
load(file = "precip_preprocess.RData")
load(file = "ebf-25-all.RData")
ns <- nrow(Y)
nt <- ncol(Y)
basis.weight <- colSums(B.ebf) / ns
plot(cumsum(basis.weight), xlab = "Knot", ylim = c(0, 1),
     ylab = "Cumulative relative contribution",
     main = "Precipitation analysis (25 knots)",
     cex.main = 1.5, cex.axis = 1.5, cex.lab = 1.5)
dev.print(device = pdf, file = "plots/precipv-25.pdf",
          width = 4.5, height = 4.5)

p1 <- map.heatmap(lat = s[, 2], lon = s[, 1], data = B.ebf[, 1],
                  mainTitle = "Basis function 1 (of 25)")
p2 <- map.heatmap(lat = s[, 2], lon = s[, 1], data = B.ebf[, 2],
                  mainTitle = "Basis function 2 (of 25)")
p3 <- map.heatmap(lat = s[, 2], lon = s[, 1], data = B.ebf[, 3],
                  mainTitle = "Basis function 3 (of 25)")
p4 <- map.heatmap(lat = s[, 2], lon = s[, 1], data = B.ebf[, 4],
                  mainTitle = "Basis function 4 (of 25)")
p5 <- map.heatmap(lat = s[, 2], lon = s[, 1], data = B.ebf[, 5],
                  mainTitle = "Basis function 5 (of 25)")
p6 <- map.heatmap(lat = s[, 2], lon = s[, 1], data = B.ebf[, 6],
                  mainTitle = "Basis function 6 (of 25)")

layout.mtx = matrix(1:6, nrow = 2, ncol = 3, byrow = TRUE)
panel <- arrangeGrob(p1, p2, p3, p4, p5, p6, ncol = 3, layout_matrix = layout.mtx)
ggsave(filename = "./plots/precip-ebf-panel.pdf", panel, device = pdf,
       width = 13.5, height = 9)

ggsave(filename = "plots/precip-ebf-1.pdf", p1, device = pdf,
       width = 4.5, height = 4.5)
ggsave(filename = "plots/precip-ebf-2.pdf", p2, device = pdf,
       width = 4.5, height = 4.5)
ggsave(filename = "plots/precip-ebf-3.pdf", p3, device = pdf,
       width = 4.5, height = 4.5)
ggsave(filename = "plots/precip-ebf-4.pdf", p4, device = pdf,
       width = 4.5, height = 4.5)
ggsave(filename = "plots/precip-ebf-5.pdf", p5, device = pdf,
       width = 4.5, height = 4.5)
ggsave(filename = "plots/precip-ebf-6.pdf", p6, device = pdf,
       width = 4.5, height = 4.5)


# Eigenvectors
rm(list = ls())
source(file = "./package_load.R", chdir = T)
load(file = "precip_preprocess.RData")
load(file = "ebf-25-all.RData")

# for correlation want ns x ns, so need cor(t(Y))
Y.mean <- apply(Y, 1, mean)
Y.center <- Y - Y.mean
tY.center <- t(Y.center)

Y.eigvec <- eigen(cor(tY.center))$vectors
Y.eigval <- eigen(cor(tY.center))$values

# standardize eigenvalues
Y.eigval <- cumsum(Y.eigval) / sum(Y.eigval)

plot(Y.eigval[1:25], xlab = "Eigenvalue contribution", ylim = c(0, 1),
     ylab = "Cumulative relative contribution",
     main = "Precipitation analysis (25 eigenvalues)",
     cex.main = 1.5, cex.axis = 1.5, cex.lab = 1.5)
dev.print(device = pdf, file = "./plots/preciplambda-25.pdf", width = 6, height = 6)


e1 <- map.heatmap(lat = s[, 2], lon = s[, 1], data = Y.eigvec[, 1],
                  mainTitle = "Principal Component 1")
e2 <- map.heatmap(lat = s[, 2], lon = s[, 1], data = Y.eigvec[, 2],
                  mainTitle = "Principal Component 2")
e3 <- map.heatmap(lat = s[, 2], lon = s[, 1], data = Y.eigvec[, 3],
                  mainTitle = "Principal Component 3")
e4 <- map.heatmap(lat = s[, 2], lon = s[, 1], data = Y.eigvec[, 4],
                  mainTitle = "Principal Component 4")
e5 <- map.heatmap(lat = s[, 2], lon = s[, 1], data = Y.eigvec[, 5],
                  mainTitle = "Principal Component 5")
e6 <- map.heatmap(lat = s[, 2], lon = s[, 1], data = Y.eigvec[, 6],
                  mainTitle = "Principal Component 6")

layout.mtx = matrix(1:6, nrow = 2, ncol = 3, byrow = TRUE)
panel <- arrangeGrob(e1, e2, e3, e4, e5, e6, ncol = 3, layout_matrix = layout.mtx)
ggsave(filename = "./plots/precip-eig-panel.pdf", panel, device = pdf,
       width = 13.5, height = 9)

ggsave(filename = "plots/precip-eig-1.pdf", e1, device = pdf,
       width = 4.5, height = 4.5)
ggsave(filename = "plots/precip-eig-2.pdf", e2, device = pdf,
       width = 4.5, height = 4.5)
ggsave(filename = "plots/precip-eig-3.pdf", e3, device = pdf,
       width = 4.5, height = 4.5)
ggsave(filename = "plots/precip-eig-4.pdf", e4, device = pdf,
       width = 4.5, height = 4.5)
ggsave(filename = "plots/precip-eig-5.pdf", e5, device = pdf,
       width = 4.5, height = 4.5)
ggsave(filename = "plots/precip-eig-6.pdf", e6, device = pdf,
       width = 4.5, height = 4.5)

rm(list = ls())
library(gridExtra)
load("precip_preprocess.RData")
load("ebf-25-all.RData")

# NYC is at (-74.0059, 40.7128)
this.nyc.1 <- which(s[, 1] == -74.0625)
this.nyc.2 <- which(s[, 2] == s[600, 2])  # theres a weird rounding issue
nyc <- intersect(this.nyc.1, this.nyc.2)

# ATL is at (-84.3880, 33.7490)
this.atl.1 <- which(s[, 1] == -84.6875)
this.atl.2 <- which(s[, 2] == s[434, 2])
atl <- intersect(this.atl.1, this.atl.2)

# IAD is at (-77.0369, 38.9072)
this.iad.1 <- which(s[, 1] == -77.8125)
this.iad.2 <- which(s[, 2] == s[502, 2])
iad <- intersect(this.iad.1, this.iad.2)

# KNX is at (-83.9207, 35.9606)
this.knx.1 <- which(s[, 1] == -84.0625)
this.knx.2 <- which(s[, 2] == s[496, 2])
knx <- intersect(this.knx.1, this.knx.2)

# get pairwise ECs for NYC, ATL, and
nyc.ebf <- atl.ebf <- iad.ebf <- knx.ebf <- rep(0, nrow(s))
for (i in 1:nrow(s)) {
  nyc.ebf[i] <- sum((B.ebf[nyc, ]^(1 / alpha) + B.ebf[i, ]^(1 / alpha))^alpha)
  atl.ebf[i] <- sum((B.ebf[atl, ]^(1 / alpha) + B.ebf[i, ]^(1 / alpha))^alpha)
  iad.ebf[i] <- sum((B.ebf[iad, ]^(1 / alpha) + B.ebf[i, ]^(1 / alpha))^alpha)
  knx.ebf[i] <- sum((B.ebf[knx, ]^(1 / alpha) + B.ebf[i, ]^(1 / alpha))^alpha)
}
nyc.ebf[nyc] <- NA
atl.ebf[atl] <- NA
iad.ebf[iad] <- NA
knx.ebf[knx] <- NA

p1 <- map.heatmap(lat = s[, 2], lon = s[, 1], data = nyc.ebf, midpoint = 1.5,
                  mainTitle = "New York City, NY", zlim = c(1, 2))
p2 <- map.heatmap(lat = s[, 2], lon = s[, 1], data = atl.ebf, midpoint = 1.5,
                  mainTitle = "Atlanta, GA", zlim = c(1, 2))
p3 <- map.heatmap(lat = s[, 2], lon = s[, 1], data = iad.ebf, midpoint = 1.5,
                  mainTitle = "Washington, DC", zlim = c(1, 2))
p4 <- map.heatmap(lat = s[, 2], lon = s[, 1], data = knx.ebf, midpoint = 1.5,
                  mainTitle = "Knoxville, TN", zlim = c(1, 2))

layout.mtx = matrix(1:4, nrow = 2, ncol = 2, byrow = TRUE)
panel <- arrangeGrob(p1, p2, p3, p4, ncol = 2, layout_matrix = layout.mtx)
ggsave(filename = "./plots/pairwise-ec.pdf", panel, device = pdf,
       width = 9, height = 9)