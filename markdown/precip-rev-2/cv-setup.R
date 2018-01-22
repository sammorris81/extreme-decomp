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

# scale the sites
s.scale        <- s
s.scale.factor <- min(diff(range(s[, 1])), diff(range(s[, 2])))
s.min          <- apply(s, 2, min)
s.scale[, 1]   <- (s[, 1] - s.min[1]) / s.scale.factor
s.scale[, 2]   <- (s[, 2] - s.min[2]) / s.scale.factor
save(Y, s, s.scale, elev, year, file = "precip_preprocess.RData")

################################################################################
#### Preprocess locations and data and setup cross-validation ##################
################################################################################
# compute the distances
d <- rdist(s.scale)
diag(d) <- 0
ns <- nrow(s.scale)
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

ec.hat <- ec.smooth <- vector(mode = "list", length = nfolds)
alpha.hats <-rep(0, nfolds)
for (fold in 1:nfolds) {
  Y.tst <- Y
  Y.tst[cv.idx[[fold]]] <- NA

  # build ec matrix: ns x ns
  ec.hat[[fold]] <- get.ec.fmad(Y = Y.tst, s = s.scale)
  ec.smooth[[fold]] <- KsmoothCV(ec.hat[[fold]], s.scale)$EC
  alpha.hats[fold] <- EstimateAlpha(ec.hat = ec.hat[[fold]], d = d, n0 = 50)

  cat("finished fold:", fold, "\n")
}
save(cv.idx, ec.hat, alpha.hats, ec.smooth, file = "cv-extcoef.RData")

#### Try to precalculate the basis functions #########
#### Hoping to save a little time in the analysis ####
load("precip_preprocess.RData")
load("cv-extcoef.RData")
nfolds <- length(cv.idx)
# openblas.set.num.threads(4)
cents.grid     <- s.scale

nknots <- c(5, 10, 15, 20, 25, 30, 35, 40)

for (L in nknots) {
  # Empirical basis functions
  cat("Starting estimation of empirical basis functions \n")
  B.ebf <- B.gsk <- vector(mode = "list", length = nfolds)

  for (fold in 1:nfolds) {
    out <- get.factors.EC(EC.smooth = ec.smooth[[fold]],
                          alpha.hat = alpha.hats[fold], L = L, s = s.scale,
                          maxit = 2000, n_starts = 3)
    B.ebf[[fold]] <- out$est
    cat("  Finished fold ", fold, " of ", nfolds, " for ebf. \n", sep = "")
  }

  filename <- paste("basis_functions/ebf-", L, ".RData", sep = "")
  save(B.ebf, file = filename)

  # Gaussian kernel functions
  set.seed(5687 + L)  # knots + L
  cat("Starting estimation of Gaussian kernels \n")
  knots <- cover.design(cents.grid, nd = L)$design
  B.gsk <- vector(mode = "list", length = nfolds)
  for (fold in 1:nfolds) {
    out <- get.rho(EC.smooth = ec.smooth[[fold]], alpha.hat = alpha.hats[fold],
                   s = s.scale, knots = knots,
                   init.rho = 0.3)
    B.gsk[[fold]] <- getW(rho = out$rho, dw2 = out$dw2)
    cat("  Finished fold ", fold, " of ", nfolds, " for gsk. \n", sep = "")
  }

  filename <- paste("basis_functions/gsk-", L, ".RData", sep = "")
  save(B.gsk, knots, file = filename)

  cat("Finished L = ", L, ".\n", sep = "")
}

# Plot pairwise extremal coefficients
quartz(width = 8, height = 8)
fmadogram(data = t(Y), coord = s, which = "ext", n.bins = 500)
dev.print(device = pdf, width = 8, height = 8, file = "plots/bin-500-fmad.pdf")
dev.off()

quartz(width = 8, height = 8)
fmadogram(data = t(Y), coord = s, which = "ext", n.bins = 1000)
dev.print(device = pdf, width = 8, height = 8, file = "plots/bin-1000-fmad.pdf")
dev.off()

quartz(width = 8, height = 8)
fmadogram(data = t(Y), coord = s, which = "ext", n.bins = 2500)
dev.print(device = pdf, width = 8, height = 8, file = "plots/bin-2500-fmad.pdf")
dev.off()

quartz(width = 8, height = 8)
fmadogram(data = t(Y), coord = s, which = "ext", n.bins = 5000)
dev.print(device = pdf, width = 8, height = 8, file = "plots/bin-5000-fmad.pdf")
dev.off()

#### Generate basis function maps ####
load(file = "precip_preprocess.RData")
load(file = "./basis_functions/ebf-10-all.RData")
ns <- nrow(Y)
nt <- ncol(Y)
basis.weight <- colSums(B.ebf) / ns
plot(cumsum(basis.weight), xlab = "Knot", ylim = c(0, 1),
     ylab = "Cumulative relative contribution",
     main = "Precipitation analysis (10 knots)",
     cex.main = 1.5, cex.axis = 1.5, cex.lab = 1.5)
dev.print(device = pdf, file = "plots/precipv-10.pdf",
          width = 4.5, height = 4.5)

p1 <- map.heatmap(lat = s[, 2], lon = s[, 1], data = B.ebf[, 1],
                  mainTitle = "Basis function 1 (of 10)")
p2 <- map.heatmap(lat = s[, 2], lon = s[, 1], data = B.ebf[, 2],
                  mainTitle = "Basis function 2 (of 10)")
p3 <- map.heatmap(lat = s[, 2], lon = s[, 1], data = B.ebf[, 3],
                  mainTitle = "Basis function 3 (of 10)")
p4 <- map.heatmap(lat = s[, 2], lon = s[, 1], data = B.ebf[, 4],
                  mainTitle = "Basis function 4 (of 10)")
p5 <- map.heatmap(lat = s[, 2], lon = s[, 1], data = B.ebf[, 5],
                  mainTitle = "Basis function 5 (of 10)")
p6 <- map.heatmap(lat = s[, 2], lon = s[, 1], data = B.ebf[, 6],
                  mainTitle = "Basis function 6 (of 10)")

layout.mtx = matrix(1:6, nrow = 2, ncol = 3, byrow = TRUE)
panel <- arrangeGrob(p1, p2, p3, p4, p5, p6, ncol = 3, layout_matrix = layout.mtx)
ggsave(filename = "./plots/precip-ebf-panel.pdf", panel, device = pdf,
       width = 13.5, height = 9)

ggsave(filename = "./plots/precip-ebf-1.pdf", p1, device = pdf,
       width = 4.5, height = 4.5)
ggsave(filename = "./plots/precip-ebf-2.pdf", p2, device = pdf,
       width = 4.5, height = 4.5)
ggsave(filename = "./plots/precip-ebf-3.pdf", p3, device = pdf,
       width = 4.5, height = 4.5)
ggsave(filename = "./plots/precip-ebf-4.pdf", p4, device = pdf,
       width = 4.5, height = 4.5)
ggsave(filename = "./plots/precip-ebf-5.pdf", p5, device = pdf,
       width = 4.5, height = 4.5)
ggsave(filename = "./plots/precip-ebf-6.pdf", p6, device = pdf,
       width = 4.5, height = 4.5)

# Eigenvectors
rm(list = ls())
source(file = "./package_load.R", chdir = T)
options(warn = 0)
load(file = "./precip_preprocess.RData")
load(file = "./basis_functions/ebf-10-all.RData")

# for correlation want ns x ns, so need cor(t(Y))
Y.mean <- apply(Y, 1, mean)
Y.center <- Y - Y.mean
tY.center <- t(Y.center)

Y.eigvec <- eigen(cor(tY.center))$vectors
Y.eigval <- eigen(cor(tY.center))$values

# standardize eigenvalues
Y.eigval <- cumsum(Y.eigval) / sum(Y.eigval)

plot(Y.eigval[1:10], xlab = "Eigenvalue contribution", ylim = c(0, 1),
     ylab = "Cumulative relative contribution",
     main = "Precipitation analysis (10 eigenvalues)",
     cex.main = 1.5, cex.axis = 1.5, cex.lab = 1.5)
dev.print(device = pdf, file = "./plots/preciplambda-10.pdf", width = 6, height = 6)


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

ggsave(filename = "./plots/precip-eig-1.pdf", e1, device = pdf,
       width = 4.5, height = 4.5)
ggsave(filename = "./plots/precip-eig-2.pdf", e2, device = pdf,
       width = 4.5, height = 4.5)
ggsave(filename = "./plots/precip-eig-3.pdf", e3, device = pdf,
       width = 4.5, height = 4.5)
ggsave(filename = "./plots/precip-eig-4.pdf", e4, device = pdf,
       width = 4.5, height = 4.5)
ggsave(filename = "./plots/precip-eig-5.pdf", e5, device = pdf,
       width = 4.5, height = 4.5)
ggsave(filename = "./plots/precip-eig-6.pdf", e6, device = pdf,
       width = 4.5, height = 4.5)