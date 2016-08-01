rm(list=ls())

#### 2d basis functions ####
library(fields)
library(splines)
library(Rcpp)
source(file = "package_load.R")
load(file = "../../code/analysis/fire/georgia_preprocess/georgia_map.RData")
load(file = "../../code/analysis/fire/georgia_preprocess/fire_data.RData")

ns <- ncol(Y)
nt <- nrow(Y)

d.act <- rdist(cents)
diag(d.act) <- 0
n <- nrow(cents)
thresh <- rep(0, ns)
counties <- rownames(cents)

#### standardize the locations ####
s <- cents
s.scale <- min(diff(range(s[, 1])), diff(range(s[, 2])))
s[, 1] <- (s[, 1] - min(s[, 1])) / s.scale
s[, 2] <- (s[, 2] - min(s[, 2])) / s.scale
neighbors <- 5
d.scale <- rdist(s)
diag(d.scale) <- 0

df.1 <- 5
df.2 <- 6
B.1 <- bs(s[, 1], df = df.1, Boundary.knots = c(-0.1, 1.1))
B.2 <- bs(s[, 2], df = df.2, Boundary.knots = c(-0.1, 1.1))

B <- matrix(NA, nrow = ns, ncol = df.1 * df.2)
for (i in 1:ncol(B.1)) {
  for (j in 1:ncol(B.2)) {
    B[, (i - 1) * df.2 + j] <- B.1[, i] * B.2[, j]
  }
}

# plots of the basis functions suggest that when colSums <  4, the basis function
# only averages over a small number of counties.
keep.bases <- which(colSums(B) > 5)
B <- B[, keep.bases]

p1 <- map.ga.ggplot(Y = B[, 1], counties = counties,
                    main = "2d b-spline 1",
                    midpoint = diff(range(B[, 1])) / 2,
                    limits = c(0, max(B[, 1])))
p2 <- map.ga.ggplot(Y = B[, 2], counties = counties,
                    main = "2d b-spline 2",
                    midpoint = diff(range(B[, 2])) / 2,
                    limits = c(0, max(B[, 2])))
p3 <- map.ga.ggplot(Y = B[, 3], counties = counties,
                    main = "2d b-spline 3",
                    midpoint = diff(range(B[, 3])) / 2,
                    limits = c(0, max(B[, 3])))
p4 <- map.ga.ggplot(Y = B[, 4], counties = counties,
                    main = "2d b-spline 4",
                    midpoint = diff(range(B[, 4])) / 2,
                    limits = c(0, max(B[, 4])))

p5 <- map.ga.ggplot(Y = B[, 5], counties = counties,
                    main = "2d b-spline 5",
                    midpoint = diff(range(B[, 5])) / 2,
                    limits = c(0, max(B[, 5])))
p6 <- map.ga.ggplot(Y = B[, 6], counties = counties,
                    main = "2d b-spline 6",
                    midpoint = diff(range(B[, 6])) / 2,
                    limits = c(0, max(B[, 6])))
p7 <- map.ga.ggplot(Y = B[, 7], counties = counties,
                    main = "2d b-spline 7",
                    midpoint = diff(range(B[, 7])) / 2,
                    limits = c(0, max(B[, 7])))
p8 <- map.ga.ggplot(Y = B[, 8], counties = counties,
                    main = "2d b-spline 8",
                    midpoint = diff(range(B[, 8])) / 2,
                    limits = c(0, max(B[, 8])))

p9 <- map.ga.ggplot(Y = B[, 9], counties = counties,
                    main = "2d b-spline 9",
                    midpoint = diff(range(B[, 9])) / 2,
                    limits = c(0, max(B[, 9])))
p10 <- map.ga.ggplot(Y = B[, 10], counties = counties,
                     main = "2d b-spline 10",
                     midpoint = diff(range(B[, 10])) / 2,
                     limits = c(0, max(B[, 10])))
p11 <- map.ga.ggplot(Y = B[, 11], counties = counties,
                     main = "2d b-spline 11",
                     midpoint = diff(range(B[, 11])) / 2,
                     limits = c(0, max(B[, 11])))
p12 <- map.ga.ggplot(Y = B[, 12], counties = counties,
                     main = "2d b-spline 12",
                     midpoint = diff(range(B[, 12])) / 2,
                     limits = c(0, max(B[, 12])))

p13 <- map.ga.ggplot(Y = B[, 13], counties = counties,
                     main = "2d b-spline 13",
                     midpoint = diff(range(B[, 13])) / 2,
                     limits = c(0, max(B[, 13])))

p14 <- map.ga.ggplot(Y = B[, 14], counties = counties,
                     main = "2d b-spline 14",
                     midpoint = diff(range(B[, 14])) / 2,
                     limits = c(0, max(B[, 14])))

multiplot(p1, p2, p3, p4, p5, p6, cols = 2)
dev.print(device = pdf, file = "plots/ga-bspline-1.pdf")

multiplot(p7, p8, p9, p10, p11, p12, cols = 2)
dev.print(device = pdf, file = "plots/ga-bspline-2.pdf")

multiplot(p13, p14, cols = 1)
dev.print(device = pdf, file = "plots/ga-bspline-3.pdf")

save(B, file = "basis")


rm(list=ls())

#### Plot EBF functions ####
library(fields)
library(splines)
library(Rcpp)
library(ggplot2)
library(gridExtra)
source(file = "package_load.R")
load(file = "../../code/analysis/fire/georgia_preprocess/georgia_map.RData")
load(file = "../../code/analysis/fire/georgia_preprocess/fire_data.RData")
load(file = "ebf-25-all.RData")
L <- 25

ns <- ncol(Y)
nt <- nrow(Y)

d.act <- rdist(cents)
diag(d.act) <- 0
n <- nrow(cents)
thresh <- rep(0, ns)
counties <- rownames(cents)

#### standardize the locations ####
s <- cents
s.scale <- min(diff(range(s[, 1])), diff(range(s[, 2])))
s[, 1] <- (s[, 1] - min(s[, 1])) / s.scale
s[, 2] <- (s[, 2] - min(s[, 2])) / s.scale
neighbors <- 5
d.scale <- rdist(s)
diag(d.scale) <- 0

v <- colSums(B.ebf) / ns
plot(1:L, cumsum(v), ylim = c(0, 1),
     main = paste("Fire analysis (", L, " knots)", sep = ""),
     ylab = "Cumulative relative contribution",
     xlab = "Knot")
plotname <- paste("plots/firev-", L, ".pdf", sep = "")
dev.print(device = pdf, file = plotname,
          width = 6, height = 6)

p1 <- map.ga.ggplot(Y = B.ebf[, 1], counties = counties,
                    main = "Basis function 1 (of 25)",
                    midpoint = diff(range(B.ebf[, 1])) / 2,
                    limits = c(0, max(B.ebf[, 1])))
p2 <- map.ga.ggplot(Y = B.ebf[, 2], counties = counties,
                    main = "Basis function 2 (of 25)",
                    midpoint = diff(range(B.ebf[, 2])) / 2,
                    limits = c(0, max(B.ebf[, 2])))
p3 <- map.ga.ggplot(Y = B.ebf[, 3], counties = counties,
                    main = "Basis function 3 (of 25)",
                    midpoint = diff(range(B.ebf[, 3])) / 2,
                    limits = c(0, max(B.ebf[, 3])))
p4 <- map.ga.ggplot(Y = B.ebf[, 4], counties = counties,
                    main = "Basis function 4 (of 25)",
                    midpoint = diff(range(B.ebf[, 4])) / 2,
                    limits = c(0, max(B.ebf[, 4])))
p5 <- map.ga.ggplot(Y = B.ebf[, 5], counties = counties,
                    main = "Basis function 5 (of 25)",
                    midpoint = diff(range(B.ebf[, 5])) / 2,
                    limits = c(0, max(B.ebf[, 5])))
p6 <- map.ga.ggplot(Y = B.ebf[, 6], counties = counties,
                    main = "Basis function 6 (of 25)",
                    midpoint = diff(range(B.ebf[, 6])) / 2,
                    limits = c(0, max(B.ebf[, 6])))

layout.mtx = matrix(1:6, nrow = 2, ncol = 3, byrow = TRUE)
panel <- arrangeGrob(p1, p2, p3, p4, p5, p6, ncol = 3, layout_matrix = layout.mtx)
ggsave(filename = "plots/fire-ebf-panel.pdf", panel, device = pdf, width = 13.5, height = 9)

ggsave(filename = "plots/fire-ebf-1.pdf", p1, device = pdf,
       width = 4.5, height = 4.5)
ggsave(filename = "plots/fire-ebf-2.pdf", p2, device = pdf,
       width = 4.5, height = 4.5)
ggsave(filename = "plots/fire-ebf-3.pdf", p3, device = pdf,
       width = 4.5, height = 4.5)
ggsave(filename = "plots/fire-ebf-4.pdf", p4, device = pdf,
       width = 4.5, height = 4.5)
ggsave(filename = "plots/fire-ebf-5.pdf", p5, device = pdf,
       width = 4.5, height = 4.5)
ggsave(filename = "plots/fire-ebf-6.pdf", p6, device = pdf,
       width = 4.5, height = 4.5)

# Eigenvectors
# Y is nt x ns so transposing
Y <- t(Y)
Y.mean <- apply(Y, 1, mean)
Y.center <- Y - Y.mean
Y <- t(Y)  # transpose back

# for correlation want ns x ns, so need cor(t(Y))
tY.center <- t(Y.center)
Y.eigvec <- eigen(cor(tY.center))$vectors
Y.eigval <- eigen(cor(tY.center))$values

# standardize eigenvalues
Y.eigval <- cumsum(Y.eigval) / sum(Y.eigval)

plot(Y.eigval[1:25], xlab = "Eigenvalue contribution", ylim = c(0, 1),
     ylab = "Cumulative relative contribution",
     main = "Precipitation analysis (25 eigenvalues)",
     cex.main = 1.5, cex.axis = 1.5, cex.lab = 1.5)
dev.print(device = pdf, file = "plots/firelambda-25.pdf")

e1 <- map.ga.ggplot(Y = Y.eigvec[, 1], counties = counties,
                    main = "Principal Component 1", midpoint = 0)
e2 <- map.ga.ggplot(Y = Y.eigvec[, 2], counties = counties,
                    main = "Principal Component 2", midpoint = 0)
e3 <- map.ga.ggplot(Y = Y.eigvec[, 3], counties = counties,
                    main = "Principal Component 3", midpoint = 0)
e4 <- map.ga.ggplot(Y = Y.eigvec[, 4], counties = counties,
                    main = "Principal Component 4", midpoint = 0)
e5 <- map.ga.ggplot(Y = Y.eigvec[, 5], counties = counties,
                    main = "Principal Component 5", midpoint = 0)
e6 <- map.ga.ggplot(Y = Y.eigvec[, 6], counties = counties,
                    main = "Principal Component 6", midpoint = 0)

layout.mtx = matrix(1:6, nrow = 2, ncol = 3, byrow = TRUE)
panel <- arrangeGrob(e1, e2, e3, e4, e5, e6, ncol = 3, layout_matrix = layout.mtx)
ggsave(filename = "plots/fire-eig-panel.pdf", panel, device = pdf, width = 13.5, height = 9)

ggsave(filename = "plots/fire-eig-1.pdf", e1, device = pdf,
       width = 4.5, height = 4.5)
ggsave(filename = "plots/fire-eig-2.pdf", e2, device = pdf,
       width = 4.5, height = 4.5)
ggsave(filename = "plots/fire-eig-3.pdf", e3, device = pdf,
       width = 4.5, height = 4.5)
ggsave(filename = "plots/fire-eig-4.pdf", e4, device = pdf,
       width = 4.5, height = 4.5)
ggsave(filename = "plots/fire-eig-5.pdf", e5, device = pdf,
       width = 4.5, height = 4.5)
ggsave(filename = "plots/fire-eig-6.pdf", e6, device = pdf,
       width = 4.5, height = 4.5)


# Eigenvectors of log(Y)
# Y is nt x ns so transposing
Y <- t(Y)
shift.Y <- Y + 0.01  # to avoid -Inf = log(0)
lY.mean <- apply(log(shift.Y), 1, mean)
lY.center <- log(shift.Y) - lY.mean
Y <- t(Y)  # transpose it back

# for correlation want ns x ns, so need cor(t(Y))
tlY.center <- t(lY.center)
lY.eigvec <- eigen(cor(tlY.center))$vectors
lY.eigval <- eigen(cor(tlY.center))$values

# standardize eigenvalues
lY.eigval <- cumsum(lY.eigval) / sum(lY.eigval)

plot(lY.eigval[1:25], xlab = "Eigenvalue contribution", ylim = c(0, 1),
     ylab = "Cumulative relative contribution",
     main = "Precipitation analysis (25 eigenvalues)",
     cex.main = 1.5, cex.axis = 1.5, cex.lab = 1.5)
dev.print(device = pdf, file = "plots/firelambda-log-25.pdf")

e1 <- map.ga.ggplot(Y = Y.eigvec[, 1], counties = counties,
                    main = "Eigenvector 1", midpoint = 0)
e2 <- map.ga.ggplot(Y = Y.eigvec[, 2], counties = counties,
                    main = "Eigenvector 2", midpoint = 0)
e3 <- map.ga.ggplot(Y = Y.eigvec[, 3], counties = counties,
                    main = "Eigenvector 3", midpoint = 0)
e4 <- map.ga.ggplot(Y = Y.eigvec[, 4], counties = counties,
                    main = "Eigenvector 4", midpoint = 0)
e5 <- map.ga.ggplot(Y = Y.eigvec[, 5], counties = counties,
                    main = "Eigenvector 5", midpoint = 0)
e6 <- map.ga.ggplot(Y = Y.eigvec[, 6], counties = counties,
                    main = "Eigenvector 6", midpoint = 0)

multiplot(e1, e2, e3, e4, e5, e6, cols = 3)
dev.print(device = pdf, width = 12, height = 8,
          file = "plots/fire-eig-log-panel.pdf")

e1
dev.print(device = pdf, width = 6, height = 6,
          file = "plots/fire-eig-log-1.pdf")
e2
dev.print(device = pdf, width = 6, height = 6,
          file = "plots/fire-eig-log-2.pdf")
e3
dev.print(device = pdf, width = 6, height = 6,
          file = "plots/fire-eig-log-3.pdf")
e4
dev.print(device = pdf, width = 6, height = 6,
          file = "plots/fire-eig-log-4.pdf")
e5
dev.print(device = pdf, width = 6, height = 6,
          file = "plots/fire-eig-log-5.pdf")
e6
dev.print(device = pdf, width = 6, height = 6,
          file = "plots/fire-eig-log-6.pdf")