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
B <- bspline.2d(s = cents, scale = TRUE, df.x = 5, df.y = 5)

map.ga.ggplot(Y = B[, 25], counties = counties, 
              main = "2d b-spline 1", 
              midpoint = diff(range(B[, 25])) / 2, 
              limits = c(0, max(B[, 25])))

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

save(B, file = "2d-bspline-ga.RData")
