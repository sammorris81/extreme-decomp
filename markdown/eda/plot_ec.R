rm(list=ls())

library(splines)
library(maps)
library(maptools)
library(fields)
library(ggplot2)
library(gridExtra)
library(rapport)
source(file = "../../code/analysis/fire/adj.R", chdir = TRUE)
source(file = "../../code/R/auxfunctions.R", chdir = TRUE)
source(file = "../../code/R/PCAX.R", chdir = TRUE)
# load(file = "../code/analysis/fire/gaCntyFires.RData")
load(file = "../../code/analysis/fire/georgia_preprocess/fire_data.RData")

# get the Georgia map and coordinates
# from georgia_preprocess in code/analysis/fire
load(file = "../../code/analysis/fire/georgia_preprocess/georgia_map.RData")
d <- rdist(cents)
diag(d) <- 0
n <- nrow(cents)

# chi[i, j] is the chi statistic between sites i and j
load(file = "../../code/analysis/fire/georgia_preprocess/chi.RData")
chi.hat <- ifelse(chi <= 0, 0, chi)
ec.hat  <- 2 - chi.hat
image.plot(1:n, 1:n, chi.hat, main = "estimated chi")
image.plot(1:n, 1:n, ec.hat, main = "estimated EC")

# plot the estimated extremal coefficients for a few counties
these.counties <- sample(x = 159, size = 10, replace = FALSE)

p.1 <- map.ga.ggplot(Y = ec.hat[, 1], 
                     main = paste("Extremal Coefficients for", 
                                  subregion[these.counties[1]], "county"),
                     fill.legend = "EC")

p.2 <- map.ga.ggplot(Y = ec.hat[, 2], 
                     main = paste("Extremal Coefficients for", 
                                  subregion[these.counties[2]], "county"),
                     fill.legend = "EC")

################################################################################
#### Sort and smooth:
#### I don't know that we need to sort, but it makes the plot of the ECs look 
#### better
################################################################################
# sort the counties so it makes sense geographically

# sort the sites by longitude
sites.order.long <- order(cents[, 1])
cents.order.long <- cents[sites.order.for, ]
image.plot(1:n, 1:n, rdist(cents.order.long))

# # more complicated sorting - doesn't appear to improve plot at all
# # dade (#41 the northwest county) is the reference
# sites.first <- 41  # Dade
# sites.last  <- 20  # Camden
# sites.order.for <- order(d[sites.first, ])
# sites.order.rev <- order(d[sites.last, ])
# sites.order <- rep(NA, n)
# sites.order[1:71] <- sites.order.for[1:71]
# sites.order[(n-70):n] <- sites.order.rev[71:1]
# 
# idx <- 72
# for (i in 72:n) {
#   if (!(sites.order.for[i] %in% sites.order)) {
#     sites.order[idx] <- sites.order.for[i] 
#     idx <- idx + 1
#   }
# }
# length(unique(sites.order))
# cents.order <- cents[sites.order, ]
# image.plot(1:n, 1:n, rdist(cents.order))

# Do smoothing
L <- 12

out       <- get.factors.EC(ec.hat,L=L,s=cents)
B.est     <- out$est
alphahat  <- out$alpha
ec.smooth <- out$EC.smooth
ec.est    <- make.EC(B.est, alphahat)

print(out$pct)

# Plot some of the results
p.B1 <- map.ga.ggplot(Y = B.est[, 1], main = "Basis function 1",
                      fill.legend = "Basis value")
p.B2 <- map.ga.ggplot(Y = B.est[, 2], main = "Basis function 2",
                      fill.legend = "Basis value")
p.B3 <- map.ga.ggplot(Y = B.est[, 3], main = "Basis function 3",
                      fill.legend = "Basis value")
p.B4 <- map.ga.ggplot(Y = B.est[, 4], main = "Basis function 4",
                      fill.legend = "Basis value")

grid.arrange(p.B1, p.B2, p.B3, p.B4, ncol = 2, widths = c(1.2, 1.5), 
             main = "Basis functions 1 -- 4")

p.B5 <- map.ga.ggplot(Y = B.est[, 5], main = "Basis function 5",
                      fill.legend = "Basis value")
p.B6 <- map.ga.ggplot(Y = B.est[, 6], main = "Basis function 6",
                      fill.legend = "Basis value")
p.B7 <- map.ga.ggplot(Y = B.est[, 7], main = "Basis function 7",
                      fill.legend = "Basis value")
p.B8 <- map.ga.ggplot(Y = B.est[, 8], main = "Basis function 8",
                      fill.legend = "Basis value")

grid.arrange(p.B5, p.B6, p.B7, p.B8, ncol = 2, widths = c(1.2, 1.5), 
             main = "Basis functions 5 -- 8")

p.B9 <- map.ga.ggplot(Y = B.est[, 9], main = "Basis function 9",
                      fill.legend = "Basis value")
p.B10 <- map.ga.ggplot(Y = B.est[, 10], main = "Basis function 10",
                       fill.legend = "Basis value")
p.B11 <- map.ga.ggplot(Y = B.est[, 11], main = "Basis function 11",
                       fill.legend = "Basis value")
p.B12 <- map.ga.ggplot(Y = B.est[, 12], main = "Basis function 12",
                       fill.legend = "Basis value")

grid.arrange(p.B9, p.B10, p.B11, p.B12, ncol = 2, widths = c(1.2, 1.5), 
             main = "Basis functions 9 -- 12")