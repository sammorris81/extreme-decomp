rm(list=ls())

library(splines)
library(maps)
library(maptools)
library(fields)
source(file = "../code/analysis/fire/adj.R", chdir = TRUE)
source(file = "../code/R/auxfunctions.R", chdir = TRUE)
load(file = "../code/analysis/fire/gaCntyFires.RData")
load(file = "../code/analysis/fire/fire_data.RData")

# get the Georgia map and coordinates
load(file = "../code/analysis/fire/georgia_centroids.RData")
d <- rdist(cents)
diag(d) <- 0
n <- nrow(cents)

# chi[i, j] is the chi statistic between sites i and j
load(file = "../code/analysis/fire/chi.RData")
chi.hat <- ifelse(chi <= 0, 0, chi)
ec.hat  <- 2 - chi.hat
image.plot(1:n, 1:n, chi.hat, main = "estimated chi")
image.plot(1:n, 1:n, ec.hat, main = "estimated EC")

# sort the counties so it makes sense geographically
# dade (#41 the northwest county) is the reference
sites.first <- 41  # Dade
sites.last  <- 20  # Camden
sites.order.for <- order(d[sites.first, ])
sites.order.rev <- order(d[sites.last, ])
sites.order <- rep(NA, n)
sites.order[1:71] <- sites.order.for[1:71]
sites.order[(n-70):n] <- sites.order.rev[71:1]

idx <- 72
for (i in 72:n) {
  if (!(sites.order.for[i] %in% sites.order)) {
    sites.order[idx] <- sites.order.for[i] 
    idx <- idx + 1
  }
}
length(unique(sites.order))
cents.order <- cents[sites.order, ]
image.plot(1:n, 1:n, rdist(cents.order))

sites.order.long <- order(cents[, 1])
cents.order.long <- cents[sites.order.for, ]
image.plot(1:n, 1:n, rdist(cents.order.long))




L <- 10

out       <- get.factors.EC(ec.hat,L=L,s=cents)
B.est     <- out$est
alphahat  <- out$alpha
ec.smooth <- out$EC.smooth
ec.est    <- make.EC(B.est, alphahat)

print(out$pct)

# Plot the results

par(mfrow=c(2,2))
matplot(B.est,type="l",main="Estimated B")
# image.plot(1:n,1:n,EC.true,main="True EC")
image.plot(1:n,1:n,ec.hat,main="Initial EC estimate (theta-hat)")
image.plot(1:n,1:n,ec.smooth,main="Smoothed EC (theta-tilde)")
image.plot(1:n,1:n,ec.est,main="Final EC estimate")

# map the B functions over GA
map.ga(B.est[, 1])
map.ga(B.est[, 2])
map.ga(B.est[, 3])
map.ga(B.est[, 4])
map.ga(B.est[, 5])

library(ggplot2)
theme_clean <- function(base_size = 12) {
  require(grid)
  theme_grey(base_size) %+replace%
    theme(
      axis.title      =   element_blank(),
      axis.text       =   element_blank(),
      panel.background    =   element_blank(),
      panel.grid      =   element_blank(),
      axis.ticks.length   =   unit(0,"cm"),
      axis.ticks.margin   =   unit(0,"cm"),
      panel.margin    =   unit(0,"lines"),
      plot.margin     =   unit(c(0,0,0,0),"lines"),
      complete = TRUE
    )
}

subregion <- sapply(strsplit(georgia$names, ","), function(x) x[2])
county_map <- map_data(map = "county", region = "georgia")

basis.1 <- B.est[, 1]
basis.1 <- data.frame(basis.1, subregion)
extcoef_map.1 <- merge(county_map, basis.1, all.x=TRUE)

p.1 <- ggplot(extcoef_map.1, aes(x=long, y=lat, group=group, fill=basis.1))
p.1 <- p.1 + geom_polygon(colour="grey", aes(fill=basis.1))
p.1 <- p.1 + expand_limits(x = extcoef_map.1$long, y = extcoef_map.1$lat) 
p.1 <- p.1 + coord_map("polyconic") 
p.1 <- p.1 + labs(title = "1st basis function", fill = "B") 
p.1 <- p.1 + scale_fill_gradient2(low = "dodgerblue4", high = "firebrick4", mid = "#ffffff", 
                              midpoint = median(B.est[, 1]))
p.1 <- p.1 + theme_clean()

basis.2 <- B.est[, 2]
basis.2 <- data.frame(basis.2, subregion)
extcoef_map.2 <- merge(county_map, basis.2, all.x=TRUE)

p.2 <- ggplot(extcoef_map.2, aes(x=long, y=lat, group=group, fill=basis.2))
p.2 <- p.2 + geom_polygon(colour="grey", aes(fill=basis.2))
p.2 <- p.2 + expand_limits(x = extcoef_map.2$long, y = extcoef_map.2$lat) 
p.2 <- p.2 + coord_map("polyconic") 
p.2 <- p.2 + labs(title = "2nd basis function", fill = "B") 
p.2 <- p.2 + scale_fill_gradient2(low = "dodgerblue4", high = "firebrick4", mid = "#ffffff", 
                                  midpoint = median(B.est[, 2]))
p.2 <- p.2 + theme_clean()

basis.3 <- B.est[, 3]
basis.3 <- data.frame(basis.3, subregion)
extcoef_map.3 <- merge(county_map, basis.3, all.x=TRUE)

p.3 <- ggplot(extcoef_map.3, aes(x=long, y=lat, group=group, fill=basis.3))
p.3 <- p.3 + geom_polygon(colour="grey", aes(fill=basis.3))
p.3 <- p.3 + expand_limits(x = extcoef_map.3$long, y = extcoef_map.3$lat) 
p.3 <- p.3 + coord_map("polyconic") 
p.3 <- p.3 + labs(title = "3rd basis function", fill = "B") 
p.3 <- p.3 + scale_fill_gradient2(low = "dodgerblue4", high = "firebrick4", mid = "#ffffff", 
                                  midpoint = median(B.est[, 3]))
p.3 <- p.3 + theme_clean()

basis.4 <- B.est[, 4]
basis.4 <- data.frame(basis.4, subregion)
extcoef_map.4 <- merge(county_map, basis.4, all.x=TRUE)

p.4 <- ggplot(extcoef_map.4, aes(x=long, y=lat, group=group, fill=basis.4))
p.4 <- p.4 + geom_polygon(colour="grey", aes(fill=basis.4))
p.4 <- p.4 + expand_limits(x = extcoef_map.4$long, y = extcoef_map.4$lat) 
p.4 <- p.4 + coord_map("polyconic") 
p.4 <- p.4 + labs(title = "4th basis function", fill = "B") 
p.4 <- p.4 + scale_fill_gradient2(low = "dodgerblue4", high = "firebrick4", mid = "#ffffff", 
                                  midpoint = median(B.est[, 4]))
p.4 <- p.4 + theme_clean()

grid.arrange(p.1, p.2, p.3, p.4, nrow = 2,
             ncol = 2, widths = c(2, 2), 
             main = "Basis function for extremal coefficients")


basis.4 <- B.est[, 5]
basis.4 <- data.frame(basis.4, subregion)
extcoef_map.4 <- merge(county_map, basis.4, all.x=TRUE)

p.4 <- ggplot(extcoef_map.4, aes(x=long, y=lat, group=group, fill=basis.4))
p.4 <- p.4 + geom_polygon(colour="grey", aes(fill=basis.4))
p.4 <- p.4 + expand_limits(x = extcoef_map.4$long, y = extcoef_map.4$lat) 
p.4 <- p.4 + coord_map("polyconic") 
p.4 <- p.4 + labs(title = "4th basis function", fill = "B") 
p.4 <- p.4 + scale_fill_gradient2(low = "dodgerblue4", high = "firebrick4", mid = "#ffffff", 
                                  midpoint = median(B.est[, 5]))
p.4 <- p.4 + theme_clean()

rm(p.4)
basis.4 <- B.est[, 6]
basis.4 <- data.frame(basis.4, subregion)
extcoef_map.4 <- merge(county_map, basis.4, all.x=TRUE)

p.4 <- ggplot(extcoef_map.4, aes(x=long, y=lat, group=group, fill=basis.4))
p.4 <- p.4 + geom_polygon(colour="grey", aes(fill=basis.4))
p.4 <- p.4 + expand_limits(x = extcoef_map.4$long, y = extcoef_map.4$lat) 
p.4 <- p.4 + coord_map("polyconic") 
p.4 <- p.4 + labs(title = "4th basis function", fill = "B") 
p.4 <- p.4 + scale_fill_gradient2(low = "dodgerblue4", high = "firebrick4", mid = "#ffffff", 
                                  midpoint = median(B.est[, 6]))
p.4 <- p.4 + theme_clean()

rm(p.4)

map.ga(colMeans(Y))
