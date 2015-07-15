# setup
library(maps)
library(maptools)
library(fields)
source(file = "../code/analysis/fire/adj.R", chdir = TRUE)
source(file = "../code/R/auxfunctions.R", chdir = TRUE)
load(file = "../code/analysis/fire/gaCntyFires.RData")
load(file = "../code/analysis/fire/fire_data.RData")
load(file = "../code/analysis/fire/chi.RData")

georgia <- map("county", "georgia", fill = TRUE, col = "transparent",
               plot = FALSE)
range(georgia$x, na.rm = TRUE)
range(georgia$y, na.rm = TRUE)
georgia$names
IDs <- sapply(strsplit(georgia$names, ","), function(x) x[2])
# maptools
georgia_sp <- map2SpatialPolygons(georgia, IDs = IDs,
                                  proj4string=CRS("+proj=longlat +datum=WGS84"))
sapply(slot(georgia_sp, "polygons"), function(x) length(slot(x, "Polygons")))
plot(georgia_sp, col="grey", axes = TRUE)
cents <- coordinates(georgia_sp)
ns <- nrow(cents)
points(cents)

d <- rdist(cents)
diag(d) <- 0

# library(evd)
# chi <- matrix(NA, ncol(Y), ncol(Y))
# for (i in 1:ncol(Y)) {
#   for (j in i:ncol(Y)) {
#     chi.ij <- mean(chiplot(Y[, c(i, j)], which = 1, ask = FALSE)$chi[95:100, 2])
#     chi[i, j] <- chi[j, i] <- chi.ij
#     if (j %% 50 == 0) {
#       print(paste("j:", j))
#     }
#   }
#   if (i %% 10 == 0) {
#     print(paste("i:", i))
#   }
# }

# chi[i, j] is the chi statistic between sites i and j
chi <- ifelse(chi <= 0, 0, chi)
image.plot(1:ncol(Y), 1:ncol(Y), chi, main = "estimated chi")
image.plot(1:ncol(Y), 1:ncol(Y), 2 - chi, main = "estimated EC")
map.ga(colMeans(Y))
map.ga(2 - chi[1, ], main = "EC with site 1")
points(cents[1, ])

# maptools
georgia_sp <- map2SpatialPolygons(georgia, IDs = IDs,
                                  proj4string=CRS("+proj=longlat +datum=WGS84"))
georgia_df <- data.frame(value = chi[, 1], row.names = IDs)
georgia_sp_df <- SpatialPolygonsDataFrame(Sr = georgia_sp, data = georgia_df)
georgia_cood_df <- SpatialPointsDataFrame(coords = coordinates(georgia_sp), 
                                          data = georgia_df, match.ID = FALSE,
                                          proj4string=CRS("+proj=longlat +datum=WGS84"))
spplot(obj = georgia_sp_df, scales = list(draw = TRUE))

these.counties <- sample(x = 159, size = 10, replace = FALSE)

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

extcoef <- 2 - chi[, these.counties[1]]
subregion <- sapply(strsplit(georgia$names, ","), function(x) x[2])
extcoef <- data.frame(extcoef, subregion)
county_map <- map_data(map = "county", region = "georgia")
extcoef_map <- merge(county_map, extcoef, all.x=TRUE)

rm(p)
p <- ggplot(extcoef_map, aes(x=long, y=lat, group=group, fill=extcoef))
p <- p + geom_polygon(colour="grey", aes(fill=extcoef))
p <- p + expand_limits(x = extcoef_map$long, y = extcoef_map$lat) 
p <- p + coord_map("polyconic") 
p <- p + labs(title = paste("Extremal Coefficients for", 
                            subregion[these.counties[1]], "county"), 
              fill = "Extremal Coefficient") 
p <- p + scale_fill_gradient2(low = "firebrick4", high = "dodgerblue4", mid = "#ffffff", 
                              midpoint = 1.5)
p <- p + theme_clean()
p

extcoef <- 2 - chi[, these.counties[2]]
subregion <- sapply(strsplit(georgia$names, ","), function(x) x[2])
extcoef <- data.frame(extcoef, subregion)
county_map <- map_data(map = "county", region = "georgia")
extcoef_map <- merge(county_map, extcoef, all.x=TRUE)

rm(p)
p <- ggplot(extcoef_map, aes(x=long, y=lat, group=group, fill=extcoef))
p <- p + geom_polygon(colour="grey", aes(fill=extcoef))
p <- p + expand_limits(x = extcoef_map$long, y = extcoef_map$lat) 
p <- p + coord_map("polyconic") 
p <- p + labs(title = paste("Extremal Coefficients for", 
                            subregion[these.counties[2]], "county"), 
              fill = "Extremal Coefficient") 
p <- p + scale_fill_gradient2(low = "firebrick4", high = "dodgerblue4", mid = "#ffffff", 
                              midpoint = 1.5)
p <- p + theme_clean()
p

# extremal coefficient
# still need to redo this loess smoother to smooth on s1, s2
ec <- 2 - chi
ec.lower <- lower.tri(ec)
diag(ec.lower) <- TRUE
ec.lower <- ec * ec.lower  # gets the bottom half of the triangle
these <- which(ec.lower != 0)
ls <- loess(as.vector(ec[these]) ~ as.vector(d[these]))
plot(as.vector(d[these]), fitted(ls))
smoothed <- matrix(NA, ns, ns)
this <- 0
for (i in 1:ns) {
  for (j in i:ns) {
    this <- this + 1
    smoothed[i, j] <- smoothed[j, i] <- fitted(ls)[this]
  }
}
plot(d[, 1], smoothed[, 1])

# do smoothing
temp <- smooth.ec(ec, cents)

# make plots
library(gridExtra)
library(rapport)

# county 1
this.county <- 1
extcoef <- temp$smoothed[, this.county]
county <- tocamel(subregion[this.county], upper = TRUE)
df <- data.frame(dist = d[, this.county], smooth = temp$smoothed[, this.county])
p1a <- ggplot(df, aes(x = dist, y = smooth)) + geom_point(size = 3)
p1a <- p1a + geom_smooth(method = loess)
p1a <- p1a + coord_fixed(ratio=6)
p1a <- p1a + xlab("distance") + ylab(expression(theta))

subregion <- sapply(strsplit(georgia$names, ","), function(x) x[2])
extcoef <- data.frame(extcoef, subregion)
county_map <- map_data(map = "county", region = "georgia")
extcoef_map <- merge(county_map, extcoef, all.x=TRUE)

p1b <- ggplot(extcoef_map, aes(x=long, y=lat, group=group, fill=extcoef))
p1b <- p1b + geom_polygon(colour="grey", aes(fill=extcoef))
p1b <- p1b + expand_limits(x = extcoef_map$long, y = extcoef_map$lat) 
p1b <- p1b + coord_map("polyconic") 
p1b <- p1b + labs(fill = "Extremal Coefficient") 
p1b <- p1b + scale_fill_gradient2(low = "firebrick4", high = "dodgerblue4", mid = "#ffffff", 
                              midpoint = 1.5, limits = c(1, 2))
p1b <- p1b + theme_clean()

grid.arrange(p1a, p1b, ncol = 2, widths = c(1.2, 1.5), 
             main = textGrob(paste("Extremal coefficients for ", county, 
                                   " county", sep = ""), 
                             vjust = 5, gp = gpar(fontface = "bold", cex = 1.5)))

# county 6
this.county <- 6
extcoef <- temp$smoothed[, this.county]
county <- tocamel(subregion[this.county], upper = TRUE)
df <- data.frame(dist = d[, this.county], smooth = temp$smoothed[, this.county])
p1a <- ggplot(df, aes(x = dist, y = smooth)) + geom_point(size = 3)
p1a <- p1a + geom_smooth(method = loess)
p1a <- p1a + coord_fixed(ratio=6)
p1a <- p1a + xlab("distance") + ylab(expression(theta))

subregion <- sapply(strsplit(georgia$names, ","), function(x) x[2])
extcoef <- data.frame(extcoef, subregion)
county_map <- map_data(map = "county", region = "georgia")
extcoef_map <- merge(county_map, extcoef, all.x=TRUE)

p1b <- ggplot(extcoef_map, aes(x=long, y=lat, group=group, fill=extcoef))
p1b <- p1b + geom_polygon(colour="grey", aes(fill=extcoef))
p1b <- p1b + expand_limits(x = extcoef_map$long, y = extcoef_map$lat) 
p1b <- p1b + coord_map("polyconic") 
p1b <- p1b + labs(fill = "Extremal Coefficient") 
p1b <- p1b + scale_fill_gradient2(low = "firebrick4", high = "dodgerblue4", mid = "#ffffff", 
                                  midpoint = 1.5, limits = c(1, 2))
p1b <- p1b + theme_clean()

grid.arrange(p1a, p1b, ncol = 2, widths = c(1.2, 1.5), 
             main = textGrob(paste("Extremal coefficients for ", county, 
                                   " county", sep = ""), 
                             vjust = 5, gp = gpar(fontface = "bold", cex = 1.5)))

# county 15
this.county <- 15
extcoef <- temp$smoothed[, this.county]
county <- tocamel(subregion[this.county], upper = TRUE)
df <- data.frame(dist = d[, this.county], smooth = temp$smoothed[, this.county])
p1a <- ggplot(df, aes(x = dist, y = smooth)) + geom_point(size = 3)
p1a <- p1a + geom_smooth(method = loess)
p1a <- p1a + coord_fixed(ratio=6)
p1a <- p1a + xlab("distance") + ylab(expression(theta))

subregion <- sapply(strsplit(georgia$names, ","), function(x) x[2])
extcoef <- data.frame(extcoef, subregion)
county_map <- map_data(map = "county", region = "georgia")
extcoef_map <- merge(county_map, extcoef, all.x=TRUE)

p1b <- ggplot(extcoef_map, aes(x=long, y=lat, group=group, fill=extcoef))
p1b <- p1b + geom_polygon(colour="grey", aes(fill=extcoef))
p1b <- p1b + expand_limits(x = extcoef_map$long, y = extcoef_map$lat) 
p1b <- p1b + coord_map("polyconic") 
p1b <- p1b + labs(fill = "Extremal Coefficient") 
p1b <- p1b + scale_fill_gradient2(low = "firebrick4", high = "dodgerblue4", mid = "#ffffff", 
                                  midpoint = 1.5, limits = c(1, 2))
p1b <- p1b + theme_clean()

grid.arrange(p1a, p1b, ncol = 2, widths = c(1.2, 1.5), 
             main = textGrob(paste("Extremal coefficients for ", county, 
                                   " county", sep = ""), 
                             vjust = 5, gp = gpar(fontface = "bold", cex = 1.5)))

# county 20
this.county <- 20
extcoef <- temp$smoothed[, this.county]
county <- tocamel(subregion[this.county], upper = TRUE)
df <- data.frame(dist = d[, this.county], smooth = temp$smoothed[, this.county])
p1a <- ggplot(df, aes(x = dist, y = smooth)) + geom_point(size = 3)
p1a <- p1a + geom_smooth(method = loess)
p1a <- p1a + coord_fixed(ratio=6)
p1a <- p1a + xlab("distance") + ylab(expression(theta))

subregion <- sapply(strsplit(georgia$names, ","), function(x) x[2])
extcoef <- data.frame(extcoef, subregion)
county_map <- map_data(map = "county", region = "georgia")
extcoef_map <- merge(county_map, extcoef, all.x=TRUE)

p1b <- ggplot(extcoef_map, aes(x=long, y=lat, group=group, fill=extcoef))
p1b <- p1b + geom_polygon(colour="grey", aes(fill=extcoef))
p1b <- p1b + expand_limits(x = extcoef_map$long, y = extcoef_map$lat) 
p1b <- p1b + coord_map("polyconic") 
p1b <- p1b + labs(fill = "Extremal Coefficient") 
p1b <- p1b + scale_fill_gradient2(low = "firebrick4", high = "dodgerblue4", mid = "#ffffff", 
                                  midpoint = 1.5, limits = c(1, 2))
p1b <- p1b + theme_clean()

grid.arrange(p1a, p1b, ncol = 2, widths = c(1.2, 1.5), 
             main = textGrob(paste("Extremal coefficients for ", county, 
                                   " county", sep = ""), 
                             vjust = 5, gp = gpar(fontface = "bold", cex = 1.5)))

# county 24
this.county <- 24
extcoef <- temp$smoothed[, this.county]
county <- tocamel(subregion[this.county], upper = TRUE)
df <- data.frame(dist = d[, this.county], smooth = temp$smoothed[, this.county])
p1a <- ggplot(df, aes(x = dist, y = smooth)) + geom_point(size = 3)
p1a <- p1a + geom_smooth(method = loess)
p1a <- p1a + coord_fixed(ratio=6)
p1a <- p1a + xlab("distance") + ylab(expression(theta))

subregion <- sapply(strsplit(georgia$names, ","), function(x) x[2])
extcoef <- data.frame(extcoef, subregion)
county_map <- map_data(map = "county", region = "georgia")
extcoef_map <- merge(county_map, extcoef, all.x=TRUE)

p1b <- ggplot(extcoef_map, aes(x=long, y=lat, group=group, fill=extcoef))
p1b <- p1b + geom_polygon(colour="grey", aes(fill=extcoef))
p1b <- p1b + expand_limits(x = extcoef_map$long, y = extcoef_map$lat) 
p1b <- p1b + coord_map("polyconic") 
p1b <- p1b + labs(fill = "Extremal Coefficient") 
p1b <- p1b + scale_fill_gradient2(low = "firebrick4", high = "dodgerblue4", mid = "#ffffff", 
                                  midpoint = 1.5, limits = c(1, 2))
p1b <- p1b + theme_clean()

grid.arrange(p1a, p1b, ncol = 2, widths = c(1.2, 1.5), 
             main = textGrob(paste("Extremal coefficients for ", county, 
                                   " county", sep = ""), 
                             vjust = 5, gp = gpar(fontface = "bold", cex = 1.5)))

# county 148
this.county <- 148
extcoef <- temp$smoothed[, this.county]
county <- tocamel(subregion[this.county], upper = TRUE)
df <- data.frame(dist = d[, this.county], smooth = temp$smoothed[, this.county])
p1a <- ggplot(df, aes(x = dist, y = smooth)) + geom_point(size = 3)
p1a <- p1a + geom_smooth(method = loess)
p1a <- p1a + coord_fixed(ratio=6)
p1a <- p1a + xlab("distance") + ylab(expression(theta))

subregion <- sapply(strsplit(georgia$names, ","), function(x) x[2])
extcoef <- data.frame(extcoef, subregion)
county_map <- map_data(map = "county", region = "georgia")
extcoef_map <- merge(county_map, extcoef, all.x=TRUE)

p1b <- ggplot(extcoef_map, aes(x=long, y=lat, group=group, fill=extcoef))
p1b <- p1b + geom_polygon(colour="grey", aes(fill=extcoef))
p1b <- p1b + expand_limits(x = extcoef_map$long, y = extcoef_map$lat) 
p1b <- p1b + coord_map("polyconic") 
p1b <- p1b + labs(fill = "Extremal Coefficient") 
p1b <- p1b + scale_fill_gradient2(low = "firebrick4", high = "dodgerblue4", mid = "#ffffff", 
                                  midpoint = 1.5, limits = c(1, 2))
p1b <- p1b + theme_clean()

grid.arrange(p1a, p1b, ncol = 2, widths = c(1.2, 1.5), 
             main = textGrob(paste("Extremal coefficients for ", county, 
                                   " county", sep = ""), 
                             vjust = 5, gp = gpar(fontface = "bold", cex = 1.5)))