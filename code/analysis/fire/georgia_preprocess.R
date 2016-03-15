################################################################################
#### Get updated data from GA
################################################################################
setwd("~/repos-git/extreme-decomp/code/analysis/fire/county_fire/")
ncounties <- 159
nyears    <- 58  # 1957 - 2014

for (i in 1:ncounties) {
  if (i < 10) {
    filename <- paste("http://weather.gfc.state.ga.us/FireData/CT00", i, 
                      "MONCYA.TXT", sep = "")
  } else if (i < 100) {
    filename <- paste("http://weather.gfc.state.ga.us/FireData/CT0", i, 
                      "MONCYA.TXT", sep = "")
  } else {
    filename <- paste("http://weather.gfc.state.ga.us/FireData/CT", i, 
                      "MONCYA.TXT", sep = "")
  }
  cmd <- paste("curl -O", filename)
  system(cmd)
}

county <- rep(NA, ncounties)
Y <- matrix(NA, nyears, ncounties)  # to maintain compatibility with existing code
for (i in 1:ncounties) {
  if (i < 10) {
    filename <- paste("CT00", i, "MONCYA.TXT", sep = "")
  } else if (i < 100) {
    filename <- paste("CT0", i, "MONCYA.TXT", sep = "")
  } else {
    filename <- paste("CT", i, "MONCYA.TXT", sep = "")
  }
  
  county[i] <- read.table(file = filename, nrows = 1, 
                          stringsAsFactors = FALSE)$V1
  Y[, i]    <- read.table(file = filename, skip = 1, header = TRUE)[1:nyears, 2]
}

colnames(Y) <- county
rownames(Y) <- seq(1957, 2014)

save(Y, file = "fire_data.RData")

################################################################################
#### Yearly amounts from each county
################################################################################
rm(list=ls())

# Y[i,j] = # acres burned in year i and county j.
# y: monthly data
# year gives the years, county gives the county names.
# Map.ga is a function to map the data, so to plot the average for each county
#   map.ga(colMeans(Y))

source("adj.R")
load("gaCntyFires.RData")

y     <- matrix(0, 876, 159)
year  <- sort(rep(1957:2029, 12))
month <- rep(1:12, 73)

county <- colnames(cnty.fire)
n      <- length(county)
for(j in 1:n){
  y[, j]<-cnty.fire[, j]
}

y     <- y[year <= 2013, ]
month <- month[year <= 2013]
year  <- year[year <= 2013]

map.ga <- function(Y, main = "", low = NULL, high = NULL){
  library(maps)
  
  if(is.null(low)){low   <- min(Y, na.rm = TRUE)}
  if(is.null(high)){high <- max(Y, na.rm = TRUE)}
  
  library(maps)
  library(fields)
  
  colorfx <- function(z){
    z[is.na(z)] <- 0
    z
  }
  Z <- colorfx((Y - low) / (high - low))
  map("county", "georgia", fill = TRUE, col = gray(1 - Z))
  
  image.plot(low + (high - low) * diag(2), legend.only = TRUE,
             col = gray(1 - colorfx(seq(0, 1, 0.01))),
             smallplot = c(0.7, 0.74, 0.6, 0.85))
  
  title(main)
}

YEAR <- unique(year)
m    <- length(YEAR)
Y    <- matrix(0, m, n)

for (j in 1:m) {
  for (k in 1:n) {
    Y[j, k] <- sum(y[year == YEAR[j], k])
  }
}

rm(year, month, cnty.fire, j, k)

adj <- NULL
for (i in 1:n) { for (j in 1:n) {
  if (i < j) { if (ADJ[i, j] == 1) {
    adj <- rbind(adj, c(i, j))
  } }
} }

adj_list <- apply(ADJ == 1, 1, which)

year <- YEAR

rm(adj, i, j, m, n, YEAR)

load("./county_fire/fire_data.RData")
Y <- Y[10:nrow(Y), ]  # years 1957 - 1964 are likely not reliable
save.image("georgia_preprocess/fire_data.RData")

################################################################################
#### Get the centroids
################################################################################
rm(list=ls())
library(maps)
library(maptools)
# gpclibPermit()
georgia <- map("county", "georgia", fill = TRUE, col = "transparent",
               plot = FALSE)
range(georgia$x, na.rm = TRUE)
range(georgia$y, na.rm = TRUE)
georgia$names
IDs <- sapply(strsplit(georgia$names, ","), function(x) x[2])

# maptools
georgia_sp <- map2SpatialPolygons(georgia, IDs = IDs,
                                  proj4string=CRS("+proj=longlat +datum=WGS84"))
cents <- coordinates(georgia_sp)
save(cents, georgia, IDs, georgia_sp, 
     file="georgia_preprocess/georgia_map.RData")

################################################################################
#### Get the estimate for chi to be used in the extremal coefficient
################################################################################
rm(list=ls())
load(file = "georgia_preprocess/fire_data.RData")
library(evd)
chi <- matrix(NA, ncol(Y), ncol(Y))
for (i in 1:ncol(Y)) {
  for (j in i:ncol(Y)) {
    chi.ij <- mean(chiplot(Y[, c(i, j)], which = 1, ask = FALSE)$chi[95:100, 2])
    chi[i, j] <- chi[j, i] <- chi.ij
    if (j %% 50 == 0) {
      print(paste("j:", j))
    }
  }
  if (i %% 10 == 0) {
    print(paste("i:", i))
  }
}
save(chi, file="georgia_preprocess/chi.RData")
