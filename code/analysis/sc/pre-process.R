# get counties
sc.counties <- read.csv(file = "sc_fires.csv", header = TRUE, 
                        stringsAsFactors = FALSE)$County
sc.counties <- sc.counties[sc.counties != ""]
sc.counties <- sc.counties[-length(sc.counties)]

# when reading in the data, 
# column 1 is the county name
# column 2 is whether it's count or acres burned
# data are fiscal years year 1 = 1946-47
#                       year T = 2014-15
sc.data.csv <- read.csv(file = "sc_fires.csv", header = FALSE, skip = 1, 
                        stringsAsFactors = FALSE)[1:92, 3:71]
acre.idx <- seq(2, 92, by = 2)  # odds are counts, evens are acres
sc.data <- matrix(NA, 46, 69)
for (t in 1:ncol(sc.data)) {
  # extract the acreage burned for each year
  sc.data[, t] <- as.numeric(gsub(",","", sc.data.csv[acre.idx, t]))
}
rownames(sc.data) <- sc.counties
colnames(sc.data) <- 1946:2014
Y <- t(sc.data)
save(Y, sc.counties, file = "sc_preprocess/fire_data.RData")

################################################################################
#### Get the centroids
################################################################################
rm(list=ls())
library(maps)
library(maptools)
# gpclibPermit()
sc <- map("county", "south carolina", fill = TRUE, col = "transparent",
               plot = FALSE)
range(sc$x, na.rm = TRUE)
range(sc$y, na.rm = TRUE)
sc$names
IDs <- sapply(strsplit(sc$names, ","), function(x) x[2])

# maptools
sc_sp <- map2SpatialPolygons(sc, IDs = IDs,
                                  proj4string=CRS("+proj=longlat +datum=WGS84"))
cents <- coordinates(sc_sp)
save(cents, sc, IDs, sc_sp, file="sc_preprocess/sc_map.RData")

