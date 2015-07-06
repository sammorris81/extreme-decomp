# load_data.R is the code I used to load the data.
# It stores the data in the fire_data workspace.

# Y[i,j] = # acres burned in year i and county j.

# year gives the years, county gives the county names.
# I don't have lat/long of the counties, you might have to look that up or use
# the adjacency matrix (ADJ) to do the smoothing.

# Map.ga is a function to map the data, so

#   map.ga(colMeans(Y))

# plots the average for each county.

rm(list=ls())
# setwd("S:\\Documents\\PCAX\\Fire")

source("adj.R")
load("gaCntyFires.RData")


y     <- matrix(0,876,159)
year  <- sort(rep(1957:2029,12))
month <- rep(1:12,73)

county <- colnames(cnty.fire)
n      <- length(county)
for(j in 1:n){
   y[,j]<-cnty.fire[,j]
}


y     <- y[year<=2013,]
month <- month[year<=2013]
year  <- year[year<=2013]

map.ga<-function(Y,main="",low=NULL,high=NULL){
   library(maps)

   if(is.null(low)){ low  <- min(Y,na.rm=TRUE)}
   if(is.null(high)){high <- max(Y,na.rm=TRUE)}

   library(maps)
   library(fields)

   colorfx<-function(z){
       z[is.na(z)]<-0
   z}
   Z<-colorfx((Y-low)/(high-low))
   map("county","georgia",fill=TRUE,col=gray(1-Z))


   image.plot(low+(high-low)*diag(2),legend.only=TRUE,
              col=gray(1-colorfx(seq(0,1,.01))),
              smallplot=c(.7,.74,0.6,0.85))

   title(main)
}


YEAR <- unique(year)
m    <- length(YEAR)
Y    <- matrix(0,m,n)

for(j in 1:m){for(k in 1:n){
   Y[j,k]<-sum(y[year==YEAR[j],k])
}}

rm(year,y,month,cnty.fire,j,k)


adj<-NULL
for(i in 1:n){for(j in 1:n){if(i<j){
   if(ADJ[i,j]==1){adj<-rbind(adj,c(i,j))}
}}}

adj_list<-apply(ADJ==1,1,which)

year<-YEAR

rm(adj,i,j,m,n,YEAR)

save.image("fire_data.RData")

map.ga(colMeans(Y))

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

save(chi, file="chi.RData")
load("chi.RData")
image.plot(1:ncol(Y), 1:ncol(Y), chi, main = "estimated chi")
image.plot(1:ncol(Y), 1:ncol(Y), ADJ, main = "adjacency")

# sort the adjacency matrix
order <- rep(NA, ncol(Y))
order[1] <- 1
for (i in 1:length(order)) {
  this <- order[i]
  neighbors <- which(ADJ[this, ] == 1)
  for (k in 1:length(neighbors)) {
    if (!(neighbors[k] %in% unique(order))) {
      this <- min(which(is.na(order)))
      order[this] <- neighbors[k]
    }
  }
}

chi.sort <- matrix(NA, ncol(Y), ncol(Y))
for (i in 1:ncol(Y)) {
  for (j in i:ncol(Y)) {
    pair.i <- order[i]
    pair.j <- order[j]
    chi.sort[i, j] <- chi.sort[j, i] <- chi[pair.i, pair.j]
  }
}
par(mfrow=c(1, 2))
chi.sort <- ifelse(chi.sort < 0, 0, chi.sort)
# image.plot(1:ncol(Y), 1:ncol(Y), chi, main = "estimated chi")
image.plot(1:ncol(Y), 1:ncol(Y), chi.sort, main = "estimated chi")
image.plot(1:ncol(Y), 1:ncol(Y), 2 - chi.sort, main = "estimated EC")

# lag-1 dependence
nyears <- nrow(Y)
plot(Y[-nyears, 1], Y[-1, 1], main = "year total: lag 1, site 1",
     xlab = bquote(y[t]), ylab = bquote(y[t + 1]))
chiplot(cbind(Y[-nyears, 1], Y[-1, 1]), which = 1, ylim1 = c(0, 1),
        main1 = bquote(paste(chi, "-plot: lag 1, site1")))
plot(Y[-nyears, 2], Y[-1, 2], main = "year total: lag 1, site 2",
     xlab = bquote(y[t]), ylab = bquote(y[t + 1]))

nmonths <- nrow(y)
par(mfrow = c(2, 2))
plot(y[-nmonths, 1], y[-1, 1], main = "monthly total: lag 1, site 1",
     xlab = bquote(y[t]), ylab = bquote(y[t + 1]))
chiplot(cbind(y[-nmonths, 1], y[-1, 1]), which = 1, ylim1 = c(0, 1),
        main1 = bquote(paste(chi,"-plot: lag 1, site 1")))
plot(y[-nmonths, 2], y[-1, 2], main = "monthly total: lag 1, site 2",
     xlab = bquote(y[t]), ylab = bquote(y[t + 1]))
chiplot(cbind(y[-nmonths, 2], y[-1, 2]), which = 1, ylim1 = c(0, 1),
        main1 = bquote(paste(chi,"-plot: lag 1, site 2")))

plot(y[-nmonths, 3], y[-1, 3], main = "monthly total: lag 1, site 3",
     xlab = bquote(y[t]), ylab = bquote(y[t + 1]))
chiplot(cbind(y[-nmonths, 3], y[-1, 3]), which = 1, ylim1 = c(0, 1),
        main1 = bquote(paste(chi,"-plot: lag 1, site 3")))
plot(y[-nmonths, 4], y[-1, 4], main = "monthly total: lag 1, site 4",
     xlab = bquote(y[t]), ylab = bquote(y[t + 1]))
chiplot(cbind(y[-nmonths, 4], y[-1, 4]), which = 1, ylim1 = c(0, 1),
        main1 = bquote(paste(chi,"-plot: lag 1, site 4")))

plot(y[-nmonths, 5], y[-1, 5], main = "monthly total: lag 1, site 5",
     xlab = bquote(y[t]), ylab = bquote(y[t + 1]))
chiplot(cbind(y[-nmonths, 5], y[-1, 5]), which = 1, ylim1 = c(0, 1),
        main1 = bquote(paste(chi,"-plot: lag 1, site 5")))
plot(y[-nmonths, 6], y[-1, 6], main = "monthly total: lag 1, site 6",
     xlab = bquote(y[t]), ylab = bquote(y[t + 1]))
chiplot(cbind(y[-nmonths, 6], y[-1, 6]), which = 1, ylim1 = c(0, 1),
        main1 = bquote(paste(chi,"-plot: lag 1, site 6")))

all.sites <- cbind(y[-nmonths, 1], y[-1, 1])
for (i in 2:ncol(y)) {
  all.sites <- rbind(all.sites, cbind(y[-nmonths, i], y[-1, i]))
}
par(mfrow = c(1, 2))
plot(all.sites[, 1], all.sites[, 2], main = "monthly total: lag 1, all sites",
     xlab = bquote(y[t]), ylab = bquote(y[t + 1]))
chiplot(all.sites, which = 1, ylim1 = c(0, 1),
        main1 = bquote(paste(chi, "-plot: lag 1, all sites")))

