# load_data.R is the code I used to load the data.  It stores the data in the fire_data workspace.

#   Y[i,j] = # acres burned in year i and county j.

#   year gives the years, county gives the county names.   I don't have lat/long of the counties, you might have to look that up or use the adjacency matrix (ADJ) to do the smoothing.  Map.ga is a function to map the data, so

#   map.ga(colMeans(Y))

# plots the average for each county.

rm(list=ls())
setwd("S:\\Documents\\PCAX\\Fire")

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




