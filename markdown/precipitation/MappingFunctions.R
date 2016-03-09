
library(ggplot2)
library(maps)

# Include functions to generate latitude/longitude pairs to cover the US
# These functions were written by Neal Grantham
generate_grid <- function(x = 100, y = 100){
  # purpose: generate grid of points over the US at which to make predictions
  # returns: two column matrix, longitude and latitude of all grid points
  library(maps)  # for map.where
  corners <- enclose_USA()
  # create grid
  grid <- expand.grid(seq(corners[1], corners[2], length = x), 
                      seq(corners[3], corners[4], length = y))
  # retain only points that fall over US soil
  inUSA <- !is.na(map.where("usa", x = grid[, 1], y = grid[, 2]))
  grid <- as.matrix(grid[inUSA, ])
  # name the columns and rows
  colnames(grid) <- c("lon", "lat")
  rownames(grid) <- 1:nrow(grid)
  return(grid)
}

enclose_USA <- function() {
  # purpose: geographic coordinate limits within the continental US
  # returns: vector of corners of USA
  max.lat <- 49.384472    # most northern
  min.lat <- 24.520833    # most southern
  max.lon <- -66.947028   # most eastern
  min.lon <- -124.733056  # most western
  corners <- c(min.lon, max.lon, min.lat, max.lat)
}

# Mapping functions
map.points <- function (lat, lon, data, 
                        color_low="white",color_high="darkred",color_na=gray(0.9),zeroiswhite=FALSE,
                        xlim=NULL, ylim=NULL, zlim=NULL,
                        mainTitle=NULL, legendTitle="") {

  # Created by Susheela Singh!
  
  # Store the base data of the underlying map
  baseData <- map_data("state")
  
  # Combine the data into a dataframe
  dfMap <- as.data.frame(cbind(lon, lat, data))
  colnames(dfMap) <- c("lon", "lat", "Value")
  
  # Set limits for x, y, z if not specified as parameters
  if (is.null(xlim)) { xlim <- range( lon,na.rm=TRUE) }
  if (is.null(ylim)) { ylim <- range( lat,na.rm=TRUE) }
  if (is.null(zlim)) { zlim <- range(data,na.rm=TRUE) }
    
  # Create the plot
  p <- ggplot(dfMap, aes(x=lon, y=lat)) + theme_bw()
  p <- p + theme(plot.title = element_text(size = rel(1.5)))
  p <- p + geom_point(aes(colour = Value))
  p <- p + geom_polygon(data=baseData, aes(x=long, y=lat, group=group), 
                        colour="black", fill="white", alpha=0) 
  p <- p + coord_fixed(ratio=1.1, xlim=xlim, ylim=ylim)
  if(zeroiswhite){
    p <- p + scale_colour_gradient2(low=color_low, 
                                    high=color_high,
                                    na.value=color_na,
                                    limits=zlim,
                                    name=legendTitle) 
  }
  if(!zeroiswhite){
    p <- p + scale_colour_gradient(low=color_low, 
                                   high=color_high,
                                   na.value=color_na,
                                   limits=zlim,
                                   name=legendTitle) 
  }

  return(p)  
}


map.heatmap <- function (lat, lon, data, 
                         color_low="white",color_high="darkred",color_na=gray(0.9),zeroiswhite=FALSE,
                         xlim=NULL, ylim=NULL, zlim=NULL,
                         mainTitle="", legendTitle="") {
  
  # Created by Susheela Singh!

  # Store the base data of the underlying map
  baseData <- map_data("state")


  
  # Combine the data into a dataframe
  dfMap <- as.data.frame(cbind(lon, lat, data))
  colnames(dfMap) <- c("lon", "lat", "Value")
    
  # Set limits for x, y, z if not specified as parameters
  if (is.null(xlim)) { xlim <- range( lon,na.rm=TRUE) }
  if (is.null(ylim)) { ylim <- range( lat,na.rm=TRUE) }
  if (is.null(zlim)) { zlim <- range(data,na.rm=TRUE) }

  # Create the plot
  p <- ggplot(dfMap, aes(x=lon, y=lat, fill=Value)) + theme_bw()
  p <- p + geom_tile()
  p <- p + geom_polygon(data=baseData, aes(x=long, y=lat, group=group), 
                        colour="black", fill="white", alpha=0) 
  p <- p + labs(title=paste(mainTitle,"\n",sep=""), x="", y="")
  p <- p + theme(plot.title = element_text(size = rel(1.5))) 
  p <- p + coord_fixed(ratio=1.1, xlim=xlim, ylim=ylim)

  if(zeroiswhite){
    p <- p + scale_fill_gradient2(low=color_low, 
                                  high=color_high,
                                  na.value=color_na,
                                  limits=zlim,
                                  name=legendTitle) 
  }
  if(!zeroiswhite){
    p <- p + scale_fill_gradient(low=color_low, 
                                 high=color_high,
                                 na.value=color_na,
                                 limits=zlim,
                                 name=legendTitle) 
  }
  
  return(p)  
}

if(TRUE){

 # Test the mapping functions
  mapGrid <- generate_grid(200,200) 
  y       <- sin(2*pi*mapGrid[,1]/10)+
             cos(2*pi*mapGrid[,2]/10)
  y <- rank(y)/length(y)

  these<-runif(length(y))<0.3

  map.points(lat=mapGrid[these, 2], lon=mapGrid[these, 1], y[these],
             color_low="pink",color_high="black",
             xlim=c(-85, -65), zlim=c(0, 1), mainTitle="Eastern United States")

  map.heatmap(lat=mapGrid[, 2], lon=mapGrid[, 1], y,color_low="white",color_high="black",
              xlim=c(-85, -65), zlim=c(0, 1), mainTitle="Eastern United States")
}



