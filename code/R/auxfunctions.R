smooth.ec.loess <- function(ec, s, diag.keep = FALSE) {
  ns <- nrow(s)
  
  # function to smooth the extremal coefficients
  lower <- lower.tri(ec)
  if (diag.keep) {
    diag(lower) <- TRUE
    rows <- 1:ns
  } else {
    rows <- 1:(ns - 1)
  }
  
  # gets the bottom half of the triangle
  these <- which(lower != 0)
  
  # create the vectors of locations
  s11 <- s12 <- s21 <- s22 <- rep(NA, length(these))
  this <- 0
  for (i in rows) {
    if (diag.keep) {cols <- i:ns} else {cols <- (i + 1):ns}
    for (j in cols) {
      this <- this + 1 
      s11[this] <- s[i, 1]
      s12[this] <- s[i, 2]
      s21[this] <- s[j, 1]
      s22[this] <- s[j, 2]
    }
  }
  
  ls <- loess(as.vector(ec[these]) ~ s11 + s12 + s21 + s22)
  
  smoothed <- matrix(1, ns, ns)
  this <- 0
  for (i in rows) {
    if (diag.keep) { cols <- i:ns } else { cols <- (i + 1):ns }
    for (j in cols) {
      this <- this + 1
      smoothed[i, j] <- smoothed[j, i] <- fitted(ls)[this]
    }
  }
  
  results <- list(actual = ec, s = s, smoothed = smoothed)
  return(results)
}

map.ga.ggplot <- function(Y, main = "", fill.legend = "", 
                          low = NULL, high = NULL) {
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
  
  georgia <- map("county", "georgia", fill = TRUE, col = "transparent",
                 plot = FALSE)
  subregion <- sapply(strsplit(georgia$names, ","), function(x) x[2])
  county_map <- map_data(map = "county", region = "georgia")
  
  Y     <- data.frame(Y, subregion)
  Y.med <- median(Y$Y, na.rm = TRUE)
  map   <- merge(county_map, Y, all.x=TRUE)
  
  p <- ggplot(map, aes(x=long, y=lat, group=group, fill=Y))
  p <- p + geom_polygon(colour="grey", aes(fill=Y))
  p <- p + expand_limits(x = map$long, y = map$lat) 
  p <- p + coord_map("polyconic") 
  p <- p + labs(title = main, fill = fill.legend) 
  p <- p + scale_fill_gradient2(low = "dodgerblue4", high = "firebrick4", mid = "#ffffff", 
                                    midpoint = Y.med)
  p <- p + theme_clean()
  
  return(p)
}