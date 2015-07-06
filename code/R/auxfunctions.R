smooth.ec <- function(ec, s, diag.keep = FALSE) {
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