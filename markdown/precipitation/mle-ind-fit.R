set.seed(2000)
temp.mtx <- matrix(rnorm(30), 10, 3)
temp.array <- array(rnorm(120), dim = c(10, 3, 4))

loop <- rep(0, dim(temp.array)[3])
for (i in 1:length(loop)) {
  loop[i] <- sum(temp.mtx * temp.array[, , i])
}

app <- apply(as.vector(temp.mtx) * temp.array, 3, sum)


library(microbenchmark)
microbenchmark(
  for (i in 1:length(loop)) {
    loop[i] <- sum(temp.mtx * temp.array[, , i])
  },
  apply(as.vector(temp.mtx) * temp.array, 3, sum)
)

temp <- matrix(rnorm(16), 4, 4)
temp <- abs(temp)
solve(temp)
solve(temp + diag(0.0001, 4))


data(PORTw)
# get MLE estimates for log(scale)
fit0 <- fevd(TMX1, PORTw,
             location.fun = ~MTMAX + MTMIN,
             scale.fun = ~MTMAX + MTMIN,
             type = "GEV", units = "deg C", use.phi = TRUE)
fit0



ll.ind <- function(beta, X, y) {
  nt <- dim(X)[2]
  np <- dim(X)[3]

  beta1 <- beta[1:np]
  beta2 <- beta[(np + 1):(2 * np)]
  xi    <- tail(beta, 1)
  ll <- 0

  for (t in 1:nt) {
    mu <- X[, t, ] %*% beta1
    sig <- exp(X[, t, ] %*% beta2)
    if (any(sig == 0)) {
      return(Inf)
    }
    ll <- ll + sum(
      dgev(y[, t], loc = mu, scale = sig, shape = xi, log = TRUE),
      na.rm = TRUE)
  }

  return(-ll)
}
