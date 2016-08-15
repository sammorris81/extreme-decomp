# How about timing for calculating determinant from chol(Sigma)

ns <- 600
s <- cbind(runif(ns), runif(ns))
d <- rdist(s)
bw <- 0.2
Sigma <- exp(-d / bw)
Sigma.chol <- chol(Sigma)
Qb <- chol2inv(Sigma.chol)

logdet.inv.chol <- function(Sigma.chol) {
  logdet <- -2 * sum(log(diag(Sigma.chol)))
  return(logdet)
}

microbenchmark(logdet.inv.chol(Sigma.chol), logdet(Qb))

logdet(Qb)
logdet.inv.chol(Sigma.chol)

logdet(Sigma)
logdet(chol = Sigma.chol)
logdet(Qb)
-logdet(chol = Sigma.chol)
