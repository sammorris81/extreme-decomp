logpost.beta.sigma <- function(beta, beta.mn, beta.sd, X.sig, y, theta, mu, xi, 
                               thresh, alpha) {
  # loglikelihood
  logsig <- 0
  for(j in 1:p.sig){
    logsig <- logsig + X.sig[, , j] * beta[j]
  }
  
  sigma <- exp(logsig)
  mu_star  <- mu + sigma * ((theta^xi) - 1) / xi
  sig_star <- alpha * sigma * (theta^xi)
  xi_star  <- alpha * xi
  
  tx <- (1 + xi_star * (y - mu_star) / sig_star)^(-1 / xi_star)
  ll <- -tx + (y > thresh) * ((xi_star + 1) * log(tx) - log(sig_star))
  ll <- ifelse(is.na(y), 0, ll)  # maybe this handles the missing data
  ll <- ifelse(is.na(ll), -Inf, ll)
  
  # account for prior 
  ll <- sum(ll) + sum(dnorm(beta, beta.mn, beta.sd, log = TRUE))
  
  return(ll)
}

grad.logsig.beta <- function(beta, beta.mn, beta.sd, X.sig, y, theta, mu, xi,
                             thresh, alpha) {
  
  p.sig <- length(beta)
  grad <- rep(-99999, p)
  
  # loglikelihood
  logsig <- 0
  for(j in 1:p.sig){
    logsig <- logsig + X.sig[, , j] * beta[j]
  }
  
  sigma <- exp(logsig)
  mu_star  <- mu + sigma * ((theta^xi) - 1) / xi
  sig_star <- alpha * sigma * (theta^xi)
  xi_star  <- alpha * xi
  
  xi.y.sig <- xi * y * sigma
  xi.mu.sig <- xi * mu * sigma
  
  # loglikelihood:
  tx.1 <- theta^(-xi) * (xi.y.sig - xi.mu.sig)
  lltx <- (tx.1 + theta^(-xi))^(-1 / xi_star - 1) / xi_star
  for (p in 1:p.sig)
  
}