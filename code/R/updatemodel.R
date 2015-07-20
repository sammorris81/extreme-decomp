updateMu <- function(FAC, a, b, t, mu, sigma, xi, cur.ll, 
                     acc.a, att.a, mh.a,
                     acc.b, att.b, mh.b) {
  
  att.a <- att.a + 1
  att.b <- att.b + 1
  
  can.mu <- mu
  for (t in 1:nt) {
    can.mu <- sum(FAC)
  }
  
}

updateSigma <- function() {
  
  

}

updateXi <- function() {

  
}

updateLogs <- function(logs, y, mu, sigma, xi, theta, alpha, u, cur.ll, FAC, 
                       cuts, mh) {
  for (t in 1:nt) { 
    v   <- xi * exp(-sigma) * (y[, t] - mu) + 1
    v   <- v^(-1 / (alpha * xi))      
    ccc <- logd(theta[, t], v)
    
    for (l in 1:nF) {
      W           <- FAC[, l]^(1 / alpha)
      level[l, t] <- get.level(logs[l, t], cuts)
      MH1         <- mh[level[l, t]]
      can.logs    <- rnorm(1, logs[l, t], MH1)
      MH2         <- mh[get.level(can.logs, cuts)]
      
      can.theta <- theta[, t] + W * (exp(can.logs) - exp(logs[l, t]))
      can.ccc   <- logd(cantheta, v)
      R <- can.ccc - ccc +
           h1(can.logs, u[l, t], alpha, log = TRUE) - 
           h1(logs[l, t], u[l, t], alpha, log = TRUE) +
           dnorm(logs[l, t], canlogs, MH2, log = TRUE) -
           dnorm(canlogs, logs[l, t], MH1, log = TRUE)
      if (!is.na(exp(R))) { if (runif(1) < exp(R)) {
        logs[l, t] <- canlogs
        ccc        <- canccc
        theta[, t] <- cantheta
      }}
    }
    cur.ll[, t] <- loglike(y[, t], mu[, t], sigma, xi, theta[, t], alpha)
  }
  
  results <- list(logs = logs, theta = theta, cur.ll = cur.ll)
}

updateAlpha <- function(y, alpha, logs, FAC, theta, mu, sigma, xi, u, cur.ll,
                        acc, att, mh) {
  att <- att + 1
  
  alpha.star     <- qnorm(alpha)
  can.alpha.star <- rnorm(1, alpha.star, mh)
  can.alpha      <- pnorm(can.alpha.star)
  can.theta      <- make.theta(FAC, logs, can.alpha)
  can.ll         <- cur.ll  
  
  for (t in 1:nt) {
    can.ll[, t] <- loglike(y[, t], mu[, t], sigma, xi, can.theta[, t], 
                           can.alpha)
  }
  R <- sum(canll - curll) +
       sum(h1(logs, u, canalpha)) -
       sum(h1(logs, u, alpha)) +
       dbeta(canalpha, pri.alpha.a, pri.alpha.b, log = TRUE) -
       dbeta(alpha, pri.alpha.a, pri.alpha.b, log = TRUE)
  
  if (!is.na(exp(R))) { if (runif(1) < exp(R)) {
    alpha  <- can.alpha
    cur.ll <- can.ll
    theta  <- can.theta
    acc    <- acc + 1
  }}           
  
  results <- list(alpha = alpha, theta = theta, cur.ll = cur.ll, 
                  acc = acc, att = att)
  return()
}