
######################################################################
# 
# Function to estimate the B functions:
#
# Inputs:
#
#  EC      := n x n matrix of estimated pairwise extremal coefficients
#  L       := number of B functions to be estimated
#  alpha   := positive stable parameter alpha
#  init.B  := inital value
#  iters   := number of iterations in the optimization algorithm
#
# Outputs
#
#  est       := estimated value of B
#  alpha     := estimated alpha
#  EC.smooth := smoothed version of EC
# 
######################################################################

get.rho.alpha <- function(EC, s = NULL, knots = NULL, bw = NULL, alpha = NULL,
                          init.rho = NULL, verbose = TRUE){
  require(fields)
  tick   <- proc.time()[3]
  
  n      <- ncol(EC)
  if (is.null(s)) {s <- 1:n}
  if (is.null(knots)) {knots <- s}  # place knots at sites
  if (is.null(bw)) {bw <- 2 * min(dist(s))}
  
  # SMOOTHING
  
  EC  <- Ksmooth(EC, s, bw)
  ECs <- EC
  if (is.null(alpha)) {alpha <- log2(mean(diag(EC)))}
  diag(EC) <- NA
  
  # INITIAL VALUES
  if (is.null(knots)) {
    knots <- s
  }
  dw2             <- as.matrix(rdist(s, knots))^2
  dw2[dw2 < 1e-6] <- 0
  
  if (is.null(init.rho)) {
    rho <- quantile(dw2, probs = 0.15)
  }
  
  w <- getW(rho = rho, dw2 = dw2)
  
  # ESTIMATION
  
  fit <- optim(rho, fn = SSE, gr = SSE.grad, Y = EC, dw2 = dw2, 
               alpha = alpha, lower = 0, upper = 0.5 * sqrt(max(dw2)), 
               method = "L-BFGS-B", control = list(maxit = 1000))
  rho <- fit$par
    
  output <- list(est = rho, alpha = alpha, EC.smooth = ECs,
                 seconds = tock-tick)
  
  return(output)
}

# SSE
SSE <- function(rho, dw2, Y, alpha) {
  w <- getW(rho = rho, dw2 = dw2)
  w <- w^(1 / alpha)
  
  n <- ncol(Y)
  EC <- matrix(NA, n, n)
  for (i in 1:n) {
    for (j in (i+1):n) {
      EC[i, j] <- EC[j, i] <- sum(colSums(w[c(i, j), ])^alpha)
    }
  } 
  sse <- sum((Y - EC)^2, na.rm = TRUE)
}

# hack for finding the gradient
SSE.grad <- function(rho, dw2, Y, alpha){
  
  delta <- 0.001
  epsilon <- SSE(rho = rho + delta, dw2 = dw2, Y = Y, alpha = alpha) - 
             SSE(rho = rho, dw2 = dw2, Y = Y, alpha = alpha)
  
  grad <- epsilon / delta
  
  return(grad)
}


make.EC  <- function(B, alpha){
  Ba    <- B^(1 / alpha) 
  EC    <- NULL
  for(j in 1:nrow(B)){
    BB <- sweep(Ba, 2, Ba[j, ], "+")
    EC <- cbind(EC, rowSums(BB^alpha))
  }
  
  return(EC)
}

# Performs kernel smoothing of the extremal coefficient matrix.
Ksmooth <- function(ECmat, s = NULL, bw = NULL){
  
  n           <- nrow(ECmat)
  diag(ECmat) <- 0
  E1          <- ifelse(ECmat == 0, 0, 1)
  if (is.null(s)) {s <- 1:n}
  if (is.null(bw)) {bw <- 2 * min(dist(s))}
  
  
  d2       <- as.matrix(dist(s) / bw)^2
  W        <- exp(-d2)
  diag(W)  <- 0
  
  num      <- W %*% ECmat %*% W
  den      <- W %*% E1 %*% W
  
  ECsmooth <- num / den
  
  return(ECsmooth)
}

getW <- function(rho, dw2) {
  w <- stdW(makeW(dw2 = dw2, rho = rho))
}


# get the kernel weighting
makeW <- function(dw2, rho) {
  w <- exp(-0.5 * dw2 / (rho^2))

  return(w)
}

# standardize the kernel weights
stdW <- function(x, single = FALSE) {
  if (single) {x <- x / sum(x)}
  if (!single) {x <- sweep(x, 1, rowSums(x), "/")}
  return(x)
}