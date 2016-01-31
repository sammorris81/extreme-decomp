
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
                          init.rho = NULL, iters = 10, verbose = TRUE){
  
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
  
  if(is.null(init.B)){
    B <- 1 / matrix(1:L, n, L, byrow = TRUE)
  }
  if(!is.null(init.B)){
    B <- init.B
  }
  B  <- sweep(B, 1, rowSums(B), "/")
  
  # ESTIMATION
  
  Delta_B   <- rep(NA, iters)
  Delta_val <- rep(NA, iters)
  
  for (iter in 1:iters) {
    prev  <- B
    maxit <- ifelse(iter == iters | iter > 100, 100, 2 * iter + 2)
    for (i in 1:n) {
      fit <- optim(B[i, ], fn = SSE, gr = SSE.grad, Y = EC[i, ], B2 = B, 
                   alpha = alpha, lower = rep(0, L), upper = rep(1, L), 
                   method = "L-BFGS-B", control = list(maxit = maxit))
      B[i, ] <- abs(fit$par) / sum(abs(fit$par))
    }
    Delta_B[iter]   <- mean((prev-B)^2)
    Delta_val[iter] <- sum((EC-make.EC(B,alpha))^2,na.rm=TRUE)
    
    if(verbose){cat("    Done with iteration", iter, "of", iters, "\n")}
  }
  
  # REORDER THE COLUMNS
  
  if (L == 1) {
    B   <- matrix(B, n, L)
    pct <- 1
  } else {
    B   <- B[, order(-colSums(B))]
    pct <- colSums(B) / sum(B)
  }
  tock   <- proc.time()[3]
  
  output <- list(est = B, pct = pct, alpha = alpha, EC.smooth = ECs,
                 Delta.B = Delta_B, Delta.val = Delta_val,
                 seconds = tock-tick)
  
  return(output)
}

# SSE for row of Y-EC
SSE <- function(B1, B2, Y, alpha, lambda = 1000){
  
  BB  <- B1^(1 / alpha)
  B2  <- B2^(1 / alpha)
  EC  <- sweep(B2, 2, BB, "+")
  EC  <- rowSums(EC^alpha)
  sse <- sum((Y - EC)^2, na.rm = TRUE) + lambda * (sum(B1) - 1)^2
  
  return(sse)
}

SSE.grad <- function(B1, B2, Y, alpha, lambda = 1000){
  
  BB   <- B1^(1 / alpha)
  B2   <- B2^(1 / alpha)
  
  BB   <- sweep(B2, 2, BB, "+")
  EC0  <- rowSums(BB^alpha)
  
  EC1  <- BB^(alpha - 1)
  EC1  <- sweep(EC1, 2, B1^(1 / alpha - 1), "*")
  EC1  <- sweep(EC1, 1, Y - EC0, "*")
  
  grad <- -2 * colSums(EC1, na.rm = TRUE) +
    2 * lambda * (sum(B1) - 1)
  
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
  if (is.null(a.cutoff)) {
    a.cutoff <- max(sqrt(dw2))
  }
  w <- exp(-0.5 * dw2 / (rho^2))

  return(w)
}

# standardize the kernel weights
stdW <- function(x, single = FALSE) {
  if (single) {x <- x / sum(x)}
  if (!single) {x <- sweep(x, 1, rowSums(x), "/")}
  return(x)
}