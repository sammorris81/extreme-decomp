
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
#  iters   := minimum number of iterations in the optimization algorithm.
#             the algorithm will continue until convergence is reached at 
#             all sites.
#
# Outputs
#
#  est       := estimated value of B
#  alpha     := estimated alpha
#  EC.smooth := smoothed version of EC
# 
######################################################################

get.factors.EC <- function(EC, L = 5, s = NULL, bw = NULL, alpha = NULL,
                           init.B = NULL, iters = 10, verbose = TRUE){

 tick   <- proc.time()[3]

 n      <- ncol(EC)
 if (is.null(s)) {s <- 1:n }
 if (is.null(bw)) {bw <- 2 * min(dist(s)) }

 # SMOOTHING

  EC  <- Ksmooth(EC, s, bw)  # run a kernel smoother on the pairwise estimates
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
  
  # keep running the algorithm until we get convergence at all sites.
  convergence <- 0
  iter        <- 0
  while (convergence > 0 | iter < iters) {
    convergence <- 0  # storage for optim convergence
    iter  <- iter + 1
    prev  <- B
    maxit <- ifelse(iter == iters | iter > 100, 100, 2 * iter + 2)
    
    for (i in 1:n) {
      fit <- optim(B[i, ], fn = SSE, gr = SSE.grad, Y = EC[i, ], B2 = B, 
                   alpha = alpha, lower = rep(0, L), upper = rep(1, L), 
                   method = "L-BFGS-B", control = list(maxit = maxit))
      B[i, ] <- abs(fit$par) / sum(abs(fit$par))
      if (fit$convergence != 0) {
        convergence <- convergence + 1
      }
    }
    
    Delta_B[iter]   <- mean((prev-B)^2)
    Delta_val[iter] <- sum((EC-make.EC(B,alpha))^2,na.rm=TRUE)
    if(verbose & (iter %% 5 == 0)){
      cat("  Done with iteration ", iter, ", still need to converge at ", 
          convergence, " sites \n", sep = "")
    }
  }
  
  if (iter > iters) {
    cat("  estimating basis functions took ", iter, "iterations. \n")
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
                 convergence = convergence, seconds = tock-tick)
  
  return(output)
}

# SSE for row of Y-EC
SSE <- function(B1, B2, Y, alpha, lambda = 1000){

  BB  <- B1^(1 / alpha)
  B2  <- B2^(1 / alpha)
  EC  <- sweep(B2, 2, BB, "+")
  EC  <- rowSums(EC^alpha)
  
  # penalty term is to make sure that the bases sum to 1
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

######################################################################
# 
# Function to estimate the B functions:
#
# Inputs:
#
#  EC       := n x n matrix of estimated pairwise extremal coefficients
#  s        := locations
#  knots    := knot locations
#  alpha    := positive stable parameter alpha
#  init.rho := inital value
#
# Outputs
#
#  rho       := estimated value of rho
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
  
  fit <- optim(rho, fn = SSE.rhoalpha, # gr = SSE.grad, 
               Y = EC, dw2 = dw2, 
               alpha = alpha, lower = 0, upper = 0.5 * sqrt(max(dw2)), 
               method = "L-BFGS-B")
  rho <- fit$par
  if (fit$convergence != 0) {
    cat(" Warning, optim returned convergence code", fit$convergence, "\n")
  }
  
  tock <- proc.time()[3]
  output <- list(rho = rho, alpha = alpha, EC.smooth = ECs, dw2 = dw2,
                 seconds = tock-tick)
  
  return(output)
}

# SSE
SSE.rhoalpha <- function(rho, dw2, Y, alpha) {
  w <- getW(rho = rho, dw2 = dw2)
  w <- w^(1 / alpha)
  
  n <- ncol(Y)
  EC <- matrix(NA, n, n)
  for (i in 1:(n - 1)) {
    for (j in (i+1):n) {
      EC[i, j] <- EC[j, i] <- sum(colSums(w[c(i, j), ])^alpha)
    }
  } 
  
  sse <- sum((Y - EC)^2, na.rm = TRUE)
  
  return(sse)
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


####################################################
# SIMPLE EXAMPLE
####################################################

if(FALSE){

 library(splines)
 library(fields)

 n     <- 100
 L     <- 5
 alpha <- 0.3

 #set.seed(0820)

 #Define the truth
 
  B.true  <- bs(1:n,df=L,intercept=TRUE)
  B.true  <- sweep(B.true,1,rowSums(B.true),"/")
  tot     <- colSums(B.true)
  B.true  <- B.true[,order(-tot)]
  EC.true <- make.EC(B.true,alpha)

 # The estimate (not generated in a realistic way)

  junk    <- matrix(rnorm(n^2),n,n)
  EC.hat  <- EC.true+0.01*t(junk)%*%junk
  EC.hat  <- ifelse(EC.hat<1,1,EC.hat)
  EC.hat  <- ifelse(EC.hat>2,2,EC.hat)

  diag(EC.hat) <- NA

 # Estimation

  
  out       <- get.factors.EC(EC.hat,L=L,s=1:n,bw=5)
  B.est     <- out$est
  alphahat  <- out$alpha
  EC.smooth <- out$EC.smooth
  EC.est    <- make.EC(B.est, alphahat)

  print(out$pct)

 # Plot the results

  par(mfrow=c(3,2))
  matplot(B.true,type="l",main="True B")
  matplot(B.est,type="l",main="Estimated B")
  image.plot(1:n,1:n,EC.true,main="True EC")
  image.plot(1:n,1:n,EC.hat,main="Initial EC estimate (theta-hat)")
  image.plot(1:n,1:n,EC.smooth,main="Smoothed EC (theta-tilde)")
  image.plot(1:n,1:n,EC.est,main="Final EC estimate")

 
}
