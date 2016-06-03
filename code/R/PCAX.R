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
  convergence.outer <- rep(FALSE, n)
  iter        <- 0
  while (sum(!convergence.outer) > 0) {
    iter  <- iter + 1
    prev  <- B

    convergence.inner <- rep(FALSE, n)  # storage for optim convergence
    iter.inner <- 0
    maxit <- 1000
    cat("  Starting initial convergence \n")
    while (sum(!convergence.inner) > 0) {
      iter.inner <- iter.inner + 1
      if (iter.inner %% 5 == 0) {
        maxit <- maxit + 500
      }
      B.star <- B^(1 / alpha)
      for (i in 1:n) {
        if (!convergence.inner[i]) {
          fit <- optim(B[i, ], fn = SSE, gr = SSE.grad, Y = EC[i, ],
                       B2 = B, B.star = B.star, alpha = alpha,
                       lower = rep(1.0e-7, L), upper = rep(0.9999999, L),
                       method = "L-BFGS-B", control = list(maxit = maxit))

          # B[i, ] <- abs(fit$par) / sum(abs(fit$par))
          if (fit$convergence == 0) {
            convergence.inner[i] <- TRUE
            B[i, ] <- abs(fit$par) / sum(abs(fit$par))  # change if converged
            B.star[i, ] <- B[i, ]^(1 / alpha)
          } else if (fit$convergence != 0) {
            convergence.inner[i] <- FALSE
          }
        }
        if (i %% 200 == 0) {
          cat("    Finished site ", i, " of ", n, " during inner iter ",
              iter.inner, " of outer iter ", iter, " \n", sep = "")
        }
      }
    }
    cat("  End initial convergence \n")

    cat("  Start convergence check \n")
    convergence.outer <- rep(FALSE, n)
    for (i in 1:n) {  # double check that everything has converged
      fit <- optim(B[i, ], fn = SSE, gr = SSE.grad, Y = EC[i, ],
                   B2 = B, B.star = B.star, alpha = alpha,
                   lower = rep(1.0e-7, L), upper = rep(0.9999999, L),
                   method = "L-BFGS-B", control = list(maxit = 1000))

      # B[i, ] <- abs(fit$par) / sum(abs(fit$par))
      if (fit$convergence == 0) {
        convergence.outer[i] <- TRUE
        B[i, ] <- abs(fit$par) / sum(abs(fit$par))
        B.star[i, ] <- B[i, ]^(1 / alpha)
      } else if (fit$convergence != 0) {
        convergence.outer[i] <- FALSE
      }
      if (i %% 200 == 0) {
        cat("    Finished site ", i, " of ", n, " during outer iteration ",
            iter, " \n", sep = "")
      }
    }
    cat("  End convergence check \n")

    if (iter == 1) {
      Delta_B   <- mean((prev - B)^2)
      Delta_val <- sum((EC-make.EC(B,alpha))^2,na.rm=TRUE)
    } else {
      Delta_B   <- c(Delta_B, mean((prev - B)^2))
      Delta_val <- c(Delta_val, sum((EC - make.EC(B, alpha))^2, na.rm = TRUE))
    }

    if(verbose & (iter %% 5 == 0)){
      cat("  Done with iteration ", iter, ", still need to converge at ",
          sum(!convergence.outer), " sites \n", sep = "")
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
                 convergence = convergence.outer, seconds = tock-tick)

  return(output)
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
  tick <- proc.time()[3]

  n <- ncol(EC)
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
  dw2 <- as.matrix(rdist(s, knots))^2
  dw2[dw2 < 1e-6] <- 0
  dknots2 <- rdist(knots)^2
  diag(dknots2) <- 0

  if (is.null(init.rho)) {
    rho <- max(dknots2[dknots2 != 0]) * 0.2
    print(rho)
  } else {
    rho <- init.rho
  }

  w <- getW(rho = rho, dw2 = dw2)

  # ESTIMATION

  fit <- optim(rho, fn = SSE.rhoalpha, # gr = SSE.grad,
               Y = EC, dw2 = dw2,
               alpha = alpha, lower = 1e-2, upper = 0.5 * sqrt(max(dw2)),
               method = "L-BFGS-B")
  rho <- fit$par
  if (fit$convergence != 0) {
    cat(" Warning, optim returned convergence code", fit$convergence, "\n")
    cat(" Message: ", fit$message, "\n")
  }

  tock <- proc.time()[3]
  output <- list(rho = rho, alpha = alpha, EC.smooth = ECs, dw2 = dw2,
                 seconds = tock-tick)

  return(output)
}

# ####################################################
# # SIMPLE EXAMPLE
# ####################################################
#
# if(FALSE){
#
#  library(splines)
#  library(fields)
#
#  n     <- 100
#  L     <- 5
#  alpha <- 0.3
#
#  #set.seed(0820)
#
#  #Define the truth
#
#   B.true  <- bs(1:n,df=L,intercept=TRUE)
#   B.true  <- sweep(B.true,1,rowSums(B.true),"/")
#   tot     <- colSums(B.true)
#   B.true  <- B.true[,order(-tot)]
#   EC.true <- make.EC(B.true,alpha)
#
#  # The estimate (not generated in a realistic way)
#
#   junk    <- matrix(rnorm(n^2),n,n)
#   EC.hat  <- EC.true+0.01*t(junk)%*%junk
#   EC.hat  <- ifelse(EC.hat<1,1,EC.hat)
#   EC.hat  <- ifelse(EC.hat>2,2,EC.hat)
#
#   diag(EC.hat) <- NA
#
#  # Estimation
#
#
#   out       <- get.factors.EC(EC.hat,L=L,s=1:n,bw=5)
#   B.est     <- out$est
#   alphahat  <- out$alpha
#   EC.smooth <- out$EC.smooth
#   EC.est    <- make.EC(B.est, alphahat)
#
#   print(out$pct)
#
#  # Plot the results
#
#   par(mfrow=c(3,2))
#   matplot(B.true,type="l",main="True B")
#   matplot(B.est,type="l",main="Estimated B")
#   image.plot(1:n,1:n,EC.true,main="True EC")
#   image.plot(1:n,1:n,EC.hat,main="Initial EC estimate (theta-hat)")
#   image.plot(1:n,1:n,EC.smooth,main="Smoothed EC (theta-tilde)")
#   image.plot(1:n,1:n,EC.est,main="Final EC estimate")
#
#
# }
