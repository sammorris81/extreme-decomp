######################################################################
#
# Function to estimate the B functions:
#
# Inputs:
#
#  EC.smooth := n x n matrix of spatially smoothed pairwise ECs
#  s         := locations
#  knots     := knot locations
#  alpha     := positive stable parameter alpha
#  init.rho  := inital value
#
# Outputs
#
#  rho       := estimated value of rho
#  alpha     := estimated alpha
#  EC.smooth := smoothed version of EC
#
######################################################################
get.rho <- function(EC.smooth, alpha.hat, s = NULL, knots = NULL,
                    init.rho = NULL, verbose = TRUE){
  require(fields)
  tick <- proc.time()[3]
  n <- ncol(EC.smooth)

  if (is.null(knots)) {  # place knots at sites
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

  # Find the best rho for alpha.
  fit <- optim(rho, fn = SSE.rhoalpha, Y = EC.smooth, dw2 = dw2,
               alpha = alpha.hat, #lower = 1e-2, upper = 0.5 * sqrt(max(dw2)),
               method = "BFGS")
  rho <- fit$par
  if (fit$convergence != 0) {
    cat(" Warning, optim returned convergence code", fit$convergence, "\n")
    cat(" Message: ", fit$message, "\n")
  }

  tock <- proc.time()[3]
  output <- list(rho = rho, dw2 = dw2, seconds = tock-tick)

  return(output)
}
