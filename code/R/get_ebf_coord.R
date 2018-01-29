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
SSE.B <- function(B1, B2, B.star, Y, alpha, lambda = 1000){

  BB  <- B1^(1 / alpha)
  B2  <- B.star
  EC  <- sweepC2plus(X = B2, y = BB)
  EC  <- rowSumsC(EC^alpha)

  # penalty term is to make sure that the bases sum to 1
  sse <- sum((Y - EC)^2, na.rm = TRUE) + lambda * (sum(B1) - 1)^2

  return(sse)
}

SSE.B.grad <- function(B1, B.star, Y, alpha, lambda = 1000, exclude = 1){

  B1.star <- B1^(1 / alpha)
  B2   <- B.star

  BB <- sweepC2plus(X = B2, y = B1.star)
  EC0  <- rowSumsC(BB^alpha)

  EC1  <- BB^(alpha - 1)
  EC1  <- sweepC2times(EC1, B1.star / B1)
  EC1  <- sweepC1times(EC1, Y - EC0)

  grad <- -2 * colSums(EC1, na.rm = TRUE) +
    2 * lambda * (sum(B1) - 1)

  return(grad)
}

get.factors.EC.coord <- function(EC.smooth, alpha.hat = NULL, L = 5, s = NULL,
                           bw = NULL, init.B = NULL, iters = 10,
                           verbose = TRUE){

  tick   <- proc.time()[3]

  n      <- ncol(EC.smooth)
  if (is.null(s)) {s <- 1:n }
  if (is.null(bw)) {bw <- 2 * min(dist(s)) }

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
      B.star <- B^(1 / alpha.hat)
      for (i in 1:n) {
        if (!convergence.inner[i]) {
          fit <- optim(B[i, ], fn = SSE.B, gr = SSE.B.grad, Y = EC.smooth[i, ],
                       B.star = B.star, alpha = alpha.hat,
                       lower = rep(1.0e-7, L), upper = rep(0.9999999, L),
                       method = "L-BFGS-B", control = list(maxit = maxit))

          # B[i, ] <- abs(fit$par) / sum(abs(fit$par))
          if (fit$convergence == 0) {
            convergence.inner[i] <- TRUE
            B[i, ] <- abs(fit$par) / sum(abs(fit$par))  # change if converged
            B.star[i, ] <- B[i, ]^(1 / alpha.hat)
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
      fit <- optim(B[i, ], fn = SSE.B, gr = SSE.B.grad, Y = EC.smooth[i, ],
                   B.star = B.star, alpha = alpha.hat,
                   lower = rep(1.0e-7, L), upper = rep(0.9999999, L),
                   method = "L-BFGS-B", control = list(maxit = 1000))

      # B[i, ] <- abs(fit$par) / sum(abs(fit$par))
      if (fit$convergence == 0) {
        convergence.outer[i] <- TRUE
        B[i, ] <- abs(fit$par) / sum(abs(fit$par))
        B.star[i, ] <- B[i, ]^(1 / alpha.hat)
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
      Delta_val <- sum((EC.smooth - make.EC(B, alpha.hat))^2, na.rm = TRUE)
    } else {
      Delta_B   <- c(Delta_B, mean((prev - B)^2))
      Delta_val <- c(Delta_val, sum((EC.smooth - make.EC(B, alpha.hat))^2,
                                    na.rm = TRUE))
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

  output <- list(est = B, pct = pct, seconds = tock - tick)

  return(output)
}
