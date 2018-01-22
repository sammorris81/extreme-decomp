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

udpate_one_row <- function(theta, Bcur, Brest, alpha,
                           maxit = 100, fudge = 0.001){
  b <- optim(log(Bcur) - mean(log(Bcur)), fn = SSE1row, gr = SSE1rowgrad,
             theta = theta, Brest = Brest, alpha = alpha, fudge = fudge,
             control = list(maxit = maxit))$par
  B <- exp(b) / sum(exp(b))
  return(B)
}

SSE1row <- function(b, theta, Brest, alpha) {
  Bcur     <- exp(b) / sum(exp(b))
  p1       <- sweep(Brest^(1 / alpha), 2, Bcur^(1 / alpha), "+")
  thetahat <- rowSums(p1^alpha)
  sse      <- sum((theta - thetahat)^2)
  return(sse)
}

SSEall <- function(alpha, theta, B){
  thetahat <- B2EC_cpp(B^(1 / alpha), alpha)
  sse      <- sum((theta - thetahat)^2, na.rm = TRUE)
  return(sse)
}

B2EC <- function(B, alpha){
  n     <- nrow(B)
  B     <- B^(1 / alpha)
  theta <- 0
  for (l in 1:ncol(B)) {
   Bl    <- matrix(B[, l], n, n, byrow = TRUE) +
            matrix(B[, l], n, n, byrow = FALSE)
   theta <- theta + Bl^alpha
  }
  diag(theta) <- 2^alpha
  return(theta)
}


sse <- function(ab, EC, n, L, alpha, fudge = 0.001) {
  b <- matrix(ab, n, L)
  B <- sweep(exp(b), 1, rowSums(exp(b)), "/")
  SSEall(alpha, EC, B) / (n * L) + fudge * sum(ab^2) / (n * L)
}

sse_grad <- function(ab, EC, n, L, alpha, fudge = 0.001) {
  b <- matrix(ab, n, L)
  B <- sweep(exp(b), 1, rowSums(exp(b)), "/")

  Gb <- b
  for (j in 1:n) {
    Gb[j, ] <- SSE1rowgrad(b[j, ], EC[j, -j], B[-j, ], alpha = alpha)
  }

  G <- as.vector(Gb) + 2 * fudge * ab

  return(G / (n * L))
}

sse_grad <- function(ab, EC, n, L, alpha, fudge = 0.001) {
  b       <- matrix(ab, n, L)
  B       <- sweep(exp(b), 1, rowSums(exp(b)), "/")
  Ba      <- B^(1 / alpha)
  theta   <- B2EC_cpp(Ba, alpha)

  R       <- EC - theta
  diag(R) <- 0
  Gb      <- 0 * b
  for (k in 1:L) {
    Rd     <- GetRd_cpp(Ba = Ba[, k], alpha = alpha, R = R)
    Rd     <- rowSums(Rd)
    Gb     <- Gb + sweep(B, 1, Rd * Ba[, k], "*")
    Gb[,k] <- Gb[, k] - Rd * Ba[, k]
  }
  G <- 2 * as.vector(Gb) + 2 * fudge * ab

  return(G / (n * L))
}

SSE1rowgrad <- function(b, theta, Brest, alpha) {
  B         <- exp(b) / sum(exp(b))
  delta     <- sweep(Brest^(1 / alpha), 2, B^(1 / alpha), "+")
  thetahat  <- rowSums(delta^alpha)
  G         <- diag(B^(1 / alpha - 1)) %*% (diag(B) - outer(B,B))
  thetahat2 <- delta^(alpha - 1)
  q0        <- -2 * colSums(sweep(thetahat2, 1, theta - thetahat, "*"))
  grad      <- q0 %*% G

  return(grad)
}

B2EC <- function(B, alpha){
  n     <- nrow(B)
  B     <- B^(1 / alpha)
  theta <- 0
  for (l in 1:ncol(B)) {
    Bl    <- matrix(B[, l], n, n, byrow = TRUE) +
             matrix(B[, l], n, n, byrow = FALSE)
    theta <- theta + Bl^alpha
  }
  diag(theta) <- 2^alpha

  return(theta)
}

library(inline)
library(Rcpp)

# Rcpp function to get delta
GetDelta_code = "
  Rcpp::NumericVector this_ba(Ba);
  double a = Rcpp::as<double>(alpha);

  int n = this_ba.size();
  Rcpp::NumericMatrix delta(n, n);

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < i; j++){
      delta(i, j) = pow(this_ba(i) + this_ba(j), a - 1);
      delta(j, i) = delta(i, j);
    }
    delta(i, i) = pow(2 * this_ba(i), a - 1);
  }

  return delta;
"

GetDelta_cpp <- cxxfunction(signature(Ba = "numeric", alpha = "numeric"),
                            body = GetDelta_code,
                            plugin = "Rcpp")

GetRd_code = "
  Rcpp::NumericVector ba_c(Ba);
  double alpha_c = Rcpp::as<double>(alpha);
  Rcpp::NumericMatrix r_c(R);

  int n = ba_c.size();
  Rcpp::NumericMatrix delta(n, n);
  double this_delta_ij;

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < i; j++){
      this_delta_ij = pow(ba_c(i) + ba_c(j), alpha_c - 1);
      delta(i, j) = r_c(i, j) * this_delta_ij;
      delta(j, i) = r_c(j, i) * this_delta_ij;
    }
    delta(i, i) = r_c(i, i) * pow(2.0 * ba_c(i), alpha_c - 1);
  }

  return delta;
"

GetRd_cpp <- cxxfunction(signature(Ba = "numeric",
                                   alpha = "numeric",
                                   R = "numeric"),
                         body = GetRd_code,
                         plugin = "Rcpp")

GetDeltaR <- function(Ba, alpha) {
  n <- length(Ba)
  delta  <- (matrix(Ba, n, n, byrow = TRUE) +
               matrix(Ba, n, n, byrow = FALSE))^(alpha - 1)
  return(delta)
}

GetRdR <- function(Ba, alpha, R) {
  delta <- GetDeltaR(Ba, alpha)
  Rd <- R * delta
  return(Rd)
}

# Rcpp function to multiply a matrix by a constant

B2EC_code = "
  Rcpp::NumericMatrix Bcpp(Ba);
  double a = Rcpp::as<double>(alpha);

  int n = Bcpp.nrow();
  int L = Bcpp.ncol();
  Rcpp::NumericMatrix theta(n, n);

  for (int l = 0; l < L; l++) {
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < i; j++) {
        theta(i, j) += pow(Bcpp(i, l) + Bcpp(j, l), a);
      }
    }
  }

  for (int i = 0; i < n; i++){
    for (int j = 0; j < i; j++){
      theta(j,i) += theta(i,j);
    }
  }

  for (int i = 0; i < n; i++){
    theta(i,i) = pow(2,a);
  }

  return theta;
"

B2EC_cpp <- cxxfunction(signature(Ba = "numeric", alpha = "numeric"),
                        body = B2EC_code,
                        plugin = "Rcpp")

get.factors.EC <- function(EC.smooth, alpha.hat, L = 5, s = NULL,
                           n_starts = 10, fudge = 0.001,
                           maxit.1 = 20, maxit.2 = 100000,
                           eps = 1e-8, verbose = TRUE) {
  tick   <- proc.time()[3]
  n <- ncol(EC.smooth)
  d <- rdist(s, s)

  # Basis function estimation.
  bestval <- Inf
  # We try a few random starts to get us in the correct direction, but then
  # move on to the real minimization.
  for (rep in 1:n_starts) {
    ab <- 10 * eigen(2 - EC.smooth)$vec[, 1:L] + rnorm(n * L) / 2
    opt        <- optim(ab, fn = sse, gr = sse_grad, n = n, L = L,
                        EC = EC.smooth, alpha = alpha.hat, fudge = fudge,
                        method = "BFGS",
                        control = list(trace = ifelse(verbose, 50, 0),
                                       maxit = maxit.1, reltol = eps))
    if (opt$val < bestval) {
      b          <- opt$par
      bestval    <- opt$val
    }
  }

  opt        <- optim(b, fn = sse, gr = sse_grad, n = n, L = L,
                      EC = EC.smooth, alpha = alpha.hat, fudge = fudge,
                      method = "BFGS",
                      control = list(trace = ifelse(verbose, 50, 0),
                                     maxit = maxit.2, reltol = eps))
  b <- opt$par

  b     <- matrix(b, n, L)
  B     <- sweep(exp(b), 1, rowSums(exp(b)), "/")
  pct   <- colSums(B) / sum(B)
  tock  <- proc.time()[3]

  output <- list(est = B, pct = pct, seconds = tock - tick)
  return(output)
}
