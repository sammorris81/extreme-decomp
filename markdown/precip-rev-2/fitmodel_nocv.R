# Already set in other file:
# process:= what kind of spatial process (ebf, gsk)
# margin := how to construct marginal basis functions
# cv     := which cross-validation testing set to use
# L      := the number of basis functions to use
library(compiler)
enableJIT(3)

#### load in the data ####
load(file = "precip_preprocess.RData")
load(file = "all-extcoef.RData")

# basis functions are precomputed, so if we change cv settings, we'll
# need to rerun all of cv-setup.
basis.file   <- paste("./basis_functions/ebf-", L, "-all.RData", sep = "")
gsk.file     <- paste("./basis_functions/gsk-", L, "-all.RData", sep = "")
results.file <- paste("./cv-results/", process, "-", time, "-", L,
                      "-all.RData", sep = "")
table.file   <- paste("./cv-tables/", process, "-", time, "-", L,
                      "-all.txt", sep = "")

#### spatial setup ####
d <- rdist(s)
diag(d) <- 0
n <- nrow(s)

# standardize the locations
s.scale        <- s
s.scale.factor <- min(diff(range(s[, 1])), diff(range(s[, 2])))
s.min          <- apply(s, 2, min)
s.scale[, 1]   <- (s[, 1] - s.min[1]) / s.scale.factor
s.scale[, 2]   <- (s[, 2] - s.min[2]) / s.scale.factor

# get candidate knot grid for Gaussian kernel functions
cents.grid <- s.scale

################################################################################
#### Load in cross-validation setup ############################################
################################################################################
# load(file = "./cv-extcoef.RData")
load(file = basis.file)
load(file = gsk.file)

################################################################################
#### Get weight functions for spatial process ##################################
################################################################################
if (process == "ebf") {
  B.sp      <- B.ebf
  ec.smooth <- ec.smooth
  alpha     <- alpha.hats
} else {
  # get the knot locations
  alpha <- alpha.hats
  B.sp  <- B.gsk
}

################################################################################
#### Covariate basis functions #################################################
################################################################################
if (margin == "ebf") {
  if (process == "ebf") {  # we can just copy the empirical basis functions
    cat("B.cov = B.sp \n")
    B.cov <- B.sp
  } else {  # we need to construct the empirical basis functions
    cat("Estimating basis functions for covariates \n")
    B.cov <- B.ebf
  }
} else if (margin == "gsk") {
  if (process == "ebf") {
    B.cov <- B.gsk
  } else{
    cat("B.cov = B.sp \n")
    B.cov <- B.sp
  }
}

################################################################################
#### Run the MCMC ##############################################################
#### Use the basis functions with the MCMC
#### The response is the total acreage burned in a year
####   Y[i, t] = acres burned in county i and year t
####   X[i, t, p] = pth covariate for site i in year t
####     Using (1, time, B.cov, B.cov * time) where time = (t - nt / 2) / nt
################################################################################

ns <- nrow(Y)
nt <- ncol(Y) / 2
Y.all <- Y

## Y contains both current and future data, so subset on the relevant years
if (time == "current") {
  Y <- Y[, 1:nt]
} else {
  Y <- Y[, (nt + 1):(2 * nt)]
}

## standardize elevations
elev.std <- (elev - mean(elev)) / sd(elev)

X <- array(1, dim = c(ns, nt, 3))
for (i in 1:ns) {
  for (t in 1:nt) {
    time <- (t - nt / 2) / nt
    X[i, t, 2:3] <- c(time, elev.std[i])
  }
}

################################################################################
#### Spatially smooth threshold ################################################
################################################################################
thresh90 <- thresh95 <- thresh99 <- rep(0, ns)
neighbors <- 5
d <- rdist(s)
diag(d) <- 0

# take the 5 closest neighbors when finding the threshold
for (i in 1:ns) {
  these <- order(d[i, ])[2:(neighbors + 1)]  # the closest is always site i
  thresh90[i] <- quantile(Y[these, ], probs = 0.90, na.rm = TRUE)
  thresh95[i] <- quantile(Y[these, ], probs = 0.95, na.rm = TRUE)
  thresh99[i] <- quantile(Y[these, ], probs = 0.99, na.rm = TRUE)
}
thresh90 <- matrix(thresh90, nrow(Y), ncol(Y))
thresh95 <- matrix(thresh95, nrow(Y), ncol(Y))
thresh99 <- matrix(thresh99, nrow(Y), ncol(Y))

################################################################################
#### run the MCMC ##############################################################
################################################################################
iters  <- 25000
burn   <- 15000
update <- 500

# iters <- 100; burn <- 90; update <- 10  # for testing
A.init <- matrix(exp(2), L, nt)  # consistent with estimates of alpha
# A.init <- matrix(1, L, nt)
theta.init <- (B.sp^(1 / alpha) %*% A.init)^alpha
xi.init <- 0.1

# find the beta estimates using ml for GEV
# going to use ML but independent to get a single mu and sigma for each site
# based on xi = 0, but with A.init and alpha

beta.int.init <- matrix(0, ns, 2)
beta.time.init <- matrix(0, ns, 2)
for (i in 1:ns) {
  fit <- optim(par = c(0, 0), fn = loglike.init,
               y = Y[i, ], thresh = rep(-Inf, nt), xi = xi.init,
               theta = theta.init[i, ], alpha = alpha,
               control = list(maxit = 5000))$par
  beta.int.init[i, ] <- fit[1:2]
  if (i %% 50 == 0) {
    print(paste("site ", i, " complete", sep = ""))
  }
}

cat("Start mcmc fit \n")
set.seed(6262)  # mcmc

# fit the model using the training data
fit <- ReShMCMC(y = Y, s = s.scale, thresh = -Inf, B = B.sp, alpha = alpha,
                beta.int = beta.int.init, canbeta.int.sd = 0.5,
                beta.time = beta.time.init, canbeta.time.sd = 0.5,
                xi = xi.init, bw.init = 0.2, A = A.init,
                iters = iters, burn = burn, update = update,
                iterplot = FALSE)
cat("Finished fit and predict \n")

results <- c(fit$timing, Sys.info()["nodename"])
names(results) <- c("timing", "system")

write.table(results, file = table.file)

if (do.upload) {
  upload.pre <- paste("samorris@hpc.stat.ncsu.edu:~/repos-git/extreme-decomp/",
                      "markdown/precipitation/cv-tables/", sep = "")
  upload.cmd <- paste("scp ", table.file, " ", upload.pre, sep = "")
  system(upload.cmd)
}

save(B.sp, knots, thresh90, thresh95, thresh99,
     alpha, fit, results, file = results.file)