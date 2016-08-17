# Already set in other file:
# process:= what kind of spatial process (ebf, gsk)
# margin := how to construct marginal basis functions
# cv     := which cross-validation testing set to use
# L      := the number of basis functions to use
library(compiler)
enableJIT(3)

#### load in the data ####
load(file = "precip_preprocess.RData")

# basis functions are precomputed, so if we change cv settings, we'll
# need to rerun all of cv-setup.
basis.file   <- paste("./ebf-", L, ".RData", sep = "")
gsk.file     <- paste("./gsk-", L, ".RData", sep = "")
results.file <- paste("./cv-results/", process, "-", time, "-", L,
                      "-", cv, ".RData", sep = "")
table.file   <- paste("./cv-tables/", process, "-", time, "-", L,
                      "-", cv, ".txt", sep = "")

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
load(file = "./cv-extcoef.RData")
load(file = basis.file)
load(file = gsk.file)

################################################################################
#### Get weight functions for spatial process ##################################
################################################################################
if (process == "ebf") {
  B.sp      <- B.ebf[[cv]]
  ec.smooth <- ec.smooth[[cv]]
  alpha     <- alphas[cv]
} else {
  # get the knot locations
  alpha <- alphas[cv]
  B.sp  <- B.gsk[[cv]]
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
    B.cov <- B.ebf[[cv]]
  }
} else if (margin == "gsk") {
  if (process == "ebf") {
    B.cov <- B.gsk[[cv]]
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
  this.cv <- cv.idx[[cv]][, 1:nt]
  Y <- Y[, 1:nt]
  Y.tst <- Y[this.cv]  # save the testing data to validate
  # Y[this.cv] <- NA  # remove the testing data
} else {
  this.cv <- cv.idx[[cv]][, (nt + 1):(2 * nt)]
  Y <- Y[, (nt + 1):(2 * nt)]
  Y.tst <- Y[this.cv]
  # Y[this.cv] <- NA
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
thresh90.tst <- thresh90[this.cv]
thresh95.tst <- thresh95[this.cv]
thresh99.tst <- thresh99[this.cv]

################################################################################
#### run the MCMC ##############################################################
################################################################################
iters  <- 25000
burn   <- 15000
update <- 500

# iters <- 100; burn <- 10; update <- 50  # for testing
A.init <- matrix(exp(2), L, nt)  # consistent with estimates of alpha
# A.init <- matrix(1, L, nt)
theta.init <- (B.sp^(1 / alpha) %*% A.init)^alpha
xi.init <- 0.1

# find the beta estimates using ml for GEV
# going to use ML but independent to get a single mu and sigma for each site
# based on xi = 0, but with A.init and alpha

# xi.init <- rep(0, ns)
beta.int.init <- matrix(0, ns, 2)
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

beta.time.init <- matrix(0, ns, 2)
fit <- ReShMCMC(y = Y, test.set = this.cv, s = s.scale, thresh = -Inf, B = B.sp,
                alpha = alpha, beta.int = beta.int.init, canbeta.int.sd = 0.5,
                beta.time = beta.time.init, canbeta.time.sd = 0.5,
                xi = xi.init, bw.init = 0.2, A = A.init,
                iters = iters, burn = burn, update = update,
                iterplot = FALSE)

cat("Finished fit and predict \n")

#### Summarize performance ####
# GG and CRPS come out directly from mcmc. Also need quantile and Brier scores
# and MAD
probs.for.qs <- c(0.95, 0.96, 0.97, 0.98, 0.99, 0.995)
qs.results <- QuantScore(preds = fit$y.pred, probs = probs.for.qs,
                         validate = Y.tst)
bs.results95 <- BrierScore(preds = fit$y.pred, validate = Y.tst,
                           thresh = thresh95.tst)
bs.results99 <- BrierScore(preds = fit$y.pred, validate = Y.tst,
                           thresh = thresh99.tst)
y.pred  <- apply(fit$y.pred, 2, median)
MAD     <- mean(abs(y.pred - Y.tst))
results <- c(qs.results, bs.results95, bs.results99, fit$GG, fit$CRPS, MAD,
             fit$timing, Sys.info()["nodename"])
names(results) <- c(probs.for.qs, "bs-95", "bs-99", "GG", "CRPS", "MAD",
                    "timing", "system")

write.table(results, file = table.file)

if (do.upload) {
  upload.pre <- paste("samorris@hpc.stat.ncsu.edu:~/repos-git/extreme-decomp/",
                      "markdown/precipitation/cv-tables/", sep = "")
  upload.cmd <- paste("scp ", table.file, " ", upload.pre, sep = "")
  system(upload.cmd)
}

save(B.sp, knots, thresh90, thresh95, thresh99, Y.tst,
     alpha, fit, cv.idx, results, file = results.file)