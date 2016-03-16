# Already set in other file:
# cv: which cross-validation testing set to use
# L: the number of basis functions to use

#### load in the data ####
load(file = "../../code/analysis/fire/georgia_preprocess/fire_data.RData")

#### spatial setup ####
# get Georgia coordinates from georgia_preprocess in code/analysis/fire
load(file = "../../code/analysis/fire/georgia_preprocess/georgia_map.RData")
d <- rdist(cents)
diag(d) <- 0
n <- nrow(cents)

# standardize the locations
s <- cents
s.scale <- min(diff(range(s[, 1])), diff(range(s[, 2])))
s[, 1] <- (s[, 1] - min(s[, 1])) / s.scale
s[, 2] <- (s[, 2] - min(s[, 2])) / s.scale

################################################################################
#### Load in cross-validation setup ############################################
################################################################################
load(file = "./cv-extcoef.RData")

################################################################################
#### Get weight functions for spatial process ##################################
################################################################################
if (process == "ebf") {
  cat("Start basis function estimation \n")
  # basis function estimates using only the training data
  out       <- get.factors.EC(ec.hat[[cv]], L = L, s = s)
  B.sp      <- out$est
  ec.smooth <- out$EC.smooth
  alpha     <- out$alpha
} else {
  cat("Start estimation of rho and alpha \n")
  # alpha and rho estimates using only the training data
  out       <- get.rho.alpha(EC = ec.hat[[cv]], s = s, knots = s)
  rho       <- out$rho
  ec.smooth <- out$EC.smooth
  alpha     <- out$alpha
  dw2       <- out$dw2
  B.sp      <- getW(rho = rho, dw2 = dw2)  # using w as basis functions in MCMC
}

if (margin == "ebf") {
  if (process == "ebf") {  # we can just copy the empirical basis functions
    cat("Using spatial basis functions for covariates \n")
    B.cov <- B.sp
  } else {  # we need to construct the empirical basis functions
    cat("Estimating basis functions for covariates \n")
    B.cov <- get.factors.EC(ec.hat[[cv]], L = L, s = s)$est
  }
} else if (margin == "2bs") {
  if (L == 4) {
    # we can't do 2 degrees of freedom for basis functions, so opting for 
    # second order linear trend
    B.cov <- cbind(s[, 1], s[, 2], s[, 1] * s[, 2], s[, 1]^2, s[, 2]^2)
  } else if (L > 4) {
    # we already scale the spatial locations, so not necessary to rescale
    B.cov <- bspline.2d(s = s, scale = FALSE, df.x = sqrt(L), df.y = sqrt(L))
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

## transpose Y because preprocessed forest fire data is Y[t, i]
Y.all <- Y <- t(Y)
Y.tst <- Y[cv.idx[[cv]]]  # save the testing data to validate
Y[cv.idx[[cv]]] <- NA  # remove the testing data

ns <- nrow(Y)
nt <- ncol(Y)
np <- 2 + L * 2  # for a single year (int, t, B1...BL, t * (B1...BL))

## create covariate matrix for training
X <- array(1, dim = c(ns, nt, np))
for (i in 1:ns) {
  for (t in 1:nt) {
    time <- (t - nt / 2) / nt
    X[i, t, 2:np] <- c(time, B.cov[i, ], B.cov[i, ] * time)
  }
}

################################################################################
#### Spatially smooth threshold ################################################
################################################################################
thresh <- thresh95 <- thresh99 <- rep(0, ns)
neighbors <- 5
d <- rdist(s)
diag(d) <- 0

# take the 5 closest neighbors when finding the threshold
for (i in 1:ns) {
  these <- order(d[i, ])[2:(neighbors + 1)]  # the closest is always site i
  thresh[i] <- quantile(Y[these, ], probs = 0.90, na.rm = TRUE)
  thresh95[i] <- quantile(Y[these, ], probs = 0.95, na.rm = TRUE)
  thresh99[i] <- quantile(Y[these, ], probs = 0.99, na.rm = TRUE)
}
thresh   <- matrix(thresh, ns, nt)
thresh95 <- matrix(thresh95, nrow(Y), ncol(Y))
thresh99 <- matrix(thresh99, nrow(Y), ncol(Y))
thresh.tst   <- thresh[cv.idx[[cv]]]
thresh95.tst <- thresh95[cv.idx[[cv]]]
thresh99.tst <- thresh99[cv.idx[[cv]]]

################################################################################
#### run the MCMC ##############################################################
################################################################################
# iters  <- 30000
# burn   <- 20000
# update <- 1000

iters <- 200; burn <- 150; update <- 10  # for testing

cat("Start mcmc fit \n")
set.seed(6262)  # mcmc
# fit the model using the training data
fit <- ReShMCMC(y = Y, X = X, thresh = thresh, B = B.sp, alpha = alpha,
                iters = iters, burn = burn, update = update, iterplot = FALSE)
cat("Finished fit and predict \n")

# # calculate the scores
# probs.for.qs <- c(0.95, 0.96, 0.97, 0.98, 0.99, 0.995)
# qs.results <- QuantScore(preds = fit$y.pred, probs = probs.for.qs,
#                          validate = Y.tst)
# bs.results95 <- BrierScore(preds = fit$y.pred, validate = Y.tst,
#                            thresh = thresh95.tst)
# bs.results99 <- BrierScore(preds = fit$y.pred, validate = Y.tst,
#                            thresh = thresh99.tst)
# results <- c(qs.results, bs.results95, bs.results99, fit$timing)
# results <- c(results, Sys.info()["nodename"])
# names(results) <- c(probs.for.qs, "bs-95", "bs-99", "timing", "system")
# 
# write.table(results, file = table.file)
# 
# upload.pre <- "samorris@hpc.stat.ncsu.edu:~/repos-git/extreme-decomp/markdown/"
# upload.pre <- paste(upload.pre, "fire-analysis/cv-tables/", sep = "")
# if (do.upload) {
#   upload.cmd <- paste("scp ", table.file, " ", upload.pre, sep = "")
#   system(upload.cmd)
# }
# save(B.sp, thresh, alpha, fit, cv.idx, results, file = results.file)