# Already set in other file:
# process:= what kind of spatial process (ebf, gsk)
# margin := how to construct marginal basis functions
# cv     := which cross-validation testing set to use
# L      := the number of basis functions to use

#### load in the data ####
load(file = "../../code/analysis/fire/georgia_preprocess/fire_data.RData")

#### spatial setup ####
# get Georgia coordinates from georgia_preprocess in code/analysis/fire
load(file = "../../code/analysis/fire/georgia_preprocess/georgia_map.RData")
d <- rdist(cents)
diag(d) <- 0
n <- nrow(cents)

# get candidate knot grid for Gaussian kernel functions
grid.x <- seq(min(cents[, 1]), max(cents[, 1]), length = 100)
grid.y <- seq(min(cents[, 2]), max(cents[, 2]), length = 100)
cents.grid <- as.matrix(expand.grid(grid.x, grid.y))
inGA <- map.where("state", x = cents.grid[, 1], y = cents.grid[, 2])
cents.grid <- cents.grid[inGA == "georgia", ]
cents.grid <- cents.grid[rowSums(is.na(cents.grid)) == 0, ]

# standardize the locations
s <- cents
s.scale <- min(diff(range(s[, 1])), diff(range(s[, 2])))
s.min   <- apply(s, 2, min)
s[, 1] <- (s[, 1] - s.min[1]) / s.scale
s[, 2] <- (s[, 2] - s.min[2]) / s.scale
cents.grid[, 1] <- (cents.grid[, 1] - s.min[1]) / s.scale
cents.grid[, 2] <- (cents.grid[, 2] - s.min[2]) / s.scale

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
  # get the knot locations
  knots <- cover.design(cents.grid, nd = L)$design
  
  cat("Start estimation of rho and alpha \n")
  # alpha and rho estimates using only the training data
  out       <- get.rho.alpha(EC = ec.hat[[cv]], s = s, knots = knots)
  rho       <- out$rho
  ec.smooth <- out$EC.smooth
  alpha     <- out$alpha
  dw2       <- out$dw2
  B.sp      <- getW(rho = rho, dw2 = dw2)  # using w as basis functions in MCMC
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
    B.cov <- get.factors.EC(ec.hat[[cv]], L = L, s = s)$est
  }
} else if (margin == "gsk") {
  if (process == "ebf") {
    # get the knot locations
    knots <- cover.design(cents.grid, nd = L)$design
    
    cat("Start estimation of Gaussian kernels for covariates \n")
    # alpha and rho estimates using only the training data
    out   <- get.rho.alpha(EC = ec.hat[[cv]], s = s, knots = knots)
    B.cov <- getW(rho = out$rho, dw2 = out$dw2)
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

## transpose Y because preprocessed forest fire data is Y[t, i]
Y.all <- Y <- t(Y)
Y.tst <- Y[cv.idx[[cv]]]  # save the testing data to validate
Y[cv.idx[[cv]]] <- NA  # remove the testing data

ns <- nrow(Y)
nt <- ncol(Y)
np <- 2 + L * 2  # for a single year (int, t, B1...BL, t * (B1...BL))

## standardize spatial basis functions
for (i in 1:L) {
  B.cov[, i] <- (B.cov[, i] - mean(B.cov[, i])) / sd(B.cov[, i])
}

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
thresh90 <- matrix(thresh90, ns, nt)
thresh95 <- matrix(thresh95, nrow(Y), ncol(Y))
thresh99 <- matrix(thresh99, nrow(Y), ncol(Y))
thresh.tst   <- thresh90[cv.idx[[cv]]]
thresh95.tst <- thresh95[cv.idx[[cv]]]
thresh99.tst <- thresh99[cv.idx[[cv]]]

################################################################################
#### initial values ############################################################
################################################################################
# the default estimate usually puts the intercept at something much larger
# than 0 thereby making it difficult to get sigma^2_beta1 to be reasonably
# sized which greatly impacts the ability to get the MCMC to mix.
# beta1.init <- rep(0, np)

################################################################################
#### run the MCMC ##############################################################
################################################################################
iters  <- 30000
burn   <- 20000
update <- 100

# iters <-2000; burn <- 500; update <- 10  # for testing

## Prior 1
# beta1.tau.a = 1, beta1.tau.b = 1
# beta2.tau.a = 1, beta2.tau.b = 1

## Prior 2
# beta1.tau.a = 0.1, beta1.tau.b = 0.1
# beta2.tau.a = 1, beta2.tau.b = 1

cat("Start mcmc fit \n")
set.seed(6262)  # mcmc
# fit the model using the training data
fit1 <- ReShMCMC(y = Y, X = X, thresh = thresh95, B = B.sp, alpha = alpha,
                 # beta1 = beta1.init,
                 beta1.tau.a = 0.1, beta1.tau.b = 0.1,
                 beta2.tau.a = 10, beta2.tau.b = 10,
                 # beta2.tau.a = 0.1, beta2.tau.b = 0.1,
                 beta1.sd = 10, beta2.sd = 1, iters = iters, burn = burn, 
                 update = update, iterplot = FALSE)
cat("Finished fit and predict \n")

cat("Start mcmc fit \n")
set.seed(6262)  # mcmc
# fit the model using the training data
fit2 <- ReShMCMC(y = Y, X = X, thresh = thresh95, B = B.sp, alpha = alpha,
                 # beta1 = beta1.init,
                 beta1.tau.a = 1, beta1.tau.b = 1,
                 beta2.tau.a = 1, beta2.tau.b = 1,
                 # beta2.tau.a = 0.1, beta2.tau.b = 0.1,
                 beta1.sd = 10, beta2.sd = 1, iters = iters, burn = burn, 
                 update = update, iterplot = FALSE)
cat("Finished fit and predict \n")

cat("Start mcmc fit \n")
set.seed(6262)  # mcmc
# fit the model using the training data
fit3 <- ReShMCMC(y = Y, X = X, thresh = thresh95, B = B.sp, alpha = alpha,
                # beta1 = beta1.init,
                beta1.tau.a = 1, beta1.tau.b = 1,
                beta2.tau.a = 10, beta2.tau.b = 10,
                # beta2.tau.a = 0.1, beta2.tau.b = 0.1,
                beta1.sd = 10, beta2.sd = 1, iters = iters, burn = burn, 
                update = update, iterplot = FALSE)
cat("Finished fit and predict \n")

par(mfrow = c(3, 4))
these.iters <- seq(1, 10000, by = 2)
for (i in 1:12) {
  plot(fit1$beta1[these.iters, i], type = "l", 
       main = "Beta mu Fit 1", ylab = bquote(beta[.(i)]))
}

for (i in 1:12) {
  plot(fit2$beta1[, i], type = "l", 
       main = "Beta mu Fit 2", ylab = bquote(beta[.(i)]))
}

for (i in 1:12) {
  plot(fit3$beta1[, i], type = "l", 
       main = "Beta mu Fit 3", ylab = bquote(beta[.(i)]))
}

for (i in 1:12) {
  plot(fit1$beta2[, i], type = "l", 
       main = "Beta sig Fit 1", ylab = bquote(beta[.(i)]))
}
for (i in 1:12) {
  plot(fit2$beta2[, i], type = "l", 
       main = "Beta sig Fit 2", ylab = bquote(beta[.(i)]))
}
for (i in 1:12) {
  plot(fit3$beta2[, i], type = "l", 
       main = "Beta sig Fit 3", ylab = bquote(beta[.(i)]))
}

## get posterior distribution for each county
mus <- logsigs <- matrix(0, nrow = (iters - burn), ncol = ns)

betamu <- apply(fit$beta1, 2, mean)
betalogsig <- apply(fit$beta2, 2, mean)

# get an intercept and time coefficient for each county
betamu.ints <- betamu[1] + B.est %*% betamu[3:(L + 2)]
betamu.time <- betamu[2] + B.est %*% betamu[(L + 3):np]
betalogsig.ints <- betalogsig[1] + B.est %*% betalogsig[3:(L + 2)]
betalogsig.time <- betalogsig[2] + B.est %*% betalogsig[(L + 3):np]




# calculate the scores
probs.for.qs <- c(0.95, 0.96, 0.97, 0.98, 0.99, 0.995)
qs.results <- QuantScore(preds = fit$y.pred, probs = probs.for.qs,
                         validate = Y.tst)
bs.results95 <- BrierScore(preds = fit$y.pred, validate = Y.tst,
                           thresh = thresh95.tst)
bs.results99 <- BrierScore(preds = fit$y.pred, validate = Y.tst,
                           thresh = thresh99.tst)
results <- c(qs.results, bs.results95, bs.results99, fit$timing)
results <- c(results, Sys.info()["nodename"])
names(results) <- c(probs.for.qs, "bs-95", "bs-99", "timing", "system")

write.table(results, file = table.file)

upload.pre <- "samorris@hpc.stat.ncsu.edu:~/repos-git/extreme-decomp/markdown/"
upload.pre <- paste(upload.pre, "fire-analysis/cv-tables/", sep = "")
if (do.upload) {
  upload.cmd <- paste("scp ", table.file, " ", upload.pre, sep = "")
  system(upload.cmd)
}
save(B.sp, out, thresh, alpha, fit, cv.idx, results, file = results.file)