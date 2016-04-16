# Already set in other file:
# process:= what kind of spatial process (ebf, gsk)
# margin := how to construct marginal basis functions
# cv     := which cross-validation testing set to use
# L      := the number of basis functions to use

#### load in the data ####
load(file = "precip_preprocess.RData")

# basis functions are precomputed, so if we change cv settings, we'll
# need to rerun all of cv-setup.
basis.file   <- paste("./ebf-", L, ".RData", sep = "")
gsk.file     <- paste("./gsk-", L, ".RData", sep = "")
results.file <- paste("./cv-results/", process, "-", margin, "-", L, "-", cv,
                      ".RData", sep = "")
table.file   <- paste("./cv-tables/", process, "-", margin, "-", L, "-", cv,
                      ".txt", sep = "")

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
  knots <- cover.design(cents.grid, nd = L)$design
  alpha     <- alphas[cv]
  B.sp      <- B.gsk[[cv]]
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
    # get the knot locations
    knots <- cover.design(cents.grid, nd = L)$design
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
  Y <- Y[, 1:nt]
  Y.tst <- Y[cv.idx[[cv]][, 1:nt]]  # save the testing data to validate
  Y[cv.idx[[cv]][, 1:nt]] <- NA  # remove the testing data
} else {
  Y <- Y[, (nt + 1):(2 * nt)]
  Y.tst <- Y[cv.idx[[cv]][, (nt + 1):(2 * nt)]]
  Y[cv.idx[[cv]][, (nt + 1):(2 * nt)]] <- NA
}

# np <- 2 + L * 2  # for a single year (int, t, B1...BL, t * (B1...BL))
np <- 3 + L  # for a single year (t, elev, log(elev), B1...BL) - No intercept

## standardize spatial basis functions
for (i in 1:L) {
  B.cov[, i] <- (B.cov[, i] - mean(B.cov[, i])) / sd(B.cov[, i])
}

## create covariate matrix for training
# X <- array(1, dim = c(ns, nt, np))
# for (i in 1:ns) {
#   for (t in 1:nt) {
#     time <- (t - nt / 2) / nt
#     X[i, t, 2:np] <- c(time, B.cov[i, ], B.cov[i, ] * time)
#   }
# }

X <- array(0, dim = c(ns, nt, np))
for (i in 1:ns) {
  for (t in 1:nt) {
    time <- (t - nt / 2) / nt
    X[i, t, ] <- c(time, elev[i], log(elev[i]), B.cov[i, ])
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
thresh90.tst <- thresh90[cv.idx[[cv]]]
thresh95.tst <- thresh95[cv.idx[[cv]]]
thresh99.tst <- thresh99[cv.idx[[cv]]]

################################################################################
#### run the MCMC ##############################################################
################################################################################
iters  <- 30000
burn   <- 20000
update <- 1000

iters <- 100; burn <- 90; update <- 10  # for testing
A.init <- 100  # consistent with estimates of alpha
beta1.init <- rep(0, np)
beta2.init <- rep(0, np)
# beta1.init[1] <- 100
# beta2.init[1] <- 3.6
# beta1.init[1] <- 65
# beta2.init[1] <- 4
# beta1.init[1] <- 120
# beta2.init[1] <- 2.5


cat("Start mcmc fit \n")
set.seed(6262)  # mcmc

# fit the model using the training data
fit <- ReShMCMC(y = Y, X = X, thresh = -Inf, B = B.sp, alpha = alpha,
                xi = 0.001, can.mu.sd = 1, can.sig.sd = 0.1,
                beta1.attempts = 50, beta2.attempts = 50, A = A.init,
                beta1 = beta1.init, # beta2 = beta2.init,
                beta1.tau.a = 0.1, beta1.tau.b = 0.1,
                beta1.sd = 100, beta1.sd.fix = FALSE,
                beta2.tau.a = 0.1, beta2.tau.b = 0.1,
                beta2.sd = 10, beta2.sd.fix = FALSE,
                beta1.block = FALSE, beta2.block = FALSE,
                # iters = iters, burn = burn, update = update, iterplot = FALSE)
                iters = iters, burn = burn, update = update, iterplot = TRUE)
cat("Finished fit and predict \n")


cat("Start mcmc fit \n")
set.seed(6262)  # mcmc

# fit the model using the training data
fit <- ReShMCMC(y = Y, X = X, thresh = -Inf, B = B.sp, alpha = alpha,
                xi = 0.001, can.mu.sd = 1, can.sig.sd = 0.1,
                beta1.attempts = 50, beta2.attempts = 50, A = A.init,
                beta1 = beta1.init, beta2 = beta2.init,
                beta1.tau.a = 0.1, beta1.tau.b = 0.1,
                beta1.sd = 50, beta1.sd.fix = TRUE,
                beta2.tau.a = 0.1, beta2.tau.b = 0.1,
                beta2.sd = 5, beta2.sd.fix = TRUE,
                beta1.block = FALSE, beta2.block = FALSE,
                # iters = iters, burn = burn, update = update, iterplot = FALSE)
                iters = iters, burn = burn, update = update, iterplot = TRUE)
cat("Finished fit and predict \n")



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
upload.pre <- paste(upload.pre, "precipitation/cv-tables/", sep = "")
if (do.upload) {
  upload.cmd <- paste("scp ", table.file, " ", upload.pre, sep = "")
  system(upload.cmd)
}
save(B.sp, B.cov, out, thresh90, thresh95, thresh99,
     alpha, fit, cv.idx, results, file = results.file)



################################################################################
#### get the MLE ###############################################################
################################################################################
Y.spatex <- t(Y)
X.spatex <- X[, 1, 3:np]
X.timeex <- t(t(X[1, , 2]))
# X.spatex <- Xmat[!is.na(Yvec), ]
# Yvec <- Yvec[!is.na(Yvec)]

# Yvec <- matrix(Yvec, 1, length(Yvec))

names <- c("elev", "logelev", paste("B", seq(1:L), sep = ""))
colnames(X.spatex) <- names

# names <- paste("year", 1:nt, sep = "")
# colnames(X.timeex) <- names
colnames(X.timeex) <- "year"


loc.form   <- Yvec ~ elev + logelev + B1 + B2 + B3 + B4 + B5 + 0
scale.form <- Yvec ~ elev + logelev + B1 + B2 + B3 + B4 + B5 + 0
shape.form <- 1

# temp.form.loc <- Yvec ~ year1 + year2 + year3 + year4 + year5 + year6 +
#   year7 + year8 + year9 + year10 + year11 + year12 +
#   year13 + year14 + year15 + year16 + year17 + year18 +
#   year19 + year20 + year21 + year22 + year23 + year24 +
#   year25 + year26 + year27 + year28 + year29 + year30 +
#   year31 + year32
#
# temp.form.scale <- Yvec ~ year1 + year2 + year3 + year4 + year5 + year6 +
#   year7 + year8 + year9 + year10 + year11 + year12 +
#   year13 + year14 + year15 + year16 + year17 + year18 +
#   year19 + year20 + year21 + year22 + year23 + year24 +
#   year25 + year26 + year27 + year28 + year29 + year30 +
#   year31 + year32
temp.form.loc   <- Y ~ year
temp.form.scale <- Y ~ year
# temp.form.shape <- Y ~ 1

fit <- fitspatgev(data = Y.spatex, covariables = X.spatex,
                  loc.form = loc.form, scale.form = scale.form,
                  shape.form = shape.form,
                  temp.cov = X.timeex,
                  temp.form.loc = temp.form.loc,
                  temp.form.scale = temp.form.scale)
# temp.form.shape = temp.form.shape)

beta <- rep(0, 2 * np + 1)
optim(par = beta, fn = ll.ind, X = X, y = Y, hessian = TRUE)