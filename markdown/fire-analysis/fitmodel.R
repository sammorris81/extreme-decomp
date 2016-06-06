# Already set in other file:
# process:= what kind of spatial process (ebf, gsk)
# margin := how to construct marginal basis functions
# cv     := which cross-validation testing set to use
# L      := the number of basis functions to use
options(warn = 2)
library(compiler)
enableJIT(3)

# basis functions are precomputed, so if we change cv settings, we'll
# need to rerun all of cv-setup.
basis.file   <- paste("./ebf-", L, ".RData", sep = "")
gsk.file     <- paste("./gsk-", L, ".RData", sep = "")

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
  alpha <- alphas[cv]
  B.sp  <- B.gsk[[cv]]
}

################################################################################
#### Covariate basis functions #################################################
################################################################################

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
np <- 2 # + L * 2  # for a single year (int, t)

## create covariate matrix for training
X <- array(1, dim = c(ns, nt, np))
for (i in 1:ns) {
  for (t in 1:nt) {
    time <- (t - nt / 2) / nt
    X[i, t, 2:np] <- time
  }
}

################################################################################
#### Spatially smooth threshold ################################################
################################################################################
thresh85 <- thresh90 <- thresh95 <- thresh99 <- rep(0, ns)
neighbors <- 5
d <- rdist(s)
diag(d) <- 0

# take the 5 closest neighbors when finding the threshold
for (i in 1:ns) {
  these <- order(d[i, ])[2:(neighbors + 1)]  # the closest is always site i
  thresh85[i] <- quantile(Y[these, ], probs = 0.85, na.rm = TRUE)
  thresh90[i] <- quantile(Y[these, ], probs = 0.90, na.rm = TRUE)
  thresh95[i] <- quantile(Y[these, ], probs = 0.95, na.rm = TRUE)
  thresh99[i] <- quantile(Y[these, ], probs = 0.99, na.rm = TRUE)
}
thresh85 <- matrix(thresh85, nrow(Y), ncol(Y))
thresh90 <- matrix(thresh90, nrow(Y), ncol(Y))
thresh95 <- matrix(thresh95, nrow(Y), ncol(Y))
thresh99 <- matrix(thresh99, nrow(Y), ncol(Y))
thresh85.tst <- thresh85[cv.idx[[cv]]]
thresh90.tst <- thresh90[cv.idx[[cv]]]
thresh95.tst <- thresh95[cv.idx[[cv]]]
thresh99.tst <- thresh99[cv.idx[[cv]]]

# Y.tmp <- Y
# Y.tmp <- ifelse(Y <= thresh85, thresh85, Y)
# colSums(Y.tmp - thresh85 != 0, na.rm = TRUE)
# Y.tmp <- ifelse(Y <= thresh90, thresh90, Y)
# colSums(Y.tmp - thresh90 != 0, na.rm = TRUE)
# Y.tmp <- ifelse(Y <= thresh95, thresh95, Y)
# colSums(Y.tmp - thresh95 != 0, na.rm = TRUE)

################################################################################
#### run the MCMC ##############################################################
################################################################################
iters  <- 30000
burn   <- 20000
update <- 1000

# iters <-200; burn <- 50;
update <- 500  # for testing

cat("Start mcmc fit \n")
set.seed(6262)  # mcmc
# fit the model using the training data
# s is scaled locations
fit <- ReShMCMC(y = Y, X = X, s = s, knots = knots,
                thresh = thresh90, B = B.sp, alpha = alpha,
                can.mu.sd = 1, can.ls.sd = 0.5,
                mu.attempts = 400, ls.attempts = 400,
                tau1.a = 1, tau1.b = 1,
                tau2.a = 0.1, tau2.b = 0.1,
                beta1.sd = 10, beta1.tau.a = 0.1, beta1.tau.b = 0.1,
                beta2.sd = 0.1, beta2.tau.a = 0.1, beta2.tau.b = 0.1,
                xi = 0, xi.min = -1, xi.max = 1, xi.mn = 0, xi.sd = 0.3,
                bw.gp.init = 0.05, can.bw.gp.sd = 0.1,
                bw.basis.init = 0.3, can.bw.basis.sd = 0.1,
                time.interact = TRUE,
                iters = iters, burn = burn, update = update, iterplot = TRUE)
                # iters = iters, burn = burn, update = update, iterplot = FALSE)
cat("Finished fit and predict \n")

# mu.post <- sig.post <- array(0, dim = c(10000, ns, nt))
# dw2 <- rdist(s, knots)^2
# dw2[dw2 < 1e-4] <- 0
# for (i in 1:10000) {
#   # update X matrix
#   B.i <- makeW(dw2 = dw2, rho = fit$bw[i])
#   X.mu <- X.sig <- add.basis.X(X = X, B.i, time.interact = TRUE)
#   for (t in 1:nt) {
#     mu.post[i, , t] <- X.mu[, t, ] %*% fit$beta1[i, ]
#     sig.post[i, , t] <- X.sig[, t, ] %*% fit$beta2[i, ]
#   }
#   if (i %% 500 == 0) {
#     print(paste(i, "finished"))
#   }
# }
#
# par(mfrow = c(5, 7))
# sites <- sample(ns, 5, replace = FALSE)
# days  <- sample(nt, 7, replace = FALSE)
#
# for (i in sites) {
#   for (t in days) {
#     plot(mu.post[, i, t], type = "l", main = bquote(paste(mu, "(", .(i), ", ", .(t), ")")))
#   }
# }
#
# for (i in sites) {
#   for (t in days) {
#     plot(sig.post[, i, t], type = "l", main = bquote(paste(sigma, "(", .(i), ", ", .(t), ")")))
#   }
# }

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
save(B.sp, knots, thresh90, thresh95, thresh99,
     alpha, fit, cv.idx, results, file = results.file)