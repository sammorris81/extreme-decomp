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
  Y[this.cv] <- NA  # remove the testing data
} else {
  this.cv <- cv.idx[[cv]][, (nt + 1):(2 * nt)]
  Y <- Y[, (nt + 1):(2 * nt)]
  Y.tst <- Y[this.cv]
  Y[this.cv] <- NA
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
#### get the MLE ###############################################################
################################################################################

# we are using an intercept, std. elev, and time to get starting values
# this seems to be a good place to start, and then we can initialize all the
# basis function coefficients to be 0.

# only fitting intercept, std.elev, and time
Y.spatex <- t(Y)
X.spatex <- matrix(X[, 1, 3], ns, 1)
X.timeex <- t(t(X[1, , 2]))

colnames(X.spatex) <- "elev"
colnames(X.timeex) <- "year"

loc.form   <- ~ elev
scale.form <- ~ elev
shape.form <- shape ~ 1

temp.form.loc   <- temp.loc ~ year
temp.form.scale <- temp.scale ~ year

options(warn = 0)
library(SpatialExtremes)
fit.mle <- fitspatgev(data = Y.spatex, covariables = X.spatex,
                      loc.form = loc.form, scale.form = scale.form,
                      shape.form = shape.form,
                      temp.cov = X.timeex,
                      temp.form.loc = temp.form.loc,
                      temp.form.scale = temp.form.scale)
options(warn = 2)

beta1.init <- c(fit.mle$fitted.values[1], fit.mle$fitted.values[6],
                fit.mle$fitted.values[2])
beta1.init <- c(beta1.init, rep(0, 2 * L))
beta2.init <- c(fit.mle$fitted.values[3], fit.mle$fitted.values[7],
                fit.mle$fitted.values[4])
beta2.init <- c(beta2.init, rep(0, 2 * L))
xi.init <- fit.mle$fitted.values[5]

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
update <- 1000

# iters <- 100; burn <- 10; update <- 10  # for testing
A.init <- exp(6)  # consistent with estimates of alpha

cat("Start mcmc fit \n")
set.seed(6262)  # mcmc
#134.94
# fit the model using the training data
# Rprof(filename = "Rprof.out", line.profiling = TRUE)
fit <- ReShMCMC(y = Y, X = X, s = s.scale, knots = knots,
                thresh = -Inf, B = B.sp, alpha = alpha,
                can.mu.sd = 0.05, can.ls.sd = 0.01,
                tau1.a = 0.1, tau1.b = 0.1,
                tau2.a = 0.1, tau2.b = 0.1,
                beta1 = beta1.init, beta1.sd = 10,
                beta1.tau.a = 0.1, beta1.tau.b = 0.1,
                beta2 = beta2.init, beta2.sd = 1,
                beta2.tau.a = 0.1, beta2.tau.b = 0.1,
                xi = 0, xi.min = -0.5, xi.max = 0.5, xi.mn = 0, xi.sd = 0.5,
                bw.gp.init = 0.3,
                A = A.init, bw.basis.init = 0.3,
                time.interact = TRUE,
                keep.sites = 5, keep.days = 5, keep.knots = 10,
                iters = iters, burn = burn, update = update, iterplot = FALSE)
               # iters = iters, burn = burn, update = update, keep.burn = TRUE,
               # iterplot = TRUE)
cat("Finished fit and predict \n")
# Rprof(filename = NULL)
# summaryRprof(filename = "Rprof.out", lines = "show")

# par(mfrow = c(7, 5))
# for (i in 1:np) {
#   plot(fit$beta1[, i], type = "l",
#        main = bquote(paste(mu, ": ", beta[.(i)])))
# }
# for (i in 1:np) {
#   plot(fit$beta2[, i], type = "l",
#        main = bquote(paste("log(", sigma, "): ", beta[.(i)])))
# }
# plot(fit.rw.noblock$xi, type = "l", main = bquote(xi))
#
# mu.post <- sig.post <- array(0, dim = c(10000, ns, nt))
# dw2 <- rdist(s.scale, knots)^2
# dw2[dw2 < 1e-4] <- 0
# for (i in 1:10000) {
#   # update X matrix
#   B.i <- makeW(dw2 = dw2, rho = fit$bw[i])
#   X.mu <- X.sig <- add.basis.X(X = X, B.i)
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

if (do.upload) {
  upload.pre <- paste("samorris@hpc.stat.ncsu.edu:~/repos-git/extreme-decomp/",
                      "markdown/precipitation/cv-tables/", sep = "")
  upload.cmd <- paste("scp ", table.file, " ", upload.pre, sep = "")
  system(upload.cmd)
}

save(B.sp, knots, thresh90, thresh95, thresh99, Y.tst,
     alpha, fit, cv.idx, results, file = results.file)

# np <- 2 + L * 2  # for a single year (int, t, B1...BL, t * (B1...BL))
# np <- 3 + L  # for a single year (t, elev, log(elev), B1...BL) - No intercept
# np <- 2 + L  # for a single year (t, elev, B1, ..., BL)
# np <- 8  # for a single year (int, t, elev, long, lat, long * lat, long^2, lat^2)
# np <- 6

# ## standardize spatial basis functions
# for (i in 1:L) {
#   B.cov[, i] <- (B.cov[, i] - mean(B.cov[, i])) / sd(B.cov[, i])
# }

# logelev.std <- (log(elev) - mean(log(elev))) / sd(log(elev))

## create covariate matrix for training
# X <- array(1, dim = c(ns, nt, np))
# for (i in 1:ns) {
#   for (t in 1:nt) {
#     time <- (t - nt / 2) / nt
#     X[i, t, 2:np] <- c(time, B.cov[i, ], B.cov[i, ] * time)
#   }
# }

# # want to try using long, lat centered and scaled
# s.shift <- s.scale * 2
# s.shift[, 1] <- s.shift[, 1] - mean(s.shift[, 1])
# s.shift[, 2] <- s.shift[, 2] - mean(s.shift[, 2])

# X <- array(1, dim = c(ns, nt, np))
# for (i in 1:ns) {
#   for (t in 1:nt) {
#     time <- (t - nt / 2) / nt
#     X[i, t, 2:np] <- c(time, elev[i], s.shift[i, 1], s.shift[i, 2],
#                        s.shift[i, 1] * s.shift[i, 2],
#                        s.shift[i, 1]^2, s.shift[i, 2]^2)
#   }
# }

################################################################################
#### get the MLE ###############################################################
################################################################################

# Y.spatex <- as.vector(Y)
# X.spatex <- X[, 1, ]
# for (t in 2:nt) {
#   X.spatex <- rbind(X.spatex, X[, t, ])
# }
# X.spatex <- X.spatex[!is.na(Y.spatex), ]
# Y.spatex <- Y.spatex[!is.na(Y.spatex)]
#
# names <- c("time", "elev")
# colnames(X.spatex) <- names
#
# X.spatex <- data.frame(X.spatex)
#
# options(warn = 0)
# fit.mle <- fevd(Y.spatex, data = X.spatex,
#                 # location.fun = loc.fun,
#                 # scale.fun = scale.fun,
#                 use.phi = TRUE)
# beta1.init <- fit.mle$results$par[1:np]
# # beta1.init[1] <- 5
# # beta1.init[2] <- -3
# beta2.init <- fit.mle$results$par[(np + 1):(2 * np)]
# # beta2.init[1] <- 0
# # beta2.init[2] <- -0.1
# # xi.init    <- tail(fit.mle$results$par, 1)
# # xi.init       <- -1
# options(warn = 2)
# #
# # mu.st <- logsig.st <- matrix(0, ns, nt)
# # for (i in 1:np) {
# #   mu.st <- mu.st + X[, , i] * beta1.init[i]
# #   logsig.st <- logsig.st + X[, , i] * beta2.init[i]
# #  }
# #
# mus <- logsigs <- rep(0, ns)
# xi <- 0
# for (i in 1:ns) {
#   fit.lmoment <- fevd(Y[i, !is.na(Y[i, ])], method = "Lmoments")
#   mus[i] <- fit.lmoment$results[1]
#   logsigs[i] <- log(fit.lmoment$results[2])
#   xi <- xi + fit.lmoment$results[3] / ns
# }
#
# fit.mu.lm <- lm(mus ~ X[, 1, 2:np] - 1)
# fit.logsig.lm <- lm(logsigs ~ X[, 1, 2:np] - 1)
#
# beta1.init <- c(0, coef(fit.mu.lm))
# beta2.init <- c(0, coef(fit.logsig.lm))
# xi.init    <- -0.001

# names <- c("time", "elev", paste("B", 1:L, sep = ""))
# names(beta1.init) <- names
# names(beta2.init) <- names

# Y.spatex <- as.vector(Y)
# X.spatex <- X[, 1, ]
# for (t in 2:nt) {
#   X.spatex <- rbind(X.spatex, X[, t, ])
# }
# X.spatex <- X.spatex[!is.na(Y.spatex), ]
# Y.spatex <- Y.spatex[!is.na(Y.spatex)]
#
# beta.init <- rep(0, np * 2)
# library(evd)
# fit1 <- optim(par = beta.init, fn = ll.ind, X = X.spatex, y = Y.spatex, xi = 0,
#       control = list(maxit = 50000))
# beta.init <- c(fit1$par, 0)
# fit1 <- optim(par = beta.init, fn = ll.ind.xi, X = X.spatex, y = Y.spatex,
#               control = list(maxit = 50000))
# beta.init <- c(fit1$par)
# beta1.init <- fit1$par[1:np]
# beta2.init <- fit1$par[(np + 1):(2 * np)]
# xi.init <- tail(fit1$par, 1)

# options(warn = 0)
# fit.lmoment <- fevd(Y[!is.na(Y)], method = "Lmoments")
# options(warn = 2)
#
# beta1.init <- beta2.init <- rep(0, np)
# beta1.init[1] <- fit.lmoment$results[1]
# beta2.init[1] <- fit.lmoment$results[2]
# xi.init <- fit.lmoment$results[3]


