# Already set in other file:
# cv: which cross-validation testing set to use
# L: the umber of basis functions to use

# load in the data
load(file = "../../code/analysis/fire/georgia_preprocess/fire_data.RData")

# get the Georgia map and coordinates
# from georgia_preprocess in code/analysis/fire
load(file = "../../code/analysis/fire/georgia_preprocess/georgia_map.RData")
d <- rdist(cents)
diag(d) <- 0
n <- nrow(cents)

# standardize the locations
s <- cents
s[, 1] <- (s[, 1] - min(s[, 1])) / diff(range(s[, 1]))
s[, 2] <- (s[, 2] - min(s[, 2])) / diff(range(s[, 2]))

# cv.idx and ec.hat were calculated ahead of time
load(file = "./cv-extcoef.RData")

################################################################################
## Estimate the basis functions
################################################################################

cat("Start basis function estimation \n")
# basis function estimates using only the training data
out       <- get.factors.EC(ec.hat[[cv]], L = L, s = s)
B.est     <- out$est
ec.smooth <- out$EC.smooth
alpha     <- out$alpha

################################################################################
#### Run the MCMC:
#### Use the basis functions with the MCMC
#### The response is the total acreage burned in a year 
####   Y[i, t] = acres burned in county i and year t 
####   X[i, t, p] = pth covariate for site i in year t
####     Using (1, time, B, B * time) where time = (t - nt / 2) / nt
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
    X[i, t, 2:np] <- c(time, B.est[i, ], B.est[i, ] * time) 
  }
}

################################################################################
## need spatially smoothed threshold - only using training data
################################################################################
thresh <- rep(0, ns)
neighbors <- 5
d <- rdist(s)
diag(d) <- 0

# take the 5 closest neighbors when finding the threshold
for (i in 1:ns) {
  these <- order(d[i, ])[2:(neighbors + 1)]  # the closest is always site i
  thresh[i] <- quantile(Y[these, ], probs = 0.90, na.rm = TRUE)
}
thresh <- matrix(thresh, ns, nt)
thresh.tst <- thresh[cv.idx[[cv]]]

################################################################################
## run the MCMC
################################################################################
# iters  <- 30000
# burn   <- 20000
# update <- 1000

iters <- 100; burn <- 90; update <- 50  # for testing

cat("Start mcmc fit \n")
set.seed(6262)  # mcmc
# fit the model using the training data
fit <- ReShMCMC(y = Y, X = X, thresh = thresh, B = B.est, alpha = alpha, 
                iters = iters, burn = burn, update = update, iterplot = FALSE)
cat("Finished fit and predict \n")

# # calculate the scores
# probs.for.qs <- c(0.95, 0.96, 0.97, 0.98, 0.99, 0.995)
# results <- QuantScore(preds = fit$y.pred, probs = probs.for.qs, 
#                       validate = Y.tst)
# 
# results <- c(results, fit$timing)
# results <- c(results, Sys.info()["nodename"])
# names(results) <- c(probs.for.qs, "timing", "system")

# write.table(results, file = table.file)
# 
# upload.pre <- "samorris@hpc.stat.ncsu.edu:~/repos-git/extreme-decomp/markdown/"
# upload.pre <- paste(upload.pre, "fire-analysis/cv-tables/", sep = "")
# if (do.upload) {
#   upload.cmd <- paste("scp ", table.file, " ", upload.pre, sep = "")
#   system(upload.cmd)
# }
# save(B.est, thresh, alpha, fit, cv.idx, results, file = results.file)