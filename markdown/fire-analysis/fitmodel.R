# Already set in other file:
# cv: which cross-validation testing set to use
# L: the umber of basis functions to use

# get the Georgia map and coordinates
# from georgia_preprocess in code/analysis/fire
load(file = "../../code/analysis/fire/georgia_preprocess/georgia_map.RData")
d <- rdist(cents)
diag(d) <- 0
n <- nrow(cents)

# chi[i, j] is the chi statistic between sites i and j
load(file = "../../code/analysis/fire/georgia_preprocess/chi.RData")
chi.hat <- ifelse(chi <= 0, 0, chi)
ec.hat  <- 2 - chi.hat

# set up the 5 fold cross validation
n <- length(county)
set.seed(28)  #cv
nfolds <- 5
cv.idx <- get.cv.test(n = n, nfolds = nfolds)

# separate the test vs training
ec.hat.trn <- ec.hat[-cv.idx[[cv]], -cv.idx[[cv]], drop = FALSE]
ec.hat.tst <- ec.hat[cv.idx[[cv]], cv.idx[[cv]], drop = FALSE]
cents.trn  <- cents[-cv.idx[[cv]], , drop = FALSE]
cents.tst  <- cents[cv.idx[[cv]], , drop = FALSE]

################################################################################
## Estimage the basis functions
################################################################################

cat("Start basis function estimation \n")

# basis function estimates using only the training data
out.trn       <- get.factors.EC(ec.hat.trn, L = L, s = cents.trn)
B.est.trn     <- out.trn$est
ec.smooth.trn <- out.trn$EC.smooth
alpha         <- out.trn$alpha  # used both for training and testing

# basis function estimates using the testing data
out.tst       <- get.factors.EC(ec.hat.tst, L = L, s = cents.tst)
B.est.tst     <- out.tst$est
ec.smooth.tst <- out.tst$EC.smooth

################################################################################
#### Run the MCMC:
#### Use the basis functions with the MCMC
#### The response is the total acreage burned in a year 
####   Y[i, t] = acres burned in county i and year t 
####   X[i, t, p] = pth covariate for site i in year t
####     Using (1, time, B, B * time) where time = (t - nt / 2) / nt
################################################################################

## transpose Y because preprocessed forest fire data is Y[t, i]
Y <- t(Y)
Y.trn <- Y[-cv.idx[[cv]], ]
Y.tst <- Y[cv.idx[[cv]], ]
ns     <- nrow(Y)
ns.trn <- nrow(Y.trn)
ns.tst <- nrow(Y.tst)

nt <- ncol(Y.trn)
np <- 2 + L * 2  # for a single year (int, t, B1...BL, t * (B1...BL))

## create covariate matrix for training
X.trn <- array(1, dim = c(ns.trn, nt, np))
for (i in 1:ns.trn) {
  for (t in 1:nt) {  # cross-validation with sites not 
    time <- (t - nt / 2) / nt
    X.trn[i, t, 2:np] <- c(time, B.est.trn[i, ], B.est.trn[i, ] * time) 
  }
}

## create covariate matrix for testing
X.tst <- array(1, dim = c(ns.tst, nt, np))
for (i in 1:ns.tst) {
  for (t in 1:nt) {  # cross-validation with sites not 
    time <- (t - nt / 2) / nt
    X.tst[i, t, 2:np] <- c(time, B.est.tst[i, ], B.est.tst[i, ] * time) 
  }
}

################################################################################
## need spatially smoothed threshold - only using training data
################################################################################
thresh <- rep(0, ns.trn)
neighbors <- 5
d.trn <- rdist(cents.trn)
diag(d.trn) <- 0

# take the 5 closest neighbors when finding the threshold
for (i in 1:ns.trn) {
  these <- order(d.trn[i, ])[2:(neighbors + 1)]  # the closest is always site i
  thresh[i] <- quantile(Y.trn[these, ], probs = 0.95)
}

################################################################################
## run the MCMC
################################################################################
iters  <- 30000
burn   <- 20000
update <- 500

cat("Start mcmc fit \n")
set.seed(6262)  # mcmc
# fit the model using the training data
fit <- ReShMCMC(y = Y.trn, X = X.trn, thresh = thresh, B = B.est.trn, 
                alpha = alpha, 
                iters = iters, burn = burn, update = update, iterplot = FALSE)

cat("Start mcmc predict \n")

# predict at other locations
y.pred <- pred.ReShMCMC(mcmcoutput = fit, X.pred = X.tst, 
                        B = B.est.tst, alpha = alpha, 
                        start = 1, end = (iters - burn), update = update)

cat("Finished fit and predict \n")

# calculate the scores
probs.for.qs <- c(0.95, 0.96, 0.97, 0.98, 0.99, 0.995)
quant.scores <- QuantScore(pred = y.pred, probs = probs.for.qs, validate = Y.tst)
write.table(quant.scores, file = table.file)

upload.pre <- "samorris@hpc.stat.ncsu.edu:~/repos-git/extreme-decomp/markdown/"
upload.pre <- paste(upload.pre, "fire-analysis/cv-tables/", sep = "")
if (do.upload) {
  upload.cmd <- paste("scp ", table.file, " ", upload.pre, sep = "")
  system(upload.cmd)
}
save(B.est.trn, B.est.tst, alpha, fit, y.pred, cv.idx, quant.scores, 
     file = results.file)