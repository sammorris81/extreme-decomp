# set up the 5 fold cross validation
n <- length(county)
set.seed(28)  #cv
nfolds <- 5
cv.idx <- get.cv.test(n = n, nfolds = nfolds)

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

# basis function estimates
out       <- get.factors.EC(ec.hat, L = L, s = cents)
B.est     <- out$est
alphahat  <- out$alpha
ec.smooth <- out$EC.smooth
ec.est    <- make.EC(B.est, alphahat)

# ################################################################################
# #### Run the MCMC:
# #### Use the basis functions with the MCMC
# #### The response is the total acreage burned in a year 
# ####   Y[i, t] = acres burned in county i and year t 
# ####   X[i, t, p] = pth covariate for site i in year t
# ####     Using (1, time, B, B * time) where time = (t - nt / 2) / nt
# ################################################################################
# 
# ## transpose Y because preprocessed forest fire data is Y[t, i]
# Y <- t(Y)
# nt <- ncol(Y)
# ns <- nrow(Y)
# np <- ncol(out$est) * 2 + 2
# 
# ## create covariate matrix
# X <- array(1, dim = c(ns, nt, np))
# for (i in 1:ns) {
#   for (t in 1:nt) {
#     time <- (t - nt / 2) / nt
#     X[i, t, 2:np] <- c(time, out$est[i, ], out$est[i, ] * time) 
#   }
# }
# 
# ## need spatially smoothed threshold
# thresh <- rep(0, ns)
# neighbors <- 5
# d <- rdist(cents)
# diag(d) <- 0
# 
# # take the 5 closest neighbors when finding the threshold
# for (i in 1:ns) {
#   these <- order(d[i, ])[2:(neighbors + 1)]  # the closest is always site i
#   thresh[i] <- quantile(Y[these, ], probs = 0.95)
# }
# 
# fit <- ReShMCMC(y = Y[train, ], X = X[train, , ], thresh = thresh[train], 
#                 B = out$est[train, , drop = FALSE], alpha = out$alpha, 
#                 iters = 300, burn = 100, update = 10, iterplot = TRUE)
# 
# y.pred <- pred.ReShMCMC(mcmcoutput = fit, X.pred = X[test, , ], 
#                         B = out$est, alpha = out$alpha, 
#                         start = 1, end = 200, update = 10)
