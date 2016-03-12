source(file = "./package_load.R", chdir = T)

# Number of bases: we can do 1 basis function, but it's not a particularly 
# interesting case because the basis function is 1 at all locations due to the 
# fact that they need to add up to 1 across all basis functions at each site.
# opting for 2, 5, 10, 15, and then looking a max stable method with fixed 
# alpha.
method <- "basis" # using kern for the results from abba
L      <- 15  # will be using this to get basis functions for covariates
results.file <- paste("./cv-results/", method, "-", L, "-1st.RData", sep = "")

# load in the data
load(file = "../../code/analysis/fire/georgia_preprocess/fire_data.RData")
# transpose Y because preprocessed forest fire data is Y[t, i]
Y  <- t(Y)

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

# find the extremal coefficients
# ec.hat <- get.pw.ec(Y = Y, qlim = c(0.90, 1), verbose = TRUE, update = 50)$ec

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

#### extract the posterior distribution for the county time coefficients
load(file = results.file)

## get posterior means
betamu <- apply(fit$beta1, 2, mean)
betalogsig <- apply(fit$beta2, 2, mean)

# get an intercept and time coefficient for each county
betamu.ints <- betamu[1] + B.est %*% betamu[3:(L + 2)]
betamu.time <- betamu[2] + B.est %*% betamu[(L + 3):np]
betalogsig.ints <- betalogsig[1] + B.est %*% betalogsig[3:(L + 2)]
betalogsig.time <- betalogsig[2] + B.est %*% betalogsig[(L + 3):np]

map.ga.ggplot(Y = betamu.ints, 
              main = bquote(paste("Basis method: location intercept by county")), 
              midpoint = 0)
map.ga.ggplot(Y = betamu.time, 
              main = bquote(paste("Basis method: location time coefficient by county")), 
              midpoint = 0)

map.ga.ggplot(Y = betalogsig.ints, 
              main = bquote(paste("Basis method log(", sigma, "): intercept by county")),
              midpoint = 0)
map.ga.ggplot(Y = betalogsig.time, 
              main = bquote(paste("Basis method log(", sigma, "): time coefficient by county")),
              midpoint = 0)

# Number of bases: we can do 1 basis function, but it's not a particularly 
# interesting case because the basis function is 1 at all locations due to the 
# fact that they need to add up to 1 across all basis functions at each site.
# opting for 2, 5, 10, 15, and then looking a max stable method with fixed 
# alpha.
method <- "kern" # using kern for the results from abba
L      <- 15  # will be using this to get basis functions for covariates
results.file <- paste("./cv-results/", method, "-", L, "-1st.RData", sep = "")

# load in the data
load(file = "../../code/analysis/fire/georgia_preprocess/fire_data.RData")
# transpose Y because preprocessed forest fire data is Y[t, i]
Y  <- t(Y)

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

# find the extremal coefficients
# ec.hat <- get.pw.ec(Y = Y, qlim = c(0.90, 1), verbose = TRUE, update = 50)$ec

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

#### extract the posterior distribution for the county time coefficients
load(file = results.file)

## get posterior means
betamu <- apply(fit$beta1, 2, mean)
betalogsig <- apply(fit$beta2, 2, mean)

# get an intercept and time coefficient for each county
betamu.ints <- betamu[1] + B.est %*% betamu[3:(L + 2)]
betamu.time <- betamu[2] + B.est %*% betamu[(L + 3):np]
betalogsig.ints <- betalogsig[1] + B.est %*% betalogsig[3:(L + 2)]
betalogsig.time <- betalogsig[2] + B.est %*% betalogsig[(L + 3):np]

map.ga.ggplot(Y = betamu.ints, 
              main = bquote(paste("Kern method: location intercept by county")), 
              midpoint = 0)
map.ga.ggplot(Y = betamu.time, 
              main = bquote(paste("Kern method: location time coefficient by county")), 
              midpoint = 0)

map.ga.ggplot(Y = betalogsig.ints, 
              main = bquote(paste("Kern method log(", sigma, "): intercept by county")),
              midpoint = 0)
map.ga.ggplot(Y = betalogsig.time, 
              main = bquote(paste("Kern method log(", sigma, "): time coefficient by county")),
              midpoint = 0)


source(file = "./package_load.R", chdir = T)

# Number of bases: we can do 1 basis function, but it's not a particularly 
# interesting case because the basis function is 1 at all locations due to the 
# fact that they need to add up to 1 across all basis functions at each site.
# opting for 2, 5, 10, 15, and then looking a max stable method with fixed 
# alpha.
method <- "basis" # using kern for the results from abba
L      <- 15  # will be using this to get basis functions for covariates
results.file <- paste("./cv-results/", method, "-", L, "-2nd.RData", sep = "")

# load in the data
load(file = "../../code/analysis/fire/georgia_preprocess/fire_data.RData")
# transpose Y because preprocessed forest fire data is Y[t, i]
Y  <- t(Y)

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

# find the extremal coefficients
# ec.hat <- get.pw.ec(Y = Y, qlim = c(0.90, 1), verbose = TRUE, update = 50)$ec

ns <- nrow(Y)
nt <- ncol(Y)
np <- 3 + L * 3  # for a single year (int, t, t^2, B1...BL, t * (B1...BL), 
                 #                    t^2 * (B1...BL))

#### extract the posterior distribution for the county time coefficients
load(file = results.file)

## get posterior means
betamu <- apply(fit$beta1, 2, mean)
betalogsig <- apply(fit$beta2, 2, mean)

# get an intercept and time coefficient for each county
betamu.ints <- betamu[1] + B.est %*% betamu[4:(L + 3)]
betamu.time <- betamu[2] + B.est %*% betamu[(L + 4):(2 * L + 3)]
betamu.t2nd <- betamu[3] + B.est %*% betamu[(2 * L + 4):np]
betalogsig.ints <- betalogsig[1] + B.est %*% betalogsig[4:(L + 3)]
betalogsig.time <- betalogsig[2] + B.est %*% betalogsig[(L + 4):(2 * L + 3)]
betalogsig.t2nd <- betalogsig[3] + B.est %*% betalogsig[(2 * L + 4):np]

map.ga.ggplot(Y = betamu.ints, 
              main = bquote(paste("Basis method: location intercept by county")), 
              midpoint = 0)
map.ga.ggplot(Y = betamu.time, 
              main = bquote(paste("Basis method: location time coefficient by county")), 
              midpoint = 0)
map.ga.ggplot(Y = betamu.t2nd, 
              main = bquote(paste("Basis method: location time"^2, " coefficient by county")), 
              midpoint = 0)

map.ga.ggplot(Y = betalogsig.ints, 
              main = bquote(paste("Basis method log(", sigma, "): intercept by county")),
              midpoint = 0)
map.ga.ggplot(Y = betalogsig.time, 
              main = bquote(paste("Basis method log(", sigma, "): time coefficient by county")),
              midpoint = 0)
map.ga.ggplot(Y = betalogsig.t2nd, 
              main = bquote(paste("Basis method log(", sigma, "): time"^2, " coefficient by county")),
              midpoint = 0)

# Number of bases: we can do 1 basis function, but it's not a particularly 
# interesting case because the basis function is 1 at all locations due to the 
# fact that they need to add up to 1 across all basis functions at each site.
# opting for 2, 5, 10, 15, and then looking a max stable method with fixed 
# alpha.
method <- "kern" # using kern for the results from abba
L      <- 15  # will be using this to get basis functions for covariates
results.file <- paste("./cv-results/", method, "-", L, "-2nd.RData", sep = "")

# load in the data
load(file = "../../code/analysis/fire/georgia_preprocess/fire_data.RData")
# transpose Y because preprocessed forest fire data is Y[t, i]
Y  <- t(Y)

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

# find the extremal coefficients
# ec.hat <- get.pw.ec(Y = Y, qlim = c(0.90, 1), verbose = TRUE, update = 50)$ec

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

#### extract the posterior distribution for the county time coefficients
load(file = results.file)

## get posterior means
betamu <- apply(fit$beta1, 2, mean)
betalogsig <- apply(fit$beta2, 2, mean)

# get an intercept and time coefficient for each county
betamu.ints <- betamu[1] + B.est %*% betamu[4:(L + 3)]
betamu.time <- betamu[2] + B.est %*% betamu[(L + 4):(2 * L + 3)]
betamu.t2nd <- betamu[3] + B.est %*% betamu[(2 * L + 4):np]
betalogsig.ints <- betalogsig[1] + B.est %*% betalogsig[4:(L + 3)]
betalogsig.time <- betalogsig[2] + B.est %*% betalogsig[(L + 4):(2 * L + 3)]
betalogsig.t2nd <- betalogsig[3] + B.est %*% betalogsig[(2 * L + 4):np]

map.ga.ggplot(Y = betamu.ints, 
              main = bquote(paste("Basis method: location intercept by county")), 
              midpoint = 0)
map.ga.ggplot(Y = betamu.time, 
              main = bquote(paste("Basis method: location time coefficient by county")), 
              midpoint = 0)
map.ga.ggplot(Y = betamu.t2nd, 
              main = bquote(paste("Basis method: location time"^2, " coefficient by county")), 
              midpoint = 0)

map.ga.ggplot(Y = betalogsig.ints, 
              main = bquote(paste("Basis method log(", sigma, "): intercept by county")),
              midpoint = 0)
map.ga.ggplot(Y = betalogsig.time, 
              main = bquote(paste("Basis method log(", sigma, "): time coefficient by county")),
              midpoint = 0)
map.ga.ggplot(Y = betalogsig.t2nd, 
              main = bquote(paste("Basis method log(", sigma, "): time"^2, " coefficient by county")),
              midpoint = 0)
