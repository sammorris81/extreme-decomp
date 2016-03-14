source(file = "./package_load.R", chdir = T)
library(rapportools)  # get camelcase function
######## Linear time trend ########

#### Basis method
method <- "basis" # using kern for the results from abba
L      <- 15  # will be using this to get basis functions for covariates
time.trend <- "linear"
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

ns <- nrow(Y)
nt <- ncol(Y)
np <- 2 + L * 2  # for a single year (int, t, B1...BL, t * (B1...BL))

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

title <- bquote(paste("Method: ", .(tocamel(method, upper=TRUE)), 
                      " (", .(time.trend), " time trend)"))
legend.title <- bquote(paste(mu, ": ", beta["int"]))
map.ga.ggplot(Y = betamu.ints, main = title, fill.legend = legend.title, 
              midpoint = 0)
pdf.name <- paste("plots/", method, "-mu-int-1st.pdf", sep = "")
dev.print(device = pdf, pdf.name)

title <- bquote(paste("Method: ", .(tocamel(method, upper=TRUE)), 
                      " (", .(time.trend), " time trend)"))
legend.title <- bquote(paste(mu, ": ", beta["time"]))
map.ga.ggplot(Y = betamu.time, main = title, fill.legend = legend.title, 
              midpoint = 0)
pdf.name <- paste("plots/", method, "-mu-t-1st.pdf", sep = "")
dev.print(device = pdf, pdf.name)

title <- bquote(paste("Method: ", .(tocamel(method, upper=TRUE)), 
                      " (", .(time.trend), " time trend)"))
legend.title <- bquote(paste("log(", sigma, "): ", beta["int"]))
map.ga.ggplot(Y = betalogsig.ints, main = title, fill.legend = legend.title,
              midpoint = 0)
pdf.name <- paste("plots/", method, "-logsig-int-1st.pdf", sep = "")
dev.print(device = pdf, pdf.name)

title <- bquote(paste("Method: ", .(tocamel(method, upper=TRUE)), 
                      " (", .(time.trend), " time trend)"))
legend.title <- bquote(paste("log(", sigma, "): ", beta["time"]))
map.ga.ggplot(Y = betalogsig.time, main = title, fill.legend = legend.title,
              midpoint = 0)
pdf.name <- paste("plots/", method, "-logsig-t-1st.pdf", sep = "")
dev.print(device = pdf, pdf.name)

# look at the posterior distribution of the slope terms
niters <- nrow(fit$beta1)
mu.int.post  <- matrix(0, niters, ns)  # niters x ns
mu.time.post <- matrix(0, niters, ns)  # niters x ns
logsig.int.post  <- matrix(0, niters, ns)  # niters x ns
logsig.time.post <- matrix(0, niters, ns)  # niters x ns

for (i in 1:niters) {
  mu.int.post[i, ]  <- fit$beta1[i, 1] + B.est %*% fit$beta1[i, 3:(L + 2)]
  mu.time.post[i, ] <- fit$beta1[i, 2] + B.est %*% fit$beta1[i, (L + 3):np]
  logsig.int.post[i, ]  <- fit$beta2[i, 1] + B.est %*% fit$beta2[i, 3:(L + 2)]
  logsig.time.post[i, ] <- fit$beta2[i, 2] + B.est %*% fit$beta2[i, (L + 3):np]
}

# county by county P(slope < 0)
p.mu.time.neg <- apply(mu.time.post < 0, 2, mean)
p.logsig.time.neg <- apply(logsig.time.post < 0, 2, mean)

title <- bquote(paste("Method: ", .(tocamel(method, upper=TRUE)), 
                      " (", .(time.trend), " time trend)"))
legend.title <- bquote(paste(mu, ": P(", beta["time"], "< 0)"))
map.ga.ggplot(Y = p.mu.time.neg, main = title, fill.legend = legend.title,
              mid = 0.5, limits = c(0, 1))
pdf.name <- paste("plots/", method, "-prob-beta-mu-t-neg-1st.pdf", sep = "")
dev.print(device = pdf, pdf.name)

title <- bquote(paste("Method: ", .(tocamel(method, upper=TRUE)), 
                      " (", .(time.trend), " time trend)"))
legend.title <- bquote(paste("log(", sigma, "): P(", beta["time"], "< 0)"))
map.ga.ggplot(Y = p.logsig.time.neg, main = title, fill.legend = legend.title, 
              mid = 0.5, limits = c(0, 1))
pdf.name <- paste("plots/", method, "-prob-beta-logsig-t-neg-1st.pdf", sep = "")
dev.print(device = pdf, pdf.name)

#### Kernel method
method <- "kern" # using kern for the results from abba
L      <- 15  # will be using this to get basis functions for covariates
time.trend <- "linear"
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

ns <- nrow(Y)
nt <- ncol(Y)
np <- 2 + L * 2  # for a single year (int, t, B1...BL, t * (B1...BL))

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

title <- bquote(paste("Method: ", .(tocamel(method, upper=TRUE)), 
                      " (", .(time.trend), " time trend)"))
legend.title <- bquote(paste(mu, ": ", beta["int"]))
map.ga.ggplot(Y = betamu.ints, main = title, fill.legend = legend.title, 
              midpoint = 0)
pdf.name <- paste("plots/", method, "-mu-int-1st.pdf", sep = "")
dev.print(device = pdf, pdf.name)

title <- bquote(paste("Method: ", .(tocamel(method, upper=TRUE)), 
                      " (", .(time.trend), " time trend)"))
legend.title <- bquote(paste(mu, ": ", beta["time"]))
map.ga.ggplot(Y = betamu.time, main = title, fill.legend = legend.title, 
              midpoint = 0)
pdf.name <- paste("plots/", method, "-mu-t-1st.pdf", sep = "")
dev.print(device = pdf, pdf.name)

title <- bquote(paste("Method: ", .(tocamel(method, upper=TRUE)), 
                      " (", .(time.trend), " time trend)"))
legend.title <- bquote(paste("log(", sigma, "): ", beta["int"]))
map.ga.ggplot(Y = betalogsig.ints, main = title, fill.legend = legend.title, 
              midpoint = 0)
pdf.name <- paste("plots/", method, "-logsig-int-1st.pdf", sep = "")
dev.print(device = pdf, pdf.name)

title <- bquote(paste("Method: ", .(tocamel(method, upper=TRUE)), 
                      " (", .(time.trend), " time trend)"))
legend.title <- bquote(paste("log(", sigma, "): ", beta["time"]))
map.ga.ggplot(Y = betalogsig.time, main = title, fill.legend = legend.title, 
              midpoint = 0)
pdf.name <- paste("plots/", method, "-logsig-t-1st.pdf", sep = "")
dev.print(device = pdf, pdf.name)

# look at the posterior distribution of the slope terms
niters <- nrow(fit$beta1)
mu.int.post  <- matrix(0, niters, ns)  # niters x ns
mu.time.post <- matrix(0, niters, ns)  # niters x ns
logsig.int.post  <- matrix(0, niters, ns)  # niters x ns
logsig.time.post <- matrix(0, niters, ns)  # niters x ns

for (i in 1:niters) {
  mu.int.post[i, ]  <- fit$beta1[i, 1] + B.est %*% fit$beta1[i, 3:(L + 2)]
  mu.time.post[i, ] <- fit$beta1[i, 2] + B.est %*% fit$beta1[i, (L + 3):np]
  logsig.int.post[i, ]  <- fit$beta2[i, 1] + B.est %*% fit$beta2[i, 3:(L + 2)]
  logsig.time.post[i, ] <- fit$beta2[i, 2] + B.est %*% fit$beta2[i, (L + 3):np]
}

# county by county P(slope < 0)
p.mu.time.neg <- apply(mu.time.post < 0, 2, mean)
p.logsig.time.neg <- apply(logsig.time.post < 0, 2, mean)

title <- bquote(paste("Method: ", .(tocamel(method, upper=TRUE)), 
                      " (", .(time.trend), " time trend)"))
legend.title <- bquote(paste(mu, ": P(", beta["time"], "< 0)"))
map.ga.ggplot(Y = p.mu.time.neg, main = title, fill.legend = legend.title,
              mid = 0.5, limits = c(0, 1))
pdf.name <- paste("plots/", method, "-prob-beta-mu-t-neg-1st.pdf", sep = "")
dev.print(device = pdf, pdf.name)

title <- bquote(paste("Method: ", .(tocamel(method, upper=TRUE)), 
                      " (", .(time.trend), " time trend)"))
legend.title <- bquote(paste("log(", sigma, "): P(", beta["time"], "< 0)"))
map.ga.ggplot(Y = p.logsig.time.neg, main = title, fill.legend = legend.title,
              mid = 0.5, limits = c(0, 1))
pdf.name <- paste("plots/", method, "-prob-beta-logsig-t-neg-1st.pdf", sep = "")
dev.print(device = pdf, pdf.name)

######## Quadratic time trend ########

#### Basis method
method <- "basis" # using kern for the results from abba
L      <- 15  # will be using this to get basis functions for covariates
time.trend <- "quadratic"
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

title <- bquote(paste("Method: ", .(tocamel(method, upper=TRUE)), 
                      " (", .(time.trend), " time trend)"))
legend.title <- bquote(paste(mu, ": ", beta["int"]))
map.ga.ggplot(Y = betamu.ints, main = title, fill.legend = legend.title,
              midpoint = 0)
pdf.name <- paste("plots/", method, "-mu-int-2nd.pdf", sep = "")
dev.print(device = pdf, pdf.name)

title <- bquote(paste("Method: ", .(tocamel(method, upper=TRUE)), 
                      " (", .(time.trend), " time trend)"))
legend.title <- bquote(paste(mu, ": ", beta["time"]))
map.ga.ggplot(Y = betamu.time, main = title, fill.legend = legend.title, 
              midpoint = 0)
pdf.name <- paste("plots/", method, "-mu-t-2nd.pdf", sep = "")
dev.print(device = pdf, pdf.name)

title <- bquote(paste("Method: ", .(tocamel(method, upper=TRUE)), 
                      " (", .(time.trend), " time trend)"))
legend.title <- bquote(paste(mu, ": ", beta["time"^2]))
map.ga.ggplot(Y = betamu.t2nd, main = title, fill.legend = legend.title, 
              midpoint = 0)
pdf.name <- paste("plots/", method, "-mu-tsq-2nd.pdf", sep = "")
dev.print(device = pdf, pdf.name)

title <- bquote(paste("Method: ", .(tocamel(method, upper=TRUE)), 
                      " (", .(time.trend), " time trend)"))
legend.title <- bquote(paste("log(", sigma, "): ", beta["int"]))
map.ga.ggplot(Y = betalogsig.ints, main = title, fill.legend = legend.title,
              midpoint = 0)
pdf.name <- paste("plots/", method, "-logsig-int-2nd.pdf", sep = "")
dev.print(device = pdf, pdf.name)

title <- bquote(paste("Method: ", .(tocamel(method, upper=TRUE)), 
                      " (", .(time.trend), " time trend)"))
legend.title <- bquote(paste("log(", sigma, "): ", beta["time"]))
map.ga.ggplot(Y = betalogsig.time, main = title, fill.legend = legend.title, 
              midpoint = 0)
pdf.name <- paste("plots/", method, "-logsig-t-2nd.pdf", sep = "")
dev.print(device = pdf, pdf.name)

title <- bquote(paste("Method: ", .(tocamel(method, upper=TRUE)), 
                      " (", .(time.trend), " time trend)"))
legend.title <- bquote(paste("log(", sigma, "): ", beta["time"^2]))
map.ga.ggplot(Y = betalogsig.t2nd, main = title, fill.legend = legend.title,
              midpoint = 0)
pdf.name <- paste("plots/", method, "-logsig-tsq-2nd.pdf", sep = "")
dev.print(device = pdf, pdf.name)

# look at the posterior distribution of the slope terms
niters <- nrow(fit$beta1)
mu.int.post  <- matrix(0, niters, ns)  # niters x ns
mu.time.post <- matrix(0, niters, ns)  # niters x ns
mu.t2nd.post <- matrix(0, niters, ns)  # niters x ns
logsig.int.post  <- matrix(0, niters, ns)  # niters x ns
logsig.time.post <- matrix(0, niters, ns)  # niters x ns
logsig.t2nd.post <- matrix(0, niters, ns)  # niters x ns

for (i in 1:niters) {
  mu.int.post[i, ]  <- fit$beta1[i, 1] + B.est %*% fit$beta1[i, 4:(L + 3)]
  mu.time.post[i, ] <- fit$beta1[i, 2] + B.est %*% fit$beta1[i, (L + 4):(2 * L + 3)]
  mu.t2nd.post[i, ] <- fit$beta1[i, 3] + B.est %*% fit$beta1[i, (2 * L + 4):np]
  logsig.int.post[i, ]  <- fit$beta2[i, 1] + B.est %*% fit$beta2[i, 3:(L + 2)]
  logsig.time.post[i, ] <- fit$beta2[i, 2] + B.est %*% fit$beta2[i, (L + 4):(2 * L + 3)]
  logsig.t2nd.post[i, ] <- fit$beta2[i, 3] + B.est %*% fit$beta2[i, (2 * L + 4):np]
}

# county by county P(slope < 0)
p.mu.time.neg <- apply(mu.time.post < 0, 2, mean)
p.mu.t2nd.neg <- apply(mu.t2nd.post < 0, 2, mean)
p.logsig.time.neg <- apply(logsig.time.post < 0, 2, mean)
p.logsig.t2nd.neg <- apply(logsig.t2nd.post < 0, 2, mean)

title <- bquote(paste("Method: ", .(tocamel(method, upper=TRUE)), 
                      " (", .(time.trend), " time trend)"))
legend.title <- bquote(paste(mu, ": P(", beta["time"], "< 0)"))
map.ga.ggplot(Y = p.mu.time.neg, main = title, fill.legend = legend.title,
              mid = 0.5, limits = c(0, 1))
pdf.name <- paste("plots/", method, "-prob-beta-mu-t-neg-2nd.pdf", sep = "")
dev.print(device = pdf, pdf.name)

title <- bquote(paste("Method: ", .(tocamel(method, upper=TRUE)), 
                      " (", .(time.trend), " time trend)"))
legend.title <- bquote(paste(mu, ": P(", beta["time"^2], "< 0)"))
map.ga.ggplot(Y = p.mu.t2nd.neg, main = title, fill.legend = legend.title, 
              mid = 0.5, limits = c(0, 1))
pdf.name <- paste("plots/", method, "-prob-beta-mu-tsq-neg-2nd.pdf", sep = "")
dev.print(device = pdf, pdf.name)

title <- bquote(paste("Method: ", .(tocamel(method, upper=TRUE)), 
                      " (", .(time.trend), " time trend)"))
legend.title <- bquote(paste("log(", sigma, "): P(", beta["time"], "< 0)"))
map.ga.ggplot(Y = p.logsig.time.neg, main = title, fill.legend = legend.title, 
              mid = 0.5, limits = c(0, 1))
pdf.name <- paste("plots/", method, "-prob-beta-logsig-t-neg-2nd.pdf", sep = "")
dev.print(device = pdf, pdf.name)

title <- bquote(paste("Method: ", .(tocamel(method, upper=TRUE)), 
                      " (", .(time.trend), " time trend)"))
legend.title <- bquote(paste("log(", sigma, "): P(", beta["time"^2], "< 0)"))
map.ga.ggplot(Y = p.logsig.t2nd.neg, main = title, fill.legend = legend.title, 
              mid = 0.5, limits = c(0, 1))
pdf.name <- paste("plots/", method, "-prob-beta-logsig-tsq-neg-2nd.pdf", sep = "")
dev.print(device = pdf, pdf.name)

#### Kernel method
method <- "kern" # using kern for the results from abba
L      <- 15  # will be using this to get basis functions for covariates
time.trend <- "quadratic"
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

ns <- nrow(Y)
nt <- ncol(Y)
np <- 3 + L * 3  # for a single year (int, t, t^2, B1...BL, t * (B1...BL), 
                 # t^2 * (B1...BL))

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

title <- bquote(paste("Method: ", .(tocamel(method, upper=TRUE)), 
                      " (", .(time.trend), " time trend)"))
legend.title <- bquote(paste(mu, ": ", beta["int"]))
map.ga.ggplot(Y = betamu.ints, main = title, fill.legend = legend.title, 
              midpoint = 0)
pdf.name <- paste("plots/", method, "-mu-int-2nd.pdf", sep = "")
dev.print(device = pdf, pdf.name)

title <- bquote(paste("Method: ", .(tocamel(method, upper=TRUE)), 
                      " (", .(time.trend), " time trend)"))
legend.title <- bquote(paste(mu, ": ", beta["time"]))
map.ga.ggplot(Y = betamu.time, main = title, fill.legend = legend.title, 
              midpoint = 0)
pdf.name <- paste("plots/", method, "-mu-t-2nd.pdf", sep = "")
dev.print(device = pdf, pdf.name)

title <- bquote(paste("Method: ", .(tocamel(method, upper=TRUE)), 
                      " (", .(time.trend), " time trend)"))
legend.title <- bquote(paste(mu, ": ", beta["time"^2]))
map.ga.ggplot(Y = betamu.t2nd, main = title, fill.legend = legend.title, 
              midpoint = 0)
pdf.name <- paste("plots/", method, "-mu-tsq-2nd.pdf", sep = "")
dev.print(device = pdf, pdf.name)

title <- bquote(paste("Method: ", .(tocamel(method, upper=TRUE)), 
                      " (", .(time.trend), " time trend)"))
legend.title <- bquote(paste("log(", sigma, "): ", beta["int"]))
map.ga.ggplot(Y = betalogsig.ints, main = title, fill.legend = legend.title, 
              midpoint = 0)
pdf.name <- paste("plots/", method, "-logsig-int-2nd.pdf", sep = "")
dev.print(device = pdf, pdf.name)

title <- bquote(paste("Method: ", .(tocamel(method, upper=TRUE)), 
                      " (", .(time.trend), " time trend)"))
legend.title <- bquote(paste("log(", sigma, "): ", beta["time"]))
map.ga.ggplot(Y = betalogsig.time, main = title, fill.legend = legend.title, 
              midpoint = 0)
pdf.name <- paste("plots/", method, "-logsig-t-2nd.pdf", sep = "")
dev.print(device = pdf, pdf.name)

title <- bquote(paste("Method: ", .(tocamel(method, upper=TRUE)), 
                      " (", .(time.trend), " time trend)"))
legend.title <- bquote(paste("log(", sigma, "): ", beta["time"^2]))
map.ga.ggplot(Y = betalogsig.t2nd, main = title, fill.legend = legend.title,
              midpoint = 0)
pdf.name <- paste("plots/", method, "-logsig-tsq-2nd.pdf", sep = "")
dev.print(device = pdf, pdf.name)

# look at the posterior distribution of the slope terms
niters <- nrow(fit$beta1)
mu.int.post  <- matrix(0, niters, ns)  # niters x ns
mu.time.post <- matrix(0, niters, ns)  # niters x ns
mu.t2nd.post <- matrix(0, niters, ns)  # niters x ns
logsig.int.post  <- matrix(0, niters, ns)  # niters x ns
logsig.time.post <- matrix(0, niters, ns)  # niters x ns
logsig.t2nd.post <- matrix(0, niters, ns)  # niters x ns

for (i in 1:niters) {
  mu.int.post[i, ]  <- fit$beta1[i, 1] + B.est %*% fit$beta1[i, 4:(L + 3)]
  mu.time.post[i, ] <- fit$beta1[i, 2] + B.est %*% fit$beta1[i, (L + 4):(2 * L + 3)]
  mu.t2nd.post[i, ] <- fit$beta1[i, 3] + B.est %*% fit$beta1[i, (2 * L + 4):np]
  logsig.int.post[i, ]  <- fit$beta2[i, 1] + B.est %*% fit$beta2[i, 3:(L + 2)]
  logsig.time.post[i, ] <- fit$beta2[i, 2] + B.est %*% fit$beta2[i, (L + 4):(2 * L + 3)]
  logsig.t2nd.post[i, ] <- fit$beta2[i, 3] + B.est %*% fit$beta2[i, (2 * L + 4):np]
}

# county by county P(slope < 0)
p.mu.time.neg <- apply(mu.time.post < 0, 2, mean)
p.mu.t2nd.neg <- apply(mu.t2nd.post < 0, 2, mean)
p.logsig.time.neg <- apply(logsig.time.post < 0, 2, mean)
p.logsig.t2nd.neg <- apply(logsig.t2nd.post < 0, 2, mean)

title <- bquote(paste("Method: ", .(tocamel(method, upper=TRUE)), 
                      " (", .(time.trend), " time trend)"))
legend.title <- bquote(paste(mu, ": P(", beta["time"], "< 0)"))
map.ga.ggplot(Y = p.mu.time.neg, main = title, fill.legend = legend.title, 
              mid = 0.5, limits = c(0, 1))
pdf.name <- paste("plots/", method, "-prob-beta-mu-t-neg-2nd.pdf", sep = "")
dev.print(device = pdf, pdf.name)

title <- bquote(paste("Method: ", .(tocamel(method, upper=TRUE)), 
                      " (", .(time.trend), " time trend)"))
legend.title <- bquote(paste(mu, ": P(", beta["time"^2], "< 0)"))
map.ga.ggplot(Y = p.mu.t2nd.neg, main = title, fill.legend = legend.title,
              mid = 0.5, limits = c(0, 1))
pdf.name <- paste("plots/", method, "-prob-beta-mu-tsq-neg-2nd.pdf", sep = "")
dev.print(device = pdf, pdf.name)

title <- bquote(paste("Method: ", .(tocamel(method, upper=TRUE)), 
                      " (", .(time.trend), " time trend)"))
legend.title <- bquote(paste("log(", sigma, "): P(", beta["time"], "< 0)"))
map.ga.ggplot(Y = p.logsig.time.neg, main = title, fill.legend = legend.title,
              mid = 0.5, limits = c(0, 1))
pdf.name <- paste("plots/", method, "-prob-beta-logsig-t-neg-2nd.pdf", sep = "")
dev.print(device = pdf, pdf.name)

title <- bquote(paste("Method: ", .(tocamel(method, upper=TRUE)), 
                      " (", .(time.trend), " time trend)"))
legend.title <- bquote(paste("log(", sigma, "): P(", beta["time"^2], "< 0)"))
map.ga.ggplot(Y = p.logsig.t2nd.neg, main = title, fill.legend = legend.title, 
              mid = 0.5, limits = c(0, 1))
pdf.name <- paste("plots/", method, "-prob-beta-logsig-tsq-neg-2nd.pdf", sep = "")
dev.print(device = pdf, pdf.name)
