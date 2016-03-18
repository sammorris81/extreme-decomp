library(colorspace)
load(file = "../../code/analysis/fire/georgia_preprocess/fire_data.RData")

# do a rank transformation of the data
# matplot makes one line per column of the dataset so we want Y as
# nt x ns
ns <- ncol(Y)
nt <- nrow(Y)

rank.transform <- Y
for (i in 1:ns) {
  rank.transform[, i] <- rank(Y[, i]) / (nt + 1)
}

matplot(rank.transform, type = "b", lty = 2, pch = 1,
        xaxt = "n", xlab = "year", ylab = "quantile",
        main = "Non-censored")
axis(1, at = seq(1, 49), labels = rownames(Y))

rank.transform.censor <- rank.transform
rank.transform.censor[rank.transform < 0.90] <- 0.90

matplot(rank.transform.censor, type = "b", lty = 2, pch = 1,
        xaxt = "n", xlab = "year", ylab = "quantile",
        main = "Censored at q(0.90)")
axis(1, at = seq(1, 49), labels = rownames(Y))

mn.noncensored <- apply(rank.transform, 1, mean)
plot(rownames(rank.transform), mn.noncensored, type = "b",
     xlab = "year", ylab = "mean of quantiles",
     main = "average quantile (no censoring)")

sd.noncensored <- apply(rank.transform, 1, sd)
plot(rownames(rank.transform), sd.noncensored, type = "b",
     xlab = "year", ylab = "log(std dev of quantiles)",
     main = "log(sd) of quantiles (no censoring)")

mn.censored <- apply(rank.transform.censor, 1, mean)
plot(rownames(rank.transform), mn.censored, type = "b",
     xlab = "year", ylab = "mean of quantiles",
     main = "average quantile (censored > 0.90)")

sd.censored    <- apply(rank.transform.censor, 1, sd)
plot(rownames(rank.transform), log(sd.censored), type = "b",
     xlab = "year", ylab = "log(std dev of quantiles > 0.90)",
     main = "log(sd) of quantiles (censored > 0.90)")

#### plot the data ####
library(fields)
load(file = "../../code/analysis/fire/georgia_preprocess/georgia_map.RData")
d <- rdist(cents)
diag(d) <- 0
n <- nrow(cents)
thresh <- rep(0, ns)

# standardize the locations
s <- cents
s[, 1] <- (s[, 1] - min(s[, 1])) / diff(range(s[, 1]))
s[, 2] <- (s[, 2] - min(s[, 2])) / diff(range(s[, 2]))
neighbors <- 5
d <- rdist(s)
diag(d) <- 0
for (i in 1:ns) {
  these <- order(d[i, ])[2:(neighbors + 1)]  # the closest is always site i
  thresh[i] <- quantile(Y[, these], probs = 0.90, na.rm = TRUE)
}
thresh <- t(matrix(thresh, ns, nt))  # note this is to keep consistent with Y
matplot(log(Y), pch = 1, type = "b", main = "all counties", xaxt = "n")
axis(1, at = seq(1, 49), labels = rownames(Y))

par(mfrow = c(2, 2))
for (i in 1:16) {
  start <- (i - 1) * 10 + 1
  end   <- min(i * 10, ns)
  matplot(log(Y[, start:end]), pch = 1, type = "b",
          main = paste("counties ", start, " -- ", end),
          xaxt = "n", ylab = "log(Y)", ylim = c(-5, 12),
          col = "dodgerblue1", lty = 3, add = FALSE)
  for (j in start:end) {
    these <- which(Y[, j] > thresh[, j])
    points(these, log(Y[these, j]), pch = 19, col = "firebrick1")
  }
  axis(1, at = seq(1, 49), labels = rownames(Y))
  if (i %% 4 == 0) {
    dev.print(device = pdf, file = paste("spag-", i / 4, ".pdf", sep = ""))
  }
}

for (i in 1:ns) {
  these <- which(Y[, i] > thresh[, i])
  if (i == 1) {
    plot(these, log(Y[these, i]), pch = 19, col = "firebrick1", type = "b",
         xaxt = "n", xlab = "", xlim = c(1, nt),
         ylim = c(3, 12), ylab = "log(Y)", lty = 3,
         main = "observations over spatial threshold")
    axis(1, at = seq(1, nt), labels = rownames(Y))
  } else {
    lines(these, log(Y[these, i]), pch = 19, col = "firebrick1",
          type = "b", lty = 3)
  }
}

#### plot time series for 10 counties ####
set.seed(1)
counties <- sample(1:ns, 10, replace = FALSE)
colors <- rainbow_hcl(n = length(counties))
# counties <- c(counties, which(Y == max(Y), arr.ind = TRUE)[2])
# counties <- c(counties, which(Y == min(Y), arr.ind = TRUE)[2])
matplot(log(Y[, counties]), type = "l",
        col = colors, cex = 1.25, lty = 1,
        ylab = "log(Y)", xaxt = "n",
        main = paste("time series at ", length(counties),
                     " randomly selected counties", sep = ""))
axis(1, at = seq(1, nt), labels = rownames(Y))
# for (i in 1:length(counties)) {
#   this.county <- counties[i]
#   these <- which(Y[, this.county] > thresh[, this.county])
#   points(these, log(Y[these, this.county]), pch = 19,
#          col = colors[i], cex = 1.25, lty = 3)
# }

#### plot time series for 25 counties ####
set.seed(1)
counties <- sample(1:ns, 25, replace = FALSE)
colors <- rainbow_hcl(n = length(counties))
# counties <- c(counties, which(Y == max(Y), arr.ind = TRUE)[2])
# counties <- c(counties, which(Y == min(Y), arr.ind = TRUE)[2])
par(mar = c(4.1, 4.1, 0.6, 0.6))
matplot(log(Y[, counties]), type = "l",
        col = colors, cex.lab = 1.5, cex.axis = 1.5, lty = 1,
        ylab = "log(acres burned)",
        xaxt = "n", xlab = "Year" #, xlim = c(0, (nt + 1))
        #main = paste("time series at ", length(counties),
        #             " randomly selected counties", sep = "")
        )

years <- rownames(Y)
years <- c("1995", years, "2015")
axis(1, at = seq(0, (nt + 1), by = 5), labels = seq(1965, 2015, by = 5),
     cex.axis = 1.5)
# for (i in 1:length(counties)) {
#   this.county <- counties[i]
#   these <- which(Y[, this.county] > thresh[, this.county])
#   points(these, log(Y[these, this.county]), pch = 19,
#          col = colors[i], cex = 1.25, lty = 3)
# }

#### plot exceedances ####
set.seed(1)
counties <- sample(1:ns, 80, replace = FALSE)
colors <- rainbow_hcl(n = length(counties))
for (i in 1:length(counties)) {
  this.county <- counties[i]
  these <- which(Y[, this.county] > thresh[, this.county])
  if (i == 1) {
    plot(these, log(Y[these, this.county]), type = "l", pch = 19,
         col = colors[i], xaxt = "n", xlab = "", xlim = c(1, nt),
         ylim = c(3, 12), ylab = "log(Y)", lty = 1, lwd = 1.5,
         main = paste("observations over spatial threshold for ",
                      length(counties), " counties", sep = ""))
  } else {
    lines(these, log(Y[these, this.county]), pch = 19,
          col = colors[i], lwd = 1.5, lty = 1)
  }
}
axis(1, at = seq(1, nt), labels = rownames(Y))

# mean residual life plot - first center and scale
library(ismev)
Y.adj <- Y <- t(Y)
qs <- apply(Y, 1, quantile, probs = c(0.25, 0.75))
for (i in 1:ns) {
  Y.adj[i, ] <- (Y[i, ] - median(Y[i, ])) / diff(qs[, i])
}

Y.adj <- Y.adj + min(Y.adj)

mrl.plot(Y.adj)
mrl.plot(Y.adj[1, ])
mrl.plot(Y.adj[2, ])

par(mfrow = c(1, 2))
mrl.plot(sqrt(Y))
mrl.plot(Y.adj, nint = 1000)
str(Y)
mean(Y)
range(Y)
mrl.plot(sqrt(Y), nint = 1000)
gpd.fitrange(Y[1, ], umin = -1.2, umax = 2, nint = 20)
gpd.fitrange(Y[1, ], umin = 10, umax = 1900, nint = 20)


mrl.plot(Y.adj, nint = 100)
title("MRL plot, full range of Y.adj")

colors <- rainbow_hcl(n = 4)
mrl.plot(Y.adj, umin = quantile(Y.adj, probs = 0.4), 
         umax = quantile(Y.adj, probs = 0.999), nint = 100, 
         xlab = "Adjusted Y")
title("MRL plot, zoom range to q(0.999)")

x <- quantile(Y.adj, probs = 0.90)
abline(v = x, lty = 1, col = colors[1])
text(x = x + 0.3, y = 60, labels = 1)

x <- quantile(Y.adj, probs = 0.95)
abline(v = x, lty = 1)
text(x = x, y = 70, labels = "q(0.95)", pos = 4, lwd = 1.25)

x <- quantile(Y.adj, probs = 0.98)
abline(v = x, lty = 1, col = colors[3])
text(x = x + 0.3, y = 60, labels = 3)

x <- quantile(Y.adj, probs = 0.99)
abline(v = x, lty = 1, col = colors[4])
text(x = x + 0.3, y = 60, labels = 4)
