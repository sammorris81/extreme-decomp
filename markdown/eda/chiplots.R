library(fields)
library(evd)
load("../../code/analysis/fire/georgia_preprocess/chi.RData")
load("../../code/analysis/fire/georgia_preprocess/fire_data.RData")
image.plot(1:ncol(Y), 1:ncol(Y), chi, main = "estimated chi")
image.plot(1:ncol(Y), 1:ncol(Y), ADJ, main = "adjacency")

# sort the adjacency matrix
order <- rep(NA, ncol(Y))
order[1] <- 1
for (i in 1:length(order)) {
  this <- order[i]
  neighbors <- which(ADJ[this, ] == 1)
  for (k in 1:length(neighbors)) {
    if (!(neighbors[k] %in% unique(order))) {
      this <- min(which(is.na(order)))
      order[this] <- neighbors[k]
    }
  }
}

chi.sort <- matrix(NA, ncol(Y), ncol(Y))
for (i in 1:ncol(Y)) {
  for (j in i:ncol(Y)) {
    pair.i <- order[i]
    pair.j <- order[j]
    chi.sort[i, j] <- chi.sort[j, i] <- chi[pair.i, pair.j]
  }
}

par(mfrow=c(1, 2))
chi.sort <- ifelse(chi.sort < 0, 0, chi.sort)
image.plot(1:ncol(Y), 1:ncol(Y), chi.sort, main = "estimated chi")
image.plot(1:ncol(Y), 1:ncol(Y), 2 - chi.sort, main = "estimated EC")

# lag-1 dependence
par(mfrow=c(2, 2))
nyears <- nrow(Y)
plot(Y[-nyears, 1], Y[-1, 1], main = "year total: lag 1, site 1",
     xlab = bquote(y[t]), ylab = bquote(y[t + 1]))
chiplot(cbind(Y[-nyears, 1], Y[-1, 1]), which = 1, ylim1 = c(-1, 1),
        main1 = bquote(paste(chi, "-plot: lag 1, site 1")))
plot(Y[-nyears, 2], Y[-1, 2], main = "year total: lag 1, site 2",
     xlab = bquote(y[t]), ylab = bquote(y[t + 1]))
chiplot(cbind(Y[-nyears, 2], Y[-1, 2]), which = 1, ylim1 = c(-1, 1),
        main1 = bquote(paste(chi, "-plot: lag 1, site 2")))

# Plots for months
nmonths <- nrow(y)
par(mfrow = c(2, 2))
plot(y[-nmonths, 1], y[-1, 1], main = "monthly total: lag 1, site 1",
     xlab = bquote(y[t]), ylab = bquote(y[t + 1]))
chiplot(cbind(y[-nmonths, 1], y[-1, 1]), which = 1, ylim1 = c(0, 1),
        main1 = bquote(paste(chi,"-plot: lag 1, site 1")))
plot(y[-nmonths, 2], y[-1, 2], main = "monthly total: lag 1, site 2",
     xlab = bquote(y[t]), ylab = bquote(y[t + 1]))
chiplot(cbind(y[-nmonths, 2], y[-1, 2]), which = 1, ylim1 = c(0, 1),
        main1 = bquote(paste(chi,"-plot: lag 1, site 2")))

plot(y[-nmonths, 3], y[-1, 3], main = "monthly total: lag 1, site 3",
     xlab = bquote(y[t]), ylab = bquote(y[t + 1]))
chiplot(cbind(y[-nmonths, 3], y[-1, 3]), which = 1, ylim1 = c(0, 1),
        main1 = bquote(paste(chi,"-plot: lag 1, site 3")))
plot(y[-nmonths, 4], y[-1, 4], main = "monthly total: lag 1, site 4",
     xlab = bquote(y[t]), ylab = bquote(y[t + 1]))
chiplot(cbind(y[-nmonths, 4], y[-1, 4]), which = 1, ylim1 = c(0, 1),
        main1 = bquote(paste(chi,"-plot: lag 1, site 4")))

plot(y[-nmonths, 5], y[-1, 5], main = "monthly total: lag 1, site 5",
     xlab = bquote(y[t]), ylab = bquote(y[t + 1]))
chiplot(cbind(y[-nmonths, 5], y[-1, 5]), which = 1, ylim1 = c(0, 1),
        main1 = bquote(paste(chi,"-plot: lag 1, site 5")))
plot(y[-nmonths, 6], y[-1, 6], main = "monthly total: lag 1, site 6",
     xlab = bquote(y[t]), ylab = bquote(y[t + 1]))
chiplot(cbind(y[-nmonths, 6], y[-1, 6]), which = 1, ylim1 = c(0, 1),
        main1 = bquote(paste(chi,"-plot: lag 1, site 6")))

all.sites <- cbind(y[-nmonths, 1], y[-1, 1])
for (i in 2:ncol(y)) {
  all.sites <- rbind(all.sites, cbind(y[-nmonths, i], y[-1, i]))
}
par(mfrow = c(1, 2))
plot(all.sites[, 1], all.sites[, 2], main = "monthly total: lag 1, all sites",
     xlab = bquote(y[t]), ylab = bquote(y[t + 1]))
chiplot(all.sites, which = 1, ylim1 = c(0, 1),
        main1 = bquote(paste(chi, "-plot: lag 1, all sites")))