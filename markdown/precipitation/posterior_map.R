rm(list = ls())
source("./package_load.R", chdir = TRUE)
library(evd)
library(gridExtra)
# trying to see if there's a difference in the estimates for gsk vs ebf

load("precip_preprocess.RData")
nx <- length(unique(s[, 1]))
ny <- length(unique(s[, 2]))

ns <- nrow(Y)
nt <- ncol(Y) / 2
niters <- 10000

time.1 <- (1 - nt / 2) / nt
time.t <- (nt - nt / 2) / nt
# mid <- "grey90"
mid <- "#FFFFFF"

#### Spatial plots of: ####
#  1. Posterior means of betas: 4
#  2. Posterior P(B_time > 0): 2
#  3. Change in q(90)

load("./cv-results/ebf-current-25-all.RData")
fit.ebf.cur <- fit
load("./cv-results/ebf-future-25-all.RData")
fit.ebf.fut <- fit
load("./cv-results/gsk-current-25-all.RData")
fit.gsk.cur <- fit
load("./cv-results/gsk-future-25-all.RData")
fit.gsk.fut <- fit

method <- "ebf"
this.fit <- fit.ebf.cur
beta.mu.int  <- apply(this.fit$beta.int[, , 1], 2, mean)
beta.mu.time <- apply(this.fit$beta.time[, , 1], 2, mean)
beta.ls.int  <- apply(this.fit$beta.int[, , 2], 2, mean)
beta.ls.time <- apply(this.fit$beta.time[, , 2], 2, mean)
prob.beta.mu.time.pos <- apply(this.fit$beta.time[, , 1] > 0, 2, mean)
prob.beta.ls.time.pos <- apply(this.fit$beta.time[, , 2] > 0, 2, mean)

main <- bquote(paste(.(toupper(method)), ": Current Posterior Mean of ",
                     beta["1, int"], sep = ""))
fill.legend <- bquote(beta["1, int"])
p <- map.heatmap(lat = s[, 2], lon = s[, 1], data = beta.mu.int,
            mainTitle = main, legendTitle = fill.legend, color_mid = mid)
ggsave(paste("./plots/precip-", method, "-cur-post-beta1int.pdf", sep = ""),
       p, width = 8, height = 8)

main <- bquote(paste(.(toupper(method)), ": Current Posterior Mean of ",
                     beta["1, time"], sep = ""))
fill.legend <- bquote(beta["1, time"])
p1.cur <- map.heatmap(lat = s[, 2], lon = s[, 1], data = beta.mu.time,
                      mainTitle = main, legendTitle = fill.legend,
                      color_mid = mid, midpoint = 0)

main <- bquote(paste(.(toupper(method)), ": Current Posterior Mean of ",
                     beta["2, time"], sep = ""))
fill.legend <- bquote(beta["2, time"])
p2.cur <- map.heatmap(lat = s[, 2], lon = s[, 1], data = beta.ls.time,
            mainTitle = main, legendTitle = fill.legend, color_mid = mid,
            midpoint = 0)

main <- bquote(paste(.(toupper(method)), ": Current Posterior P(", beta["1, time"],
                     " > 0)", sep = ""))
fill.legend <- bquote(paste("P(", beta["1, time"], " > 0)", sep = ""))
p3.cur <- map.heatmap(lat = s[, 2], lon = s[, 1], data = prob.beta.mu.time.pos,
                      mainTitle = main, legendTitle = fill.legend,
                      color_mid = mid, midpoint = 0.5)
# ggsave(paste("./plots/precip-", method, "-cur-post-beta1timepos.pdf", sep = ""),
#        p, width = 8, height = 8)

main <- bquote(paste(.(toupper(method)), ": Current Posterior P(", beta["2, time"],
                     " > 0)", sep = ""))
fill.legend <- bquote(paste("P(", beta["2, time"], " > 0)", sep = ""))
p4.cur <- map.heatmap(lat = s[, 2], lon = s[, 1], data = prob.beta.ls.time.pos,
                      mainTitle = main, legendTitle = fill.legend,
                      color_mid = mid, midpoint = 0.5)
# ggsave(paste("./plots/precip-", method, "-cur-post-beta2timepos.pdf", sep = ""),
#        p, width = 8, height = 8)

#### Get change in q.90 by county
mu.1  <- beta.mu.int + beta.mu.time * time.1
# mu.10 <- beta.mu.int + beta.mu.time * time.10
mu.t  <- beta.mu.int + beta.mu.time * time.t
ls.1  <- beta.ls.int + beta.ls.time * time.1
# ls.10 <- beta.ls.int + beta.ls.time * time.10
ls.t  <- beta.ls.int + beta.ls.time * time.t
xi    <- mean(fit$xi)
q.90.1 <- q.90.t <- rep(NA, ns)
for (i in 1:ns) {
  q.90.1[i] <- qgev(p = 0.90, loc = mu.1[i], scale = exp(ls.1[i]), shape = xi)
  # q.90.10[i] <- qgev(p = 0.90, loc = mu.10[i], scale = exp(ls.10[i]),
  #                    shape = xi)
  q.90.t[i] <- qgev(p = 0.90, loc = mu.t[i], scale = exp(ls.t[i]), shape = xi)
}
q.90.diff.nt <- q.90.t - q.90.1   # change over whole dataset
# q.90.diff.10 <- q.90.t - q.90.10  # change in last 10 years

main <- paste(toupper(method), ": Current Change in q(0.90)", sep = "")
fill.legend <- "Difference"
pdiff.ebf.cur <- map.heatmap(lat = s[, 2], lon = s[, 1], data = q.90.diff.nt,
                             mainTitle = main, legendTitle = fill.legend,
                             color_mid = mid, midpoint = 0)
# ggsave(paste("./plots/precip-", method, "-cur-q90diff.pdf", sep = ""),
#        p, width = 8, height = 8)

load(paste("./cv-results/", method, "-future-25-all.RData", sep = ""))
beta.mu.int  <- apply(fit$beta.int[, , 1], 2, mean)
beta.mu.time <- apply(fit$beta.time[, , 1], 2, mean)
beta.ls.int  <- apply(fit$beta.int[, , 2], 2, mean)
beta.ls.time <- apply(fit$beta.time[, , 2], 2, mean)
prob.beta.mu.time.pos <- apply(fit$beta.time[, , 1] > 0, 2, mean)
prob.beta.ls.time.pos <- apply(fit$beta.time[, , 2] > 0, 2, mean)

main <- bquote(paste(.(toupper(method)), ": Future Posterior Mean of ",
                     beta["1, int"], sep = ""))
fill.legend <- bquote(beta["1, int"])
p <- map.heatmap(lat = s[, 2], lon = s[, 1], data = beta.mu.int,
            mainTitle = main, legendTitle = fill.legend, color_mid = mid)
ggsave(paste("./plots/precip-", method, "-fut-post-beta1int.pdf", sep = ""),
       p, width = 8, height = 8)

main <- bquote(paste(.(toupper(method)), ": Future Posterior Mean of ",
                     beta["1, time"], sep = ""))
fill.legend <- bquote(beta["1, time"])
p1.fut <- map.heatmap(lat = s[, 2], lon = s[, 1], data = beta.mu.time,
                      mainTitle = main, legendTitle = fill.legend, color_mid = mid,
                      midpoint = 0)
# ggsave(paste("./plots/precip-", method, "-fut-post-beta1time.pdf", sep = ""),
#        p, width = 8, height = 8)

main <- bquote(paste(.(toupper(method)), ": Future Posterior Mean of ",
                     beta["2, time"], sep = ""))
fill.legend <- bquote(beta["2, time"])
p2.fut <- map.heatmap(lat = s[, 2], lon = s[, 1], data = beta.ls.time,
                      mainTitle = main, legendTitle = fill.legend,
                      color_mid = mid, midpoint = 0)
# ggsave(paste("./plots/precip-", method, "-fut-post-beta2time.pdf", sep = ""),
#        p, width = 8, height = 8)

main <- bquote(paste(.(toupper(method)), ": Future Posterior P(", beta["1, time"],
                     " > 0)", sep = ""))
fill.legend <- bquote(paste("P(", beta["1, time"], " > 0)", sep = ""))
p3.fut <- map.heatmap(lat = s[, 2], lon = s[, 1], data = prob.beta.mu.time.pos,
                      mainTitle = main, legendTitle = fill.legend,
                      color_mid = mid, midpoint = 0.5)
# ggsave(paste("./plots/precip-", method, "-cur-post-beta1timepos.pdf", sep = ""),
#        p, width = 8, height = 8)

main <- bquote(paste(.(toupper(method)), ": Future Posterior P(", beta["2, time"],
                     " > 0)", sep = ""))
fill.legend <- bquote(paste("P(", beta["2, time"], " > 0)", sep = ""))
p4.fut <- map.heatmap(lat = s[, 2], lon = s[, 1], data = prob.beta.ls.time.pos,
                      mainTitle = main, legendTitle = fill.legend,
                      color_mid = mid, midpoint = 0.5)
# ggsave(paste("./plots/precip-", method, "-cur-post-beta2timepos.pdf", sep = ""),
#        p, width = 8, height = 8)

#### Get change in q.90 by county
mu.1  <- beta.mu.int + beta.mu.time * time.1
# mu.10 <- beta.mu.int + beta.mu.time * time.10
mu.t  <- beta.mu.int + beta.mu.time * time.t
ls.1  <- beta.ls.int + beta.ls.time * time.1
# ls.10 <- beta.ls.int + beta.ls.time * time.10
ls.t  <- beta.ls.int + beta.ls.time * time.t
xi    <- mean(fit$xi)
q.90.1 <- q.90.t <- rep(NA, ns)
for (i in 1:ns) {
  q.90.1[i] <- qgev(p = 0.90, loc = mu.1[i], scale = exp(ls.1[i]), shape = xi)
  # q.90.10[i] <- qgev(p = 0.90, loc = mu.10[i], scale = exp(ls.10[i]),
  #                    shape = xi)
  q.90.t[i] <- qgev(p = 0.90, loc = mu.t[i], scale = exp(ls.t[i]), shape = xi)
}
q.90.diff.nt <- q.90.t - q.90.1   # change over whole dataset
# q.90.diff.10 <- q.90.t - q.90.10  # change in last 10 years

main <- paste(toupper(method), ": Future Change in q(0.90)", sep = "")
fill.legend <- "Difference"
pdiff.ebf.fut <- map.heatmap(lat = s[, 2], lon = s[, 1], data = q.90.diff.nt,
                             mainTitle = main, legendTitle = fill.legend,
                             color_mid = mid, midpoint = 0)
# ggsave(paste("./plots/precip-", method, "-fut-q90diff.pdf", sep = ""),
#        p, width = 8, height = 8)

layout.mtx <- matrix(1:4, nrow = 2, ncol = 2)
panel <- arrangeGrob(p1.cur, p1.fut, p2.cur, p2.fut, ncol = 2,
                     layout_matrix = layout.mtx)
ggsave(paste("./plots/precip-", method, "-post-betatime.pdf", sep = ""),
       panel, width = 9, height = 9)

panel <- arrangeGrob(p3.cur, p3.fut, p4.cur, p4.fut, ncol = 2,
                     layout_matrix = layout.mtx)
ggsave(paste("./plots/precip-", method, "-post-betatimepos.pdf", sep = ""),
       panel, width = 9, height = 9)


#### repeat for GSK ####

method <- "gsk"

load(paste("./cv-results/", method, "-current-25-all.RData", sep = ""))
beta.mu.int  <- apply(fit$beta.int[, , 1], 2, mean)
beta.mu.time <- apply(fit$beta.time[, , 1], 2, mean)
beta.ls.int  <- apply(fit$beta.int[, , 2], 2, mean)
beta.ls.time <- apply(fit$beta.time[, , 2], 2, mean)
prob.beta.mu.time.pos <- apply(fit$beta.time[, , 1] > 0, 2, mean)
prob.beta.ls.time.pos <- apply(fit$beta.time[, , 2] > 0, 2, mean)

main <- bquote(paste(.(toupper(method)), ": Current Posterior Mean of ",
                     beta["1, int"], sep = ""))
fill.legend <- bquote(beta["1, int"])
p <- map.heatmap(lat = s[, 2], lon = s[, 1], data = beta.mu.int,
                 mainTitle = main, legendTitle = fill.legend, color_mid = mid)
ggsave(paste("./plots/precip-", method, "-cur-post-beta1int.pdf", sep = ""),
       p, width = 8, height = 8)

main <- bquote(paste(.(toupper(method)), ": Current Posterior Mean of ",
                     beta["1, time"], sep = ""))
fill.legend <- bquote(beta["1, time"])
p1.cur <- map.heatmap(lat = s[, 2], lon = s[, 1], data = beta.mu.time,
                      mainTitle = main, legendTitle = fill.legend,
                      color_mid = mid, midpoint = 0)
# ggsave(paste("./plots/precip-", method, "-cur-post-beta1time.pdf", sep = ""),
#        p, width = 8, height = 8)

main <- bquote(paste(.(toupper(method)), ": Current Posterior Mean of ",
                     beta["2, time"], sep = ""))
fill.legend <- bquote(beta["2, time"])
p2.cur <- map.heatmap(lat = s[, 2], lon = s[, 1], data = beta.ls.time,
                      mainTitle = main, legendTitle = fill.legend,
                      color_mid = mid, midpoint = 0)
# ggsave(paste("./plots/precip-", method, "-cur-post-beta2time.pdf", sep = ""),
#        p, width = 8, height = 8)

main <- bquote(paste(.(toupper(method)), ": Current Posterior P(", beta["1, time"],
                     " > 0)", sep = ""))
fill.legend <- bquote(paste("P(", beta["1, time"], " > 0)", sep = ""))
p3.cur <- map.heatmap(lat = s[, 2], lon = s[, 1], data = prob.beta.mu.time.pos,
                      mainTitle = main, legendTitle = fill.legend,
                      color_mid = mid, midpoint = 0.5)
# ggsave(paste("./plots/precip-", method, "-cur-post-beta1timepos.pdf", sep = ""),
#        p, width = 8, height = 8)

main <- bquote(paste(.(toupper(method)), ": Current Posterior P(", beta["2, time"],
                     " > 0)", sep = ""))
fill.legend <- bquote(paste("P(", beta["2, time"], " > 0)", sep = ""))
p4.cur <- map.heatmap(lat = s[, 2], lon = s[, 1], data = prob.beta.ls.time.pos,
                      mainTitle = main, legendTitle = fill.legend,
                      color_mid = mid, midpoint = 0.5)
# ggsave(paste("./plots/precip-", method, "-cur-post-beta2timepos.pdf", sep = ""),
#        p, width = 8, height = 8)

#### Get change in q.90 by county
mu.1  <- beta.mu.int + beta.mu.time * time.1
# mu.10 <- beta.mu.int + beta.mu.time * time.10
mu.t  <- beta.mu.int + beta.mu.time * time.t
ls.1  <- beta.ls.int + beta.ls.time * time.1
# ls.10 <- beta.ls.int + beta.ls.time * time.10
ls.t  <- beta.ls.int + beta.ls.time * time.t
xi    <- mean(fit$xi)
q.90.1 <- q.90.t <- rep(NA, ns)
for (i in 1:ns) {
  q.90.1[i] <- qgev(p = 0.90, loc = mu.1[i], scale = exp(ls.1[i]), shape = xi)
  # q.90.10[i] <- qgev(p = 0.90, loc = mu.10[i], scale = exp(ls.10[i]),
  #                    shape = xi)
  q.90.t[i] <- qgev(p = 0.90, loc = mu.t[i], scale = exp(ls.t[i]), shape = xi)
}
q.90.diff.nt <- q.90.t - q.90.1   # change over whole dataset
# q.90.diff.10 <- q.90.t - q.90.10  # change in last 10 years

main <- paste(toupper(method), ": Current Change in q(0.90)", sep = "")
fill.legend <- "Difference"
pdiff.gsk.cur <- map.heatmap(lat = s[, 2], lon = s[, 1], data = q.90.diff.nt,
                             mainTitle = main, legendTitle = fill.legend,
                             color_mid = mid, midpoint = 0)
# ggsave(paste("./plots/precip-", method, "-cur-q90diff.pdf", sep = ""),
#        p, width = 8, height = 8)

load(paste("./cv-results/", method, "-future-25-all.RData", sep = ""))
beta.mu.int  <- apply(fit$beta.int[, , 1], 2, mean)
beta.mu.time <- apply(fit$beta.time[, , 1], 2, mean)
beta.ls.int  <- apply(fit$beta.int[, , 2], 2, mean)
beta.ls.time <- apply(fit$beta.time[, , 2], 2, mean)
prob.beta.mu.time.pos <- apply(fit$beta.time[, , 1] > 0, 2, mean)
prob.beta.ls.time.pos <- apply(fit$beta.time[, , 2] > 0, 2, mean)

main <- bquote(paste(.(toupper(method)), ": Future Posterior Mean of ",
                     beta["1, int"], sep = ""))
fill.legend <- bquote(beta["1, int"])
p <- map.heatmap(lat = s[, 2], lon = s[, 1], data = beta.mu.int,
                 mainTitle = main, legendTitle = fill.legend, color_mid = mid)
# ggsave(paste("./plots/precip-", method, "-fut-post-beta1int.pdf", sep = ""),
#        p, width = 8, height = 8)

main <- bquote(paste(.(toupper(method)), ": Future Posterior Mean of ",
                     beta["1, time"], sep = ""))
fill.legend <- bquote(beta["1, time"])
p1.fut <- map.heatmap(lat = s[, 2], lon = s[, 1], data = beta.mu.time,
                      mainTitle = main, legendTitle = fill.legend,
                      color_mid = mid, midpoint = 0)
# ggsave(paste("./plots/precip-", method, "-fut-post-beta1time.pdf", sep = ""),
#        p, width = 8, height = 8)

main <- bquote(paste(.(toupper(method)), ": Future Posterior Mean of ",
                     beta["2, time"], sep = ""))
fill.legend <- bquote(beta["2, time"])
p2.fut <- map.heatmap(lat = s[, 2], lon = s[, 1], data = beta.ls.time,
                      mainTitle = main, legendTitle = fill.legend,
                      color_mid = mid, midpoint = 0)
# ggsave(paste("./plots/precip-", method, "-fut-post-beta2time.pdf", sep = ""),
#        p, width = 8, height = 8)

main <- bquote(paste(.(toupper(method)), ": Future Posterior P(", beta["1, time"],
                     " > 0)", sep = ""))
fill.legend <- bquote(paste("P(", beta["1, time"], " > 0)", sep = ""))
p3.fut <- map.heatmap(lat = s[, 2], lon = s[, 1], data = prob.beta.mu.time.pos,
                      mainTitle = main, legendTitle = fill.legend,
                      color_mid = mid, midpoint = 0.5)
# ggsave(paste("./plots/precip-", method, "-post-beta1timepos.pdf", sep = ""),
#        p, width = 8, height = 8)

main <- bquote(paste(.(toupper(method)), ": Future Posterior P(", beta["2, time"],
                     " > 0)", sep = ""))
fill.legend <- bquote(paste("P(", beta["2, time"], " > 0)", sep = ""))
p4.fut <- map.heatmap(lat = s[, 2], lon = s[, 1], data = prob.beta.ls.time.pos,
                      mainTitle = main, legendTitle = fill.legend,
                      color_mid = mid, midpoint = 0.5)
# ggsave(paste("./plots/precip-", method, "-fut-post-beta2timepos.pdf", sep = ""),
#        p, width = 8, height = 8)

#### Get change in q.90 by county
mu.1  <- beta.mu.int + beta.mu.time * time.1
# mu.10 <- beta.mu.int + beta.mu.time * time.10
mu.t  <- beta.mu.int + beta.mu.time * time.t
ls.1  <- beta.ls.int + beta.ls.time * time.1
# ls.10 <- beta.ls.int + beta.ls.time * time.10
ls.t  <- beta.ls.int + beta.ls.time * time.t
xi    <- mean(fit$xi)
q.90.1 <- q.90.t <- rep(NA, ns)
for (i in 1:ns) {
  q.90.1[i] <- qgev(p = 0.90, loc = mu.1[i], scale = exp(ls.1[i]), shape = xi)
  # q.90.10[i] <- qgev(p = 0.90, loc = mu.10[i], scale = exp(ls.10[i]),
  #                    shape = xi)
  q.90.t[i] <- qgev(p = 0.90, loc = mu.t[i], scale = exp(ls.t[i]), shape = xi)
}
q.90.diff.nt <- q.90.t - q.90.1   # change over whole dataset
# q.90.diff.10 <- q.90.t - q.90.10  # change in last 10 years

main <- paste(toupper(method), ": Future change in q(0.90)", sep = "")
fill.legend <- "Difference"
pdiff.gsk.fut <- map.heatmap(lat = s[, 2], lon = s[, 1], data = q.90.diff.nt,
                             mainTitle = main, legendTitle = fill.legend,
                             color_mid = mid, midpoint = 0)
# ggsave(paste("./plots/precip-", method, "-fut-q90diff.pdf", sep = ""),
#        p, width = 8, height = 8)

layout.mtx <- matrix(1:4, nrow = 2, ncol = 2)
panel <- arrangeGrob(p1.cur, p1.fut, p2.cur, p2.fut, ncol = 2,
                     layout_matrix = layout.mtx)
ggsave(paste("./plots/precip-", method, "-post-betatime.pdf", sep = ""),
       panel, width = 9, height = 9)

panel <- arrangeGrob(p3.cur, p3.fut, p4.cur, p4.fut, ncol = 2,
                     layout_matrix = layout.mtx)
ggsave(paste("./plots/precip-", method, "-post-betatimepos.pdf", sep = ""),
       panel, width = 9, height = 9)


panel <- arrangeGrob(pdiff.ebf.cur, pdiff.ebf.fut, pdiff.gsk.cur, pdiff.gsk.fut,
                     ncol = 2, layout_matrix = layout.mtx)
ggsave(paste("./plots/precip-q90diff-compare.pdf", sep = ""),
       panel, width = 9, height = 9)

#### Plot set 2: ####
## 1. Delta.mu = mu_2070 - mu_2000
## 2. Delta.ls = ls_2070 - ls_2000
## 3. Delta.q90 = q90_2070 - q90_2000
## 4. P(Delta.mu > 0)
## 5. P(Delta.ls > 0)
## 6. P(Delta.q90 > 0)

method <- "ebf"
this.fit.cur <- get(paste("fit.", method, ".cur", sep = ""))
this.fit.fut <- get(paste("fit.", method, ".fut", sep = ""))
beta.mu.int.cur  <- apply(this.fit.cur$beta.int[, , 1], 2, mean)
beta.mu.int.fut  <- apply(this.fit.fut$beta.int[, , 1], 2, mean)
beta.mu.time.cur <- apply(this.fit.cur$beta.time[, , 1], 2, mean)
beta.mu.time.fut <- apply(this.fit.fut$beta.time[, , 1], 2, mean)

mu.cur <- beta.mu.int.cur + beta.mu.time.cur * time.t
mu.fut <- beta.mu.int.fut + beta.mu.time.fut * time.t
Delta.mu <- mu.fut - mu.cur

beta.ls.int.cur  <- apply(this.fit.cur$beta.int[, , 2], 2, mean)
beta.ls.int.fut  <- apply(this.fit.fut$beta.int[, , 2], 2, mean)
beta.ls.time.cur <- apply(this.fit.cur$beta.time[, , 2], 2, mean)
beta.ls.time.fut <- apply(this.fit.fut$beta.time[, , 2], 2, mean)

ls.cur <- beta.ls.int.cur + beta.ls.time.cur * time.t
ls.fut <- beta.ls.int.fut + beta.ls.time.fut * time.t
Delta.ls <- ls.fut - ls.cur

xi.cur <- mean(this.fit.cur$xi)
xi.fut <- mean(this.fit.fut$xi)
# q90.cur <- q90.fut <- rep(NA, ns)
# for (i in 1:ns) {
#   q90.cur[i] <- qgev(p = 0.90, loc = mu.cur[i], scale = exp(ls.cur[i]),
#                      shape = xi.cur)
#   q90.fut[i] <- qgev(p = 0.90, loc = mu.fut[i], scale = exp(ls.fut[i]),
#                      shape = xi.fut)
# }
# Delta.q90 <- q90.fut - q90.cur

# main <- bquote(paste(.(toupper(method)), ": ", Delta, mu, sep = ""))
main <- bquote(paste(Delta, hat(mu), sep = ""))
fill.legend <- "Difference"
plot.Delta.mu <- map.heatmap(lat = s[, 2], lon = s[, 1], data = Delta.mu,
                             mainTitle = main, legendTitle = fill.legend,
                             color_mid = mid, midpoint = 0)
ggsave(paste("./plots/precip-", method, "-delta-mu.pdf", sep = ""),
       plot.Delta.mu, width = 4.5, height = 4.5)

# main <- bquote(paste(.(toupper(method)), ": ", Delta, "log(", sigma, ")",
#                      sep = ""))
main <- bquote(paste( Delta, "log(", hat(sigma), ")", sep = ""))
fill.legend <- "Difference"
plot.Delta.ls <- map.heatmap(lat = s[, 2], lon = s[, 1], data = Delta.ls,
                             mainTitle = main, legendTitle = fill.legend,
                             color_mid = mid, midpoint = 0)
ggsave(paste("./plots/precip-", method, "-delta-ls.pdf", sep = ""),
       plot.Delta.ls, width = 4.5, height = 4.5)

## P(Delta.mu > 0)
## Get posterior of mu.cur and mu.fut
post.mu.diff <- post.ls.diff <- post.q90.diff <- matrix(NA, niters, ns)
this.q90.cur <- this.q90.fut <- rep(0, ns)

for (i in 1:niters) {
  this.mu.cur <- this.fit.cur$beta.int[i, , 1] + this.fit.cur$beta.time[i, , 1] * time.t
  this.mu.fut <- this.fit.fut$beta.int[i, , 1] + this.fit.fut$beta.time[i, , 1] * time.t
  this.ls.cur <- this.fit.cur$beta.int[i, , 2] + this.fit.cur$beta.time[i, , 2] * time.t
  this.ls.fut <- this.fit.fut$beta.int[i, , 2] + this.fit.fut$beta.time[i, , 2] * time.t
  for (j in 1:ns) {
    this.q90.cur[j] <- qgev(p = 0.90, loc = this.mu.cur[j],
                            scale = exp(this.ls.cur[j]), shape = xi.cur)
    this.q90.fut[j] <- qgev(p = 0.90, loc = this.mu.fut[j],
                            scale = exp(this.ls.fut[j]), shape = xi.fut)
  }


  post.mu.diff[i, ] <- this.mu.fut - this.mu.cur
  post.ls.diff[i, ] <- this.ls.fut - this.ls.cur
  post.q90.diff[i, ] <- this.q90.fut - this.q90.cur
  if (i %% 500 == 0) {
    print(paste("Finished iter ", i, sep = ""))
  }
}

Delta.q90 <- apply(post.q90.diff, 2, mean)
pdiff.mu.pos <- apply(post.mu.diff > 0, 2, mean)
pdiff.ls.pos <- apply(post.ls.diff > 0, 2, mean)
pdiff.q90.pos <- apply(post.q90.diff > 0, 2, mean)

# main <- bquote(paste(.(toupper(method)), ": P[", Delta, mu, " > 0]", sep = ""))
main <- bquote(paste("P[", Delta, hat(mu), " > 0]", sep = ""))
fill.legend <- bquote(paste("P[", Delta, hat(mu), " > 0]", sep = ""))
plot.pdiff.mu <- map.heatmap(lat = s[, 2], lon = s[, 1], data = pdiff.mu.pos,
                             mainTitle = main, legendTitle = fill.legend,
                             color_mid = mid, midpoint = 0.5, zlim = c(0, 1))
ggsave(paste("./plots/precip-", method, "-pdiff-mu-pos.pdf", sep = ""),
       plot.pdiff.mu, width = 6, height = 6)

# main <- bquote(paste(.(toupper(method)), ": P[", Delta, "log(", sigma, ") > 0]",
#                      sep = ""))
main <- bquote(paste("P[", Delta, "log(", hat(sigma), ") > 0]", sep = ""))
fill.legend <- bquote(paste("P[", Delta, "log(", hat(sigma), ") > 0]", sep = ""))
plot.pdiff.ls <- map.heatmap(lat = s[, 2], lon = s[, 1], data = pdiff.ls.pos,
                             mainTitle = main, legendTitle = fill.legend,
                             color_mid = mid, midpoint = 0.5, zlim = c(0, 1))
ggsave(paste("./plots/precip-", method, "-pdiff-ls-pos.pdf", sep = ""),
       plot.pdiff.ls, width = 6, height = 6)

# main <- bquote(paste(.(toupper(method)), ": ", Delta, "Q90", sep = ""))
main <- bquote(paste(Delta, "Q90", sep = ""))
fill.legend <- "Difference"
plot.Delta.q90 <- map.heatmap(lat = s[, 2], lon = s[, 1], data = Delta.q90,
                              mainTitle = main, legendTitle = fill.legend,
                              color_mid = mid, midpoint = 0)
ggsave(paste("./plots/precip-", method, "-delta-q90.pdf", sep = ""),
       plot.Delta.q90, width = 4.5, height = 4.5)

# main <- bquote(paste(.(toupper(method)), ": P[", Delta, "Q90 > 0]", sep = ""))
main <- bquote(paste("P[", Delta, "Q90 > 0]", sep = ""))
fill.legend <- bquote(paste("P[", Delta, "Q90 > 0]", sep = ""))
plot.pdiff.q90 <- map.heatmap(lat = s[, 2], lon = s[, 1], data = pdiff.q90.pos,
                              mainTitle = main, legendTitle = fill.legend,
                              color_mid = mid, midpoint = 0.5, zlim = c(0, 1))
ggsave(paste("./plots/precip-", method, "-pdiff-q90-pos.pdf", sep = ""),
       plot.pdiff.q90, width = 6, height = 6)

layout.mtx <- matrix(1:6, nrow = 3, ncol = 2)
# panel <- arrangeGrob(plot.Delta.mu, plot.Delta.ls, plot.Delta.q90,
#                      plot.pdiff.mu, plot.pdiff.ls, plot.pdiff.q90,
#                      ncol = 2, layout_matrix = layout.mtx)

plots <- list(plot.Delta.mu, plot.Delta.ls, plot.Delta.q90,
              plot.pdiff.mu, plot.pdiff.ls, plot.pdiff.q90)
g     <- lapply(plots, ggplotGrob)
col1 <- rbind(g[[1]], g[[2]], g[[3]], size = "max")
col2 <- rbind(g[[4]], g[[5]], g[[6]], size = "max")
panel <- cbind(col1, col2, size = "max")
ggsave(paste("./plots/precip-", method, "-postpanel.pdf", sep = ""),
       panel, width = 10, height = 13)



#### Plot for GSK ####
method <- "gsk"
this.fit.cur <- get(paste("fit.", method, ".cur", sep = ""))
this.fit.fut <- get(paste("fit.", method, ".fut", sep = ""))
beta.mu.int.cur  <- apply(this.fit.cur$beta.int[, , 1], 2, mean)
beta.mu.int.fut  <- apply(this.fit.fut$beta.int[, , 1], 2, mean)
beta.mu.time.cur <- apply(this.fit.cur$beta.time[, , 1], 2, mean)
beta.mu.time.fut <- apply(this.fit.fut$beta.time[, , 1], 2, mean)

mu.cur <- beta.mu.int.cur + beta.mu.time.cur * time.t
mu.fut <- beta.mu.int.fut + beta.mu.time.fut * time.t
Delta.mu <- mu.fut - mu.cur

beta.ls.int.cur  <- apply(this.fit.cur$beta.int[, , 2], 2, mean)
beta.ls.int.fut  <- apply(this.fit.fut$beta.int[, , 2], 2, mean)
beta.ls.time.cur <- apply(this.fit.cur$beta.time[, , 2], 2, mean)
beta.ls.time.fut <- apply(this.fit.fut$beta.time[, , 2], 2, mean)

ls.cur <- beta.ls.int.cur + beta.ls.time.cur * time.t
ls.fut <- beta.ls.int.fut + beta.ls.time.fut * time.t
Delta.ls <- ls.fut - ls.cur

xi.cur <- mean(this.fit.cur$xi)
xi.fut <- mean(this.fit.fut$xi)
# q90.cur <- q90.fut <- rep(NA, ns)
# for (i in 1:ns) {
#   q90.cur[i] <- qgev(p = 0.90, loc = mu.cur[i], scale = exp(ls.cur[i]),
#                      shape = xi.cur)
#   q90.fut[i] <- qgev(p = 0.90, loc = mu.fut[i], scale = exp(ls.fut[i]),
#                      shape = xi.fut)
# }
# Delta.q90 <- q90.fut - q90.cur

# main <- bquote(paste(.(toupper(method)), ": ", Delta, mu, sep = ""))
main <- bquote(paste(Delta, hat(mu), sep = ""))
fill.legend <- "Difference"
plot.Delta.mu <- map.heatmap(lat = s[, 2], lon = s[, 1], data = Delta.mu,
                             mainTitle = main, legendTitle = fill.legend,
                             color_mid = mid, midpoint = 0)
ggsave(paste("./plots/precip-", method, "-delta-mu.pdf", sep = ""),
       plot.Delta.mu, width = 4.5, height = 4.5)

# main <- bquote(paste(.(toupper(method)), ": ", Delta, "log(", sigma, ")",
#                      sep = ""))
main <- bquote(paste( Delta, "log(", hat(sigma), ")", sep = ""))
fill.legend <- "Difference"
plot.Delta.ls <- map.heatmap(lat = s[, 2], lon = s[, 1], data = Delta.ls,
                             mainTitle = main, legendTitle = fill.legend,
                             color_mid = mid, midpoint = 0)
ggsave(paste("./plots/precip-", method, "-delta-ls.pdf", sep = ""),
       plot.Delta.ls, width = 4.5, height = 4.5)

## P(Delta.mu > 0)
## Get posterior of mu.cur and mu.fut
post.mu.diff <- post.ls.diff <- post.q90.diff <- matrix(NA, niters, ns)
this.q90.cur <- this.q90.fut <- rep(0, ns)

for (i in 1:niters) {
  this.mu.cur <- this.fit.cur$beta.int[i, , 1] + this.fit.cur$beta.time[i, , 1] * time.t
  this.mu.fut <- this.fit.fut$beta.int[i, , 1] + this.fit.fut$beta.time[i, , 1] * time.t
  this.ls.cur <- this.fit.cur$beta.int[i, , 2] + this.fit.cur$beta.time[i, , 2] * time.t
  this.ls.fut <- this.fit.fut$beta.int[i, , 2] + this.fit.fut$beta.time[i, , 2] * time.t
  for (j in 1:ns) {
    this.q90.cur[j] <- qgev(p = 0.90, loc = this.mu.cur[j],
                            scale = exp(this.ls.cur[j]), shape = xi.cur)
    this.q90.fut[j] <- qgev(p = 0.90, loc = this.mu.fut[j],
                            scale = exp(this.ls.fut[j]), shape = xi.fut)
  }


  post.mu.diff[i, ] <- this.mu.fut - this.mu.cur
  post.ls.diff[i, ] <- this.ls.fut - this.ls.cur
  post.q90.diff[i, ] <- this.q90.fut - this.q90.cur
  if (i %% 500 == 0) {
    print(paste("Finished iter ", i, sep = ""))
  }
}

Delta.q90 <- apply(post.q90.diff, 2, mean)
pdiff.mu.pos <- apply(post.mu.diff > 0, 2, mean)
pdiff.ls.pos <- apply(post.ls.diff > 0, 2, mean)
pdiff.q90.pos <- apply(post.q90.diff > 0, 2, mean)

# main <- bquote(paste(.(toupper(method)), ": P[", Delta, mu, " > 0]", sep = ""))
main <- bquote(paste("P[", Delta, hat(mu), " > 0]", sep = ""))
fill.legend <- bquote(paste("P[", Delta, hat(mu), " > 0]", sep = ""))
plot.pdiff.mu <- map.heatmap(lat = s[, 2], lon = s[, 1], data = pdiff.mu.pos,
                             mainTitle = main, legendTitle = fill.legend,
                             color_mid = mid, midpoint = 0.5, zlim = c(0, 1))
ggsave(paste("./plots/precip-", method, "-pdiff-mu-pos.pdf", sep = ""),
       plot.pdiff.mu, width = 6, height = 6)

# main <- bquote(paste(.(toupper(method)), ": P[", Delta, "log(", sigma, ") > 0]",
#                      sep = ""))
main <- bquote(paste("P[", Delta, "log(", hat(sigma), ") > 0]", sep = ""))
fill.legend <- bquote(paste("P[", Delta, "log(", hat(sigma), ") > 0]", sep = ""))
plot.pdiff.ls <- map.heatmap(lat = s[, 2], lon = s[, 1], data = pdiff.ls.pos,
                             mainTitle = main, legendTitle = fill.legend,
                             color_mid = mid, midpoint = 0.5, zlim = c(0, 1))
ggsave(paste("./plots/precip-", method, "-pdiff-ls-pos.pdf", sep = ""),
       plot.pdiff.ls, width = 6, height = 6)

# main <- bquote(paste(.(toupper(method)), ": ", Delta, "Q90", sep = ""))
main <- bquote(paste(Delta, "Q90", sep = ""))
fill.legend <- "Difference"
plot.Delta.q90 <- map.heatmap(lat = s[, 2], lon = s[, 1], data = Delta.q90,
                              mainTitle = main, legendTitle = fill.legend,
                              color_mid = mid, midpoint = 0)
ggsave(paste("./plots/precip-", method, "-delta-q90.pdf", sep = ""),
       plot.Delta.q90, width = 4.5, height = 4.5)

# main <- bquote(paste(.(toupper(method)), ": P[", Delta, "Q90 > 0]", sep = ""))
main <- bquote(paste("P[", Delta, "Q90 > 0]", sep = ""))
fill.legend <- bquote(paste("P[", Delta, "Q90 > 0]", sep = ""))
plot.pdiff.q90 <- map.heatmap(lat = s[, 2], lon = s[, 1], data = pdiff.q90.pos,
                              mainTitle = main, legendTitle = fill.legend,
                              color_mid = mid, midpoint = 0.5, zlim = c(0, 1))
ggsave(paste("./plots/precip-", method, "-pdiff-q90-pos.pdf", sep = ""),
       plot.pdiff.q90, width = 6, height = 6)

layout.mtx <- matrix(1:6, nrow = 3, ncol = 2)
# panel <- arrangeGrob(plot.Delta.mu, plot.Delta.ls, plot.Delta.q90,
#                      plot.pdiff.mu, plot.pdiff.ls, plot.pdiff.q90,
#                      ncol = 2, layout_matrix = layout.mtx)

plots <- list(plot.Delta.mu, plot.Delta.ls, plot.Delta.q90,
              plot.pdiff.mu, plot.pdiff.ls, plot.pdiff.q90)
g     <- lapply(plots, ggplotGrob)
col1 <- rbind(g[[1]], g[[2]], g[[3]], size = "max")
col2 <- rbind(g[[4]], g[[5]], g[[6]], size = "max")
panel <- cbind(col1, col2, size = "max")
ggsave(paste("./plots/precip-", method, "-postpanel.pdf", sep = ""),
       panel, width = 10, height = 13)