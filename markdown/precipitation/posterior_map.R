rm(list = ls())
source("./package_load.R", chdir = TRUE)
library(evd)
# trying to see if there's a difference in the estimates for gsk vs ebf

load("precip_preprocess.RData")
nx <- length(unique(s[, 1]))
ny <- length(unique(s[, 2]))

ns <- nrow(Y)
nt <- ncol(Y) / 2

time.1 <- (1 - nt / 2) / nt
time.t <- (nt - nt / 2) / nt
# mid <- "grey90"
mid <- "#FFFFFF"

#### Spatial plots of: ####
#  1. Posterior means of betas: 4
#  2. Posterior P(B_time > 0): 2
#  3. Change in q(90)

method <- "ebf"

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