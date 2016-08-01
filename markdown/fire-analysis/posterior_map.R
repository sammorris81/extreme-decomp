rm(list = ls())
source("./package_load.R", chdir = TRUE)
library(gridExtra)
library(evd)
# trying to see if there's a difference in the estimates for gsk vs ebf

#### load in the data ####
load(file = "../../code/analysis/fire/georgia_preprocess/fire_data.RData")

#### spatial setup ####
# get Georgia coordinates from georgia_preprocess in code/analysis/fire
load(file = "../../code/analysis/fire/georgia_preprocess/georgia_map.RData")
d <- rdist(cents)
diag(d) <- 0
ns <- nrow(cents)
nt <- nrow(Y)  # Y is nt x ns
time.1  <- (1 - nt / 2) / nt
time.10 <- ((nt - 10) - nt / 2) / nt
time.t <- (nt - nt / 2) / nt

# cleanup county names for mapping
county <- tolower(county)
county[9] <- "ben hill"
county[15] <- "bryan"
county[44] <- "de kalb"
county[60] <- "fulton"
county[67] <- "gwinnett"
county[80] <- "jeff davis"
# mid <- "grey90"
mid <- "#FFFFFF"

#### County by county plots of: ####
#  1. Posterior means of betas: 4
#  2. Posterior P(B_time > 0): 2
#  3. Change in q(90)

#### EBF spatial kernel functions
method <- "ebf"
load(paste("./cv-results/", method, "-gsk-25-all.RData", sep = ""))
beta.mu.int  <- apply(fit$beta.int[, , 1], 2, mean)
beta.mu.time <- apply(fit$beta.time[, , 1], 2, mean)
beta.ls.int  <- apply(fit$beta.int[, , 2], 2, mean)
beta.ls.time <- apply(fit$beta.time[, , 2], 2, mean)
prob.beta.mu.time.pos <- apply(fit$beta.time[, , 1] > 0, 2, mean)
prob.beta.ls.time.pos <- apply(fit$beta.time[, , 2] > 0, 2, mean)

main <- bquote(paste(.(toupper(method)), ": Posterior Mean of ", beta["1, int"],
                     sep = ""))
fill.legend <- bquote(beta["1, int"])
p <- map.ga.ggplot(Y = beta.mu.int, counties = county,
                   midpoint = mean(beta.mu.int),
                   mid = mid, fill.legend = fill.legend, main = main)
# ggsave(paste("./plots/fire-", method, "-post-beta1int.pdf", sep = ""),
#        p, width = 4.5, height = 4.5)

main <- bquote(paste(.(toupper(method)), ": Posterior Mean of ", beta["2, int"],
                     sep = ""))
fill.legend <- bquote(paste(beta["2, int"]))
p <- map.ga.ggplot(Y = beta.ls.int, counties = county,
                   midpoint = mean(beta.ls.int),
                   mid = mid, fill.legend = fill.legend, main = main)
# ggsave(paste("./plots/fire-", method, "-post-beta2int.pdf", sep = ""),
#        p, width = 4.5, height = 4.5)

main <- bquote(paste(.(toupper(method)), ":Posterior Mean of ", beta["1, time"],
                     sep = ""))
fill.legend <- bquote(beta["1, time"])
p1 <- map.ga.ggplot(Y = beta.mu.time, counties = county,
                   midpoint = 0,
                   mid = mid, fill.legend = fill.legend, main = main)

main <- bquote(paste(.(toupper(method)), ": Posterior Mean of ", beta["2, time"],
                     sep = ""))
fill.legend <- bquote(paste(beta["2, time"]))
p2 <- map.ga.ggplot(Y = beta.ls.time, counties = county,
                   midpoint = 0,
                   mid = mid, fill.legend = fill.legend, main = main)

layout.mtx <- matrix(1:2, nrow = 1, ncol = 2)
panel <- arrangeGrob(p1, p2, ncol = 2, layout_matrix = layout.mtx)

ggsave(paste("./plots/fire-", method, "-post-betatime.pdf", sep = ""),
       panel, width = 9, height = 4.5)


main <- bquote(paste(.(toupper(method)), ": Posterior P(", beta["1, time"],
                       " > 0)", sep = ""))
fill.legend <- bquote(paste("P(", beta["1, time"], " > 0)", sep = ""))
p1 <- map.ga.ggplot(Y = prob.beta.mu.time.pos, counties = county, midpoint = 0.5,
                   mid = mid, fill.legend = fill.legend, main = main)

main <- bquote(paste(.(toupper(method)), ": Posterior P(", beta["2, time"],
                     " > 0)", sep = ""))
fill.legend <- bquote(paste("P(", beta["2, time"], " > 0)", sep = ""))
p2 <- map.ga.ggplot(Y = prob.beta.ls.time.pos, counties = county, midpoint = 0.5,
              mid = mid, fill.legend = fill.legend, main = main)
ggsave(paste("./plots/fire-", method, "-post-betatimepos.pdf", sep = ""),
       panel, width = 9, height = 4.5)

#### Get change in q.90 by county
mu.1  <- beta.mu.int + beta.mu.time * time.1
mu.10 <- beta.mu.int + beta.mu.time * time.10
mu.t  <- beta.mu.int + beta.mu.time * time.t
ls.1  <- beta.ls.int + beta.ls.time * time.1
ls.10 <- beta.ls.int + beta.ls.time * time.10
ls.t  <- beta.ls.int + beta.ls.time * time.t
xi    <- mean(fit$xi)
q.90.1 <- q.90.10 <- q.90.t <- rep(NA, ns)
for (i in 1:ns) {
  q.90.1[i] <- qgev(p = 0.90, loc = mu.1[i], scale = exp(ls.1[i]), shape = xi)
  q.90.10[i] <- qgev(p = 0.90, loc = mu.10[i], scale = exp(ls.10[i]),
                     shape = xi)
  q.90.t[i] <- qgev(p = 0.90, loc = mu.t[i], scale = exp(ls.t[i]), shape = xi)
}
q.90.diff.nt <- q.90.t - q.90.1   # change over whole dataset
q.90.diff.10 <- q.90.t - q.90.10  # change in last 10 years

# there are three exceptionally high changes
idx.1 <- which.max(q.90.diff.nt)
idx.2 <- which(q.90.diff.nt == max(q.90.diff.nt[-idx.1]))
idx.3 <- which(q.90.diff.nt == max(q.90.diff.nt[-c(idx.1, idx.2)]))
these <- c(idx.1, idx.2, idx.3)
cents.text <- q.90.diff.nt[these]
round(cents.text)  # c(6084, 6693, 9690) - Charlton, Clinch, Ware
cents.plot <- cents[these, ]
cents.plot[1, 1] <- -82.5
cents.plot[1, 2] <- 31.2
cents.plot[2, 2] <- 30.95
cents.plot[3, 1] <- -82.2
cents.plot[3, 2] <- 30.9
q.90.diff.nt[these] <- 2500

main <- paste(toupper(method), ": Change in q(0.90)", sep = "")
fill.legend <- "Difference"
pdiff.ebf <- map.ga.ggplot(Y = q.90.diff.nt, counties = county, midpoint = 0,
                           mid = mid, fill.legend = fill.legend, main = main,
                           cents = cents.plot, cents.text = c(1, 2, 3))


#### GSK spatial kernel functions
method <- "gsk"
load(paste("./cv-results/", method, "-gsk-25-all.RData", sep = ""))
beta.mu.int  <- apply(fit$beta.int[, , 1], 2, mean)
beta.mu.time <- apply(fit$beta.time[, , 1], 2, mean)
beta.ls.int  <- apply(fit$beta.int[, , 2], 2, mean)
beta.ls.time <- apply(fit$beta.time[, , 2], 2, mean)
prob.beta.mu.time.pos <- apply(fit$beta.time[, , 1] > 0, 2, mean)
prob.beta.ls.time.pos <- apply(fit$beta.time[, , 2] > 0, 2, mean)

main <- bquote(paste(.(toupper(method)), ": Posterior Mean of ", beta["1, int"],
                     sep = ""))
fill.legend <- bquote(beta["1, int"])
p <- map.ga.ggplot(Y = beta.mu.int, counties = county,
                   midpoint = mean(beta.mu.int),
                   mid = mid, fill.legend = fill.legend, main = main)
# ggsave(paste("./plots/fire-", method, "-post-beta1int.pdf", sep = ""),
#        p, width = 4.5, height = 4.5)

main <- bquote(paste(.(toupper(method)), ": Posterior Mean of ", beta["2, int"],
                     sep = ""))
fill.legend <- bquote(paste(beta["2, int"]))
p <- map.ga.ggplot(Y = beta.ls.int, counties = county,
                   midpoint = mean(beta.ls.int),
                   mid = mid, fill.legend = fill.legend, main = main)
# ggsave(paste("./plots/fire-", method, "-post-beta2int.pdf", sep = ""),
#        p, width = 4.5, height = 4.5)

main <- bquote(paste(.(toupper(method)), ":Posterior Mean of ", beta["1, time"],
                     sep = ""))
fill.legend <- bquote(beta["1, time"])
p1 <- map.ga.ggplot(Y = beta.mu.time, counties = county,
                   midpoint = 0,
                   mid = mid, fill.legend = fill.legend, main = main)

main <- bquote(paste(.(toupper(method)), ": Posterior Mean of ", beta["2, time"],
                     sep = ""))
fill.legend <- bquote(paste(beta["2, time"]))
p2 <- map.ga.ggplot(Y = beta.ls.time, counties = county,
                   midpoint = 0,
                   mid = mid, fill.legend = fill.legend, main = main)

layout.mtx <- matrix(1:2, nrow = 1, ncol = 2)
panel <- arrangeGrob(p1, p2, ncol = 2, layout_matrix = layout.mtx)

ggsave(paste("./plots/fire-", method, "-post-betatime.pdf", sep = ""),
       panel, width = 9, height = 4.5)

main <- bquote(paste(.(toupper(method)), ": Posterior P(", beta["1, time"],
                     " > 0)", sep = ""))
fill.legend <- bquote(paste("P(", beta["1, time"], " > 0)", sep = ""))
p1 <- map.ga.ggplot(Y = prob.beta.mu.time.pos, counties = county, midpoint = 0.5,
                   mid = mid, fill.legend = fill.legend, main = main)


main <- bquote(paste(.(toupper(method)), ": Posterior P(", beta["2, time"],
                     " > 0)", sep = ""))
fill.legend <- bquote(paste("P(", beta["2, time"], " > 0)", sep = ""))
p2 <- map.ga.ggplot(Y = prob.beta.ls.time.pos, counties = county, midpoint = 0.5,
                   mid = mid, fill.legend = fill.legend, main = main)

layout.mtx <- matrix(1:2, nrow = 1, ncol = 2)
panel <- arrangeGrob(p1, p2, ncol = 2, layout_matrix = layout.mtx)
ggsave(paste("./plots/fire-", method, "-post-betatimepos.pdf", sep = ""),
       panel, width = 9, height = 4.5)

#### Get change in q.90 by county
mu.1  <- beta.mu.int + beta.mu.time * time.1
mu.10 <- beta.mu.int + beta.mu.time * time.10
mu.t  <- beta.mu.int + beta.mu.time * time.t
ls.1  <- beta.ls.int + beta.ls.time * time.1
ls.10 <- beta.ls.int + beta.ls.time * time.10
ls.t  <- beta.ls.int + beta.ls.time * time.t
xi    <- mean(fit$xi)
q.90.1 <- q.90.10 <- q.90.t <- rep(NA, ns)
for (i in 1:ns) {
  q.90.1[i] <- qgev(p = 0.90, loc = mu.1[i], scale = exp(ls.1[i]), shape = xi)
  q.90.10[i] <- qgev(p = 0.90, loc = mu.10[i], scale = exp(ls.10[i]),
                     shape = xi)
  q.90.t[i] <- qgev(p = 0.90, loc = mu.t[i], scale = exp(ls.t[i]), shape = xi)
}
q.90.diff.nt <- q.90.t - q.90.1   # change over whole dataset
q.90.diff.10 <- q.90.t - q.90.10  # change in last 10 years

# there are three exceptionally high changes
idx.1 <- which.max(q.90.diff.nt)
idx.2 <- which(q.90.diff.nt == max(q.90.diff.nt[-idx.1]))
idx.3 <- which(q.90.diff.nt == max(q.90.diff.nt[-c(idx.1, idx.2)]))
these <- c(idx.1, idx.2, idx.3)
cents.text <- q.90.diff.nt[these]
round(cents.text)  # c(6520, 6794, 8934) - Clinch, Charlton, Ware
cents.plot <- cents[these, ]
cents.plot[1, 1] <- -82.5
cents.plot[1, 2] <- 31.2
cents.plot[2, 1] <- -82.2
cents.plot[2, 2] <- 30.9
cents.plot[3, 2] <- 30.95
q.90.diff.nt[these] <- 2500

main <- paste(toupper(method), ": Change in q(0.90)", sep = "")
fill.legend <- "Difference"
pdiff.gsk <- map.ga.ggplot(Y = q.90.diff.nt, counties = county, midpoint = 0,
                           mid = mid, fill.legend = fill.legend, main = main,
                           cents = cents.plot, cents.text = c(1, 3, 2))
layout.mtx <- matrix(1:2, nrow = 1, ncol = 2)
panel <- arrangeGrob(pdiff.ebf, pdiff.gsk, ncol = 2, layout_matrix = layout.mtx)

ggsave(paste("./plots/fire-q90diff-compare.pdf", sep = ""),
       panel, width = 9, height = 4.5)