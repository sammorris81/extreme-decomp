cty <- read.csv("county.csv", as.is = TRUE)

i1 <- which(cty[, 2] == " GA")

id1 <- substr(cty[i1, 1], 1, 6)

neigh <- NULL
for (j in 1:158) {
  ccc <- cty[(i1[j] + 1):(i1[j + 1] - 1), ]
  ccc <- ccc[ccc[, 4] == " GA", ]
  neigh[[j]] <- substr(ccc[, 3], 1, 6)
}

ADJ <- matrix(0, 159, 159)

for (j in 1:158) {
  nnn <- neigh[[j]]
  for (l in 1:length(nnn)) {
    eee <- which(id1 == nnn[l])
    ADJ[j, eee] <- ADJ[eee, j] <- 1
  }
}

diag(ADJ) <- 0
rm(cty, i1, id1, neigh, j, ccc, l, nnn, eee)