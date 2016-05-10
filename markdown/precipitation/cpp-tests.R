library(microbenchmark)

temp.mtx <- matrix(rnorm(30000), 1000, 30)
temp.vec1 <- rnorm(1000)
temp.vec2 <- rnorm(30)

mean(sweepC2plus(temp.mtx, temp.vec2) / sweepC2plus2(temp.mtx, temp.vec2))
mean(sweepC2plus(temp.mtx, temp.vec2) / sweep(temp.mtx, 2, temp.vec2, "+"))

microbenchmark(sweep(temp.mtx, 2, temp.vec2, "+"),
               sweepC2plus(temp.mtx, temp.vec2))


mean(sweepC1plus(temp.mtx, temp.vec1) / sweep(temp.mtx, 1, temp.vec1, "+"))

microbenchmark(sweep(temp.mtx, 1, temp.vec1, "+"),
               sweepC1plus(temp.mtx, temp.vec1))
