library(microbenchmark)

temp.mtx <- matrix(rnorm(30000), 1000, 30)
temp.vec1 <- rnorm(1000)
temp.vec2 <- rnorm(30)

mean(sweepC1plus(temp.mtx, temp.vec1) / sweep(temp.mtx, 1, temp.vec1, "+"))

microbenchmark(sweep(temp.mtx, 1, temp.vec1, "+"),
               sweepC1plus(temp.mtx, temp.vec1))

mean(sweepC2plus(temp.mtx, temp.vec2) / sweep(temp.mtx, 2, temp.vec2, "+"))

microbenchmark(sweep(temp.mtx, 2, temp.vec2, "+"),
               sweepC2plus(temp.mtx, temp.vec2))

mean(sweepC1times(temp.mtx, temp.vec1) / sweep(temp.mtx, 1, temp.vec1, "*"))

microbenchmark(sweep(temp.mtx, 1, temp.vec1, "*"),
               sweepC1times(temp.mtx, temp.vec1))


mean(sweepC2times(temp.mtx, temp.vec2) / sweep(temp.mtx, 2, temp.vec2, "*"))

microbenchmark(sweep(temp.mtx, 2, temp.vec2, "*"),
               sweepC2times(temp.mtx, temp.vec2))


microbenchmark(temp.mtx^0.4, powC(temp.mtx, 0.4))
