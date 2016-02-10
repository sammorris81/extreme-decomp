rm(list=ls())

source(file = "./package_load.R", chdir = T)

procs  <- c("basis", "kern")  # what process determines the spatial structure
nprocs <- length(procs)
bases  <- c(2, 5, 10, 15, 20)
nbases <- length(bases)
folds  <- 1:5
nfolds <- length(folds)
timing <- matrix(NA, nprocs * nbases, length(folds))
these.rownames <- paste(rep(procs, each = nbases), 
                        rep(bases, times = nprocs))
rownames(timing) <- these.rownames
colnames(timing) <- paste("Fold:", folds)

source.files <- c("./fitmodel_basis_time.R", "./fitmodel_kern_time.R")

for (p in 1:nprocs) {
  for (b in 1:nbases) {
    this.row <- (p - 1) * nbases + b
    for (f in 1:nfolds) {
      L <- bases[b]
      cv <- f
      source(file = source.files[p])
      timing[this.row, f] <- fit$timing
      write.table(timing, file = "./cv-tables/timing.txt")
    }
    cat("Finished", bases[b], "basis functions for", procs[p], "spatial process \n")
  }
}
