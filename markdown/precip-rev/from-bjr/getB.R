rm(list=ls())
setwd("S:\\Documents\\PCAX\\Revision\\extreme-decomp-master\\markdown\\")
load("precip_rev\\precip_preprocess.RData")
source("Code\\auxfunctions.R")
source("Code\\PCAX_fixed_alpha.R")

ns   <- nrow(Y)
nt   <- ncol(Y)

set.seed(0820)
fold <- matrix(sample(1:5,ns*nt,replace=TRUE),ns,nt)
l    <- 40

print(l)

out  <- list()


for(f in 1:6){
   print(f)
   Ym          <- Y
   Ym[fold==f] <- NA
   EC.hat      <- get.chi(Ym)
   out[[f]]    <- get.factors.EC(EC.hat,L=l,s=s,verbose=TRUE)
}

junk   <- ls()
OUTPUT <- out
L      <- l
save.image(paste0("precip_rev\\basis_L",L,".RData"))
rm(list=junk)
rm(junk)
save.image(paste0("precip_rev\\basis_L",L,".RData"))



