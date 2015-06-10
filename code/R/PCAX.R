
######################################################################
# 
# Function to estimate the B functions:
#
# Inputs:
#
#  EC      := n x n matrix of estimated pairwise extremal coefficients
#  L       := number of B functions to be estimated
#  alpha   := positive stable parameter alpha
#  init.B  := inital value
#  iters   := number of iterations in the optimization algorithm
#
# Outputs
#
#  est       := estimated value of B
#  alpha     := estimated alpha
#  EC.smooth := smoothed version of EC
# 
######################################################################

get.factors.EC <- function(EC,L=5,s=NULL,bw=NULL,alpha=NULL,init.B=NULL,iters=10,verbose=TRUE){

 tick   <- proc.time()[3]

 n      <- ncol(EC)
 if( is.null(s)){  s <- 1:n}
 if(is.null(bw)){ bw <- 2*min(dist(s))}

 # SMOOTHING

  EC  <- Ksmooth(EC,s,bw)
  ECs <- EC
  if(is.null(alpha)){alpha  <- log2(mean(diag(EC)))}
  diag(EC) <- NA

 # INITIAL VALUES

  if(is.null(init.B)){
    B <- 1/matrix(1:L,n,L,byrow=TRUE)
  }
  if(!is.null(init.B)){
    B <- init.B
  }
  B  <- sweep(B,1,rowSums(B),"/")
 
 # ESTIMATION

  Delta_B   <- rep(NA,iters)
  Delta_val <- rep(NA,iters)

  for(iter in 1:iters){
   prev  <- B
   maxit <- ifelse(iter==iters | iter>100,100,2*iter+2)
   for(i in 1:n){
    fit       <- optim(B[i,],fn=SSE,gr=SSE.grad,Y=EC[i,],B2=B,alpha=alpha,
                       lower=rep(0,L),upper=rep(1,L),method="L-BFGS-B",
                       control=list(maxit=maxit))
    B[i,] <- abs(fit$par)/sum(abs(fit$par))
   }
   Delta_B[iter]   <- mean((prev-B)^2)
   Delta_val[iter] <- sum((EC-make.EC(B,alpha))^2,na.rm=TRUE)

   if(verbose){print(paste("Done with iteration",iter,"of",iters))}
  }

 # REORDER THE COLUMNS
    
  B      <- B[,order(-colSums(B))]
  pct    <- colSums(B)/sum(B)
  tock   <- proc.time()[3]

  output <- list(est=B,pct=pct,alpha=alpha,EC.smooth=ECs,
                 Delta.B=Delta_B,Delta.val=Delta_val,
                 seconds=tock-tick)
 
return(output)}




# SSE for row of Y-EC
SSE <- function(B1,B2,Y,alpha,lambda=1000){

   BB  <- B1^(1/alpha)
   B2  <- B2^(1/alpha)
   EC  <- sweep(B2,2,BB,"+")
   EC  <- rowSums(EC^alpha)
   sse <- sum((Y-EC)^2,na.rm=TRUE)+lambda*(sum(B1)-1)^2
return(sse)}

SSE.grad <- function(B1,B2,Y,alpha,lambda=1000){

   BB   <- B1^(1/alpha)
   B2   <- B2^(1/alpha)

   BB   <- sweep(B2,2,BB,"+")
   EC0  <- rowSums(BB^alpha)

   EC1  <- BB^(alpha-1)
   EC1  <- sweep(EC1,2,B1^(1/alpha-1),"*")
   EC1  <- sweep(EC1,1,Y-EC0,"*")
 
   grad <- -2*colSums(EC1,na.rm=TRUE)+
            2*lambda*(sum(B1)-1)

return(grad)}


make.EC  <- function(B,alpha){
   Ba    <- B^(1/alpha) 
   EC    <- NULL
   for(j in 1:nrow(B)){
      BB <- sweep(Ba,2,Ba[j,],"+")
      EC <- cbind(EC,rowSums(BB^alpha))
   }
return(EC)}

# Performs kernel smoothing of the extremal coefficient matrix.
Ksmooth <- function(ECmat,s=NULL,bw=NULL){

   n           <- nrow(ECmat)
   diag(ECmat) <- 0
   E1          <- ifelse(ECmat==0,0,1)
   if(is.null(s)){  s<-1:n}
   if(is.null(bw)){bw<-2*min(dist(s))}


   d2       <- as.matrix(dist(s)/bw)^2
   W        <- exp(-d2)
   diag(W)  <- 0

   num      <- W%*%ECmat%*%W
   den      <- W%*%E1%*%W

   ECsmooth <- num/den

return(ECsmooth)}


####################################################
# SIMPLE EXAMPLE
####################################################

if(TRUE){

 library(splines)
 library(fields)

 n     <- 100
 L     <- 5
 alpha <- 0.3

 #set.seed(0820)

 #Define the truth
 
  B.true  <- bs(1:n,df=L,intercept=TRUE)
  B.true  <- sweep(B.true,1,rowSums(B.true),"/")
  tot     <- colSums(B.true)
  B.true  <- B.true[,order(-tot)]
  EC.true <- make.EC(B.true,alpha)

 # The estimate (not generated in a realistic way)

  junk    <- matrix(rnorm(n^2),n,n)
  EC.hat  <- EC.true+0.01*t(junk)%*%junk
  EC.hat  <- ifelse(EC.hat<1,1,EC.hat)
  EC.hat  <- ifelse(EC.hat>2,2,EC.hat)

  diag(EC.hat) <- NA

 # Estimation

  
  out       <- get.factors.EC(EC.hat,L=L,s=1:n,bw=5)
  B.est     <- out$est
  alphahat  <- out$alpha
  EC.smooth <- out$EC.smooth
  EC.est    <- make.EC(B.est,alphahat)

  print(out$pct)

 # Plot the results

  par(mfrow=c(3,2))
  matplot(B.true,type="l",main="True B")
  matplot(B.est,type="l",main="Estimated B")
  image.plot(1:n,1:n,EC.true,main="True EC")
  image.plot(1:n,1:n,EC.hat,main="Initial EC estimate (theta-hat)")
  image.plot(1:n,1:n,EC.smooth,main="Smoothed EC (theta-tilde)")
  image.plot(1:n,1:n,EC.est,main="Final EC estimate")

 
}
