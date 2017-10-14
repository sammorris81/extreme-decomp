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
#  iters   := minimum number of iterations in the optimization algorithm.
#             the algorithm will continue until convergence is reached at
#             all sites.
#
# Outputs
#
#  est       := estimated value of B
#  alpha     := estimated alpha
#  EC.smooth := smoothed version of EC
#
######################################################################

udpate_one_row <- function(theta,Bcur,Brest,alpha,maxit=100,fudge=0.001){
  b    <- optim(log(Bcur)-mean(log(Bcur)),fn=SSE1row,gr=SSE1rowgrad,
                theta=theta,Brest=Brest,alpha=alpha,fudge=fudge,
                control=list(maxit=maxit))$par
  B    <- exp(b)/sum(exp(b))
  return(B)}

SSE1row <- function(b,theta,Brest,alpha){
  Bcur     <- exp(b)/sum(exp(b))
  p1       <- sweep(Brest^(1/alpha),2,Bcur^(1/alpha),"+")
  thetahat <- rowSums(p1^alpha)
  sse      <- sum((theta-thetahat)^2) 
  return(sse)}



SSEall <- function(alpha,theta,B){
  #   thetahat <- B2EC(B,alpha)
  thetahat <- B2ECcpp(B^(1/alpha),alpha)
  sse      <- sum((theta-thetahat)^2,na.rm=TRUE)
  return(sse)}

B2EC <- function(B,alpha){
  n     <- nrow(B)
  B     <- B^(1/alpha)
  theta <- 0
  for(l in 1:ncol(B)){
    Bl    <- matrix(B[,l],n,n,byrow=TRUE)+
      matrix(B[,l],n,n,byrow=FALSE)
    theta <- theta + Bl^alpha   
  }
  diag(theta) <- 2^alpha
  return(theta)}


sse <- function(ab,EC,n,L,alpha,fudge=0.001){
  b     <- matrix(ab,n,L)
  B     <- sweep(exp(b),1,rowSums(exp(b)),"/")
  SSEall(alpha,EC,B)/(n*L)+fudge*sum(ab^2)/(n*L)}

sse_grad <- function(ab,EC,n,L,alpha,fudge=0.001){
  b     <- matrix(ab,n,L)
  B     <- sweep(exp(b),1,rowSums(exp(b)),"/")
  
  Gb <- b
  for(j in 1:n){
    Gb[j,] <- SSE1rowgrad(b[j,],EC[j,-j],B[-j,],alpha=alpha)
  }
  
  G  <- as.vector(Gb)+2*fudge*ab
  
  return(G/(n*L))}

sse_grad <- function(ab,EC,n,L,alpha,fudge=0.001){
  b       <- matrix(ab,n,L)
  B       <- sweep(exp(b),1,rowSums(exp(b)),"/")
  Ba      <- B^(1/alpha)
  theta   <- B2ECcpp(Ba,alpha)
  
  R       <- EC-theta  
  diag(R) <- 0   
  Gb      <- 0*b
  for(k in 1:L){
    delta  <- (matrix(Ba[,k],n,n,byrow=TRUE)+
                 matrix(Ba[,k],n,n,byrow=FALSE))^(alpha-1)
    Rd     <- rowSums(R*delta)
    Gb     <- Gb     + sweep(B,1,Rd*Ba[,k],"*")
    Gb[,k] <- Gb[,k] - Rd*Ba[,k]    
  }
  G  <- 2*as.vector(Gb)+2*fudge*ab
  return(G/(n*L))}



SSE1rowgrad <- function(b,theta,Brest,alpha){
  B         <- exp(b)/sum(exp(b))
  delta     <- sweep(Brest^(1/alpha),2,B^(1/alpha),"+")
  thetahat  <- rowSums(delta^alpha)
  G         <- diag(B^(1/alpha-1))%*%(diag(B)-outer(B,B))
  thetahat2 <- delta^(alpha-1)
  q0        <- -2*colSums(sweep(thetahat2,1,theta-thetahat,"*"))
  grad      <- q0%*%G 
  return(grad)}

B2EC <- function(B,alpha){
  n     <- nrow(B)
  B     <- B^(1/alpha)
  theta <- 0
  for(l in 1:ncol(B)){
    Bl    <- matrix(B[,l],n,n,byrow=TRUE)+
      matrix(B[,l],n,n,byrow=FALSE)
    theta <- theta + Bl^alpha   
  }
  diag(theta) <- 2^alpha
  return(theta)}




library(inline)
library(Rcpp)

# Rcpp function to multiply a matrix by a constant

B2EC_code = "

Rcpp::NumericMatrix Bcpp(Ba);
double a = Rcpp::as<double>(alpha);

int n = Bcpp.nrow();
int L = Bcpp.ncol();
Rcpp::NumericMatrix theta(n, n);

for(int l = 0; l<L; l++){
for (int i = 0; i < n; i++){
for (int j = 0; j < i; j++){
theta(i,j) += pow(Bcpp(i,l)+Bcpp(j,l),a);
}
}
}

for (int i = 0; i < n; i++){
for (int j = 0; j < i; j++){
theta(j,i) += theta(i,j);
}
}

for (int i = 0; i < n; i++){
theta(i,i) = pow(2,a);
}

return theta;
"

B2ECcpp <- cxxfunction(signature(Ba = "numeric", alpha = "numeric"), 
                       body = B2EC_code, 
                       plugin = "Rcpp")




get.factors.EC <- function(EChat, L = 5, s = NULL, 
                           bw = NULL, init.B = NULL, fudge=0.001,
                           maxit1 = 500, maxit2=500, eps=0.0001,verbose = TRUE){
  tick   <- proc.time()[3]
  n      <- ncol(EChat)
  if (is.null(bw)) {bw <- 2 * min(dist(s)) }
  
  # SMOOTHING
  
  ECs        <- KsmoothCV(EChat, s)  # run a kernel smoother on the pairwise estimates
  EC         <- ECs
  diag(EC)   <- NA
  
  # Initial value
  alpha      <- 0.5  
  alpha.grid <- alpha
  ab         <- rnorm(n*L)
  opt        <- optim(ab,fn=sse,gr=sse_grad,n=n,L=L,EC=EC,alpha=alpha.grid[1],fudge=0.001,
                      method="BFGS",control = list(trace=ifelse(verbose,50,0),maxit=maxit2,reltol=0.00000001))
  val        <- opt$val
  init.ab    <- opt$par
  ab         <- init.ab
  bestval    <- val
  
  # Refinement
  
  res        <- c(0.2,0.1,0.05,0.025)
  for(round in 1:length(res)){
    alpha.grid <- c(alpha+res[round]*c(-1,1),alpha.grid)
    alpha.grid <- ifelse(alpha.grid>0.9,0.9,alpha.grid)
    alpha.grid <- ifelse(alpha.grid<0.1,0.1,alpha.grid)
    val        <- c(rep(NA,2),val)
    for(aaa in 1:2){     
      ab  <- rnorm(n*L)
      ab  <- init.ab
      opt <- optim(ab,fn=sse,gr=sse_grad,n=n,L=L,EC=EC,alpha=alpha.grid[aaa],fudge=0.001,
                   method="BFGS",control = list(trace=ifelse(verbose,50,0),maxit=maxit1,reltol=0.00000001))
      val[aaa]    <- opt$val
      plot(alpha.grid,val,pch=ifelse(val==min(val,na.rm=TRUE),19,1))
      
      if(opt$val<bestval){
        bestval <- opt$val
        init.ab <- opt$par
      }
    }
    alpha       <- alpha.grid[which.min(val)]
  }
  
  # Now estimate b
  
  opt <- optim(init.ab,fn=sse,gr=sse_grad,n=n,L=L,EC=EC,alpha=alpha,fudge=0.001,
               method="BFGS",control = list(trace=ifelse(verbose,50,0),maxit=maxit2,reltol=0.00000001))
  
  b     <- matrix(opt$par,n,L)
  B     <- sweep(exp(b), 1, rowSums(exp(b)), "/")
  
  pct   <- colSums(B)/sum(B)
  tock  <- proc.time()[3]
  
  output <- list(est = B, pct = pct, alpha = alpha,
                 EC.smooth = ECs,seconds = tock-tick)
  
  return(output)}

