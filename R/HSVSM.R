#' HSVS-M for Multi-SNP Multi-Trait Association Test
#'
#' Test for association between a set of SNPs and multiple traits using summary statistics.
#'
#' @param Z A p*k numeric matrix for the summary statistics of p SNPs and k traits. 
#' @param cormat A k*k numeric matrix for the correlation matrix for k traits. When NULL, an empirical correlation matrix will be calculated using Z, i.e., cormat=cor(Z).
#' @param burnin A non-negative integer for the number of burnin MCMC iterations. Default is 1000.
#' @param N A positive integer for the number of non-burnin MCMC iterations. Default is 2000.
#' @param priors A list of values for hyperparameters. When NULL, default values will be used.
#' @param initials A list of values for the initial values for MCMC sampling. When NULL, default values will be used.
#' @return A list that contains:
#' \describe{
#'   \item{p}{HSVS-M's Bayesian p-value.}
#'   \item{weight}{A data frame for HSVS-M's weight estimates for each trait.}
#' }
#' @import tmvn
#' @import statmod 
#' @export
#' @examples
#' data(HSVSM.example)
#' out<-HSVSM(HSVSM.example)
#' HSVSM.p<-out$p
#' HSVSM.weight<-out$weight
HSVSM<-function(Z,cormat=NULL,burnin=1000,N=2000,priors=list(NULL),initials=list(NULL)) { #1
  #library(MASS)
  Z<-as.matrix(Z)
  if (is.null(cormat)) cormat<-cor(Z)
  z<-as.vector(t(Z))
  p<-nrow(Z) #p: number of SNPs
  k<-ncol(Z) #k: number of traits
  kp<-k*p #length of z
  # Specify priors
  if(class(priors) != "list") stop("the argument priors should be a list")
  
  if(is.null(priors$alpha)) { 	# para. in prior beta dist for alpha
    a=1; b=1
    #a=1; b=240.7174 
  } else { 
    a = priors$alpha[1]; b = priors$alpha[2]
  }
  
  if(is.null(priors$lambda1)) { # para. in prior gamma dist for lambda1
    aa1 = 1; bb1 = 2
    #aa1 = 0.1; bb1 = 0.1 #original: 1 2
  } else {
    aa1 = priors$lambda1[1]; bb1 = priors$lambda1[2]
  }
  
  if(is.null(priors$lambda2)) { # para. in prior gamma dist for lambda2
    aa2 = 1; bb2 = 2
    #aa2 = 0.1; bb2 = 0.1 #original: 1 2
  } else {
    aa2 = priors$lambda2[1]; bb2 = priors$lambda2[2]
  }
  
  # Specify initial values for MCMC sampling
  if(class(initials) != "list") stop("the argument initials should be a list")
  
  if(is.null(initials$beta)) 	beta = rep(0,p) 	else 
    if(length(initials$beta) != p) stop("the initials of beta should be a px1 vector") else beta = initials$beta
  
  if(is.null(initials$w)) 	w = rep(1,k) 	else #w = z
    if(length(initials$w) != k) stop("the initials of beta should be a kx1 vector") else w = initials$w
  
  if(is.null(initials$sigma)) 	sigma = 1 	else 
    if(length(initials$sigma)!=1) stop("the initial of sigma should be a positive scalar") else
      if(initials$sigma <= 0) stop("the initial of sigma should be a positive scalar") else sigma = initials$sigma
  
  gamma<-0
  
  if(is.null(initials$alpha)) 	alpha = 0.5  	else
    if(length(initials$alpha)!=1) stop("the initial of alpha should be a scalar within (0,1)") else
      if(initials$alpha <=0 | initial$alpha >=1) stop("the initial of alpha should be a scalar within (0,1)") else alpha = initials$alpha
  
  if(is.null(initials$tau)) 	tau = rep(1,p) 	else
    if(length(initials$tau) != p) stop("the initials of tau should be a px1 vector of positive values") else
      if(any(initials$tau <= 0)) stop("the initials of tau should be a px1 vector of positive values") else tau = initials$tau
  
  if(is.null(initials$omega)) omega = rep(1,p-1) else
    if(length(initials$omega) != p-1) stop("the initials of omega should be a (p-1)x1 vector of positive values") else
      if(any(initials$omega <= 0)) stop("the initials of omega should be a (p-1)x1 vector of positive values") else omega = initials$omega
  
  if(is.null(initials$lambda1)) 	lambda1<-0.5 	else
    if(length(initials$lambda1) != 1) stop("the initials of lambda1 should be a positive scalar") else
      if(any(initials$lambda1 <= 0)) stop("the initials of lambda1 should be a positive scalar") else lambda1 = initials$lambda1
  
  if(is.null(initials$lambda2)) 	lambda2<-0.5 	else
    if(length(initials$lambda2) != 1) stop("the initials of lambda2 should be a positive scalar") else
      if(any(initials$lambda2 <= 0)) stop("the initials of lambda2 should be a positive scalar") else lambda2 = initials$lambda2
  
  #Construct Sig_Z^-1
  cormatinv<-solve(cormat) #Sig_Zi^-1
  sigzinv<-kronecker(diag(p),cormatinv) #Sig_Z^-1
  
  #B2<-eigenMapMatMult(eigenMapMatMult(t(z),sigzinv),z)
  B2<-t(z)%*%sigzinv%*%z
  
  #Construct matrix used in building invsig
  sigmat<-matrix(c(1,-1,-1,1),2,2)
  
  #get direction
  snpeffect<-apply(Z,1,mean)
  index<-which(abs(snpeffect)>=0.5)
  up<-rep(0,k)
  l<-length(index)
  if (l>0){
    if (l>1) phe<-apply(Z[index,],2,mean) else phe<-Z[index,]
    up[which(phe<=-0.5)]<-(-1) #t=0.5
    up[which(phe>=0.5)]<-1
  }
  lo<-rep(0,k)
  ind<-which(up==0)
  lo[ind]<-(-1)
  up[ind]<-1
  ind<-which(up==-1)
  lo[ind]<-(-1)
  up[ind]<-0
  
  lo1<-lo+0.001
  up1<-up-0.001
  D<-diag(k)
  diagp<-diag(p)
  b0<-rep(0,p)
  
  # Define vectors/matrices to store MCMC samples
  MCMC.beta = array(NA,dim=c(p,N))
  MCMC.w = array(NA,dim=c(k,N))
  MCMC.sigma = array(NA,N)
  MCMC.gamma = array(NA,N)
  MCMC.alpha = array(NA,N)
  MCMC.tau = array(NA,dim=c(p,N))
  MCMC.omega = array(NA,dim=c(p-1,N))
  MCMC.lambda1 = array(NA,N)
  MCMC.lambda2 = array(NA,N)
  
  # Monte Carlo Simulation Process
  
  for (kk in 1:(burnin+N)) { #2
    
    W<-kronecker(diagp,w)
    Wsig<-t(W)%*%sigzinv #Wsig<-eigenMapMatMult(t(W),sigzinv)
    WsigW<-Wsig%*%W #WsigW<-eigenMapMatMult(Wsig,W)
    mu<-Wsig%*%z #mu<-eigenMapMatMult(Wsig,z)
    
    # Update Beta
    if(p==1) {       						# gene has only one SNP
      invsig <- 1/tau
      sig<-1/( invsig + WsigW )
      sigmu<-sig*mu
      log.c <- -log(tau)/2+log(sig)/2+0.5*mu*sig*mu/sigma
    } else {      							# gene has more than one SNP
      invsig <- diag(1/tau)
      for(i in 1:(p-1)) invsig[i:(i+1),i:(i+1)]<-invsig[i:(i+1),i:(i+1)]+1/omega[i]*sigmat
      sig<-solve(invsig+WsigW)
      sigmu<-sig%*%mu #sigmu<-eigenMapMatMult(sig,mu)
      log.c = log(det(invsig))/2+log(det(sig))/2+0.5*t(mu)%*%sigmu/sigma
      #log.c = log(getDeterminant(invsig))/2+log(getDeterminant(sig))/2+0.5*eigenMapMatMult(t(mu),sigmu)/sigma
    }
    pg <- (1-alpha)/(1-alpha+alpha*exp(log.c))	# probility of zero
    if (is.na(pg)) pg<-0
    if (runif(1)<pg) {
      beta <- b0 #b0=rep(0,p)
      gamma <- 0
    } else {
      #beta[i] = mvrnorm(1,sig%*%mu,sig*sigma); gamma[g] = 1
      beta <- as.vector(t(chol(sig*sigma))%*%rnorm(p)+sigmu) #beta <- as.vector(eigenMapMatMult(t(chol(sig*sigma)), rnorm(p))+sigmu)
      gamma <- 1
    }
    
    #Update w
    if (gamma==0) w<-runif(k,lo,up) else{
      ind<-which(abs(beta)>=0.5)
      if (length(ind)>0){
        sqbeta<-as.vector(beta[ind]%*%beta[ind])
        sbeta<-as.vector(beta[ind]%*%Z[ind,])
        tmp<-sbeta/sqbeta
        ini<-tmp
        ind<-which(ini<=lo)
        ini[ind]<-lo1[ind]
        ind<-which(ini>=up)
        ini[ind]<-up1[ind]
        w<-rtmvn(n=1,Mean=tmp,Sigma=sigma/sqbeta*cormat,D=D,lower=lo,upper=up,init=ini)
        w<-as.vector(w)
      } else w<-runif(k,lo,up)
    }
    
    # Update Sigma
    if (gamma==0){
      A<-0
      B<-B2
    }else{
      if (p==1) A<-invsig*(beta^2) else A<-t(beta)%*%invsig%*%beta #A<-eigenMapMatMult(eigenMapMatMult(t(beta),invsig),beta)
      B<-z-W%*%beta #B<-z-eigenMapMatMult(W, beta)
      B<-t(B)%*%sigzinv%*%B #B<-eigenMapMatMult(eigenMapMatMult(t(B),sigzinv),B)
    }
    sigma <- 1/rgamma(1,(kp+gamma*p)/2,rate=(A+B)/2)
    
    # Update Alpha
    alpha <- rbeta(1, gamma+a, 1-gamma+b)
    
    # Update Tau
    if (gamma==0) tau <- rexp(p,rate=1/2) else
      for (i in 1:p) tau[i] <- 1/rinvgauss(1,sqrt(lambda1*sigma/beta[i]^2),lambda1)
    
    # Update Omega
    if (p>1){
      if (gamma==0) omega <- rexp(p-1,rate=1/2) else
        for (i in 1:(p-1)) omega[i] <- 1/rinvgauss(1,sqrt(lambda2*sigma/(beta[i]-beta[i+1])^2),lambda2)
    }
    
    # Update Lambda1 & Lambda2
    lambda1 <- rgamma(1,p+aa1,rate=sum(tau)/2+bb1)
    if (p>1) lambda2 <- rgamma(1,p-1+aa2,rate=sum(omega)/2+bb2)
    
    # Store MCMC samples into vector/matrix
    if(kk > burnin) {
      kk2 = kk-burnin
      MCMC.beta[,kk2] <- beta
      MCMC.w[,kk2]<-w
      MCMC.sigma[kk2] <- sigma
      MCMC.gamma[kk2] <- gamma
      MCMC.alpha[kk2] <- alpha
      MCMC.tau[,kk2] <- tau
      MCMC.omega[,kk2] <- omega
      MCMC.lambda1[kk2] <- lambda1
      MCMC.lambda2[kk2] <- lambda2
    }
    
  } #2 end of MCMC sampling
  
  return.list<-list(MCMC.beta,MCMC.w,MCMC.sigma,MCMC.gamma,MCMC.alpha,MCMC.tau,MCMC.omega,MCMC.lambda1,MCMC.lambda2)
  return.p<-1-mean(return.list[[4]])
  return.weight<-c()
  for (i in 1:k){
    return.weight<-rbind(return.weight,c(i,mean(return.list[[2]][i,]),median(return.list[[2]][i,]),sd(return.list[[2]][i,]),unname(quantile(return.list[[2]][i,],probs=c(0.025,0.975)))))  
  }
  colnames(return.weight)<-c("Trait","Mean","Median","SD","CI0.025","CI0.975")
  out<-list()
  out$p<-return.p
  out$weight<-return.weight
  return(out)
} #1
