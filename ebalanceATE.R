# Entropy balancing for estimation of ATE
require(ebal)
ebalanceATE <- function(Ux,tr, div="KL"){

  K <- ncol(Ux)
  Ux <- cbind(1,as.matrix(Ux))
  Ux1 <- Ux[tr==1,] # balancing functions for treated
  Ux0 <- Ux[tr==0,] # balancing functions for control
  
  if(rank(Ux1) < ncol(Ux1) || rank(Ux0) < ncol(Ux0)) { 
    return('error: .... needs to have full rank...') 
  }
  
  if(div=="KL"){
    
    m <- colMeans(Ux[,-1]) # sum of combined sample
    U1 <- U0 <- Ux[,-1]
  
    U1[tr==0,] <- matrix(m, sum(1-tr), K, byrow=T)	
    U0[tr==1,] <- matrix(m, sum(tr), K, byrow=T)	
  
    w1 <- ebalance(1-tr, U1, constraint.tolerance = 0.1)
    w0 <- ebalance(tr, U0, constraint.tolerance = 0.1)
  
    return(list(w1 = w1$w, w0 = w0$w))
  
  } else {
    
    s <- colSums(Ux)
  
    kappa1 <- solve(t(Ux1) %*% Ux1, s)
    w1 <- Ux1 %*% kappa1 #treatment group 
    
    kappa0 <- solve(t(Ux0) %*% Ux0, s)
    w0 <- Ux0 %*% kappa0 #control group
    
    return(list(w1 = w1,w0 = w0))
  }
}
  
N <- 1000
X1 <- rnorm(N)
X2 <- rnorm(N)
Ux <- data.frame(X1,X2)
ps <- 1/(1+ exp(1.5*X1-X2))
tr <- rbinom(n = N,1,prob = ps)
Y1 <- 2*X1 - X2 + rnorm(N,0,sd = 1)
Y0 <- X1 - .5*X2 + rnorm(N,0,sd = 1)

ebKL <- ebalanceATE(Ux, tr, div="KL") # med KL-divergence
ebQR <- ebalanceATE(Ux, tr, div="QR") # med QR-divergence

# ATE-estimates
ATE.kl <- weighted.mean(Y1[tr==1], w=ebKL$w1) - weighted.mean(Y0[tr==0], w=ebKL$w0)

# variance estimation
lm1 <- lm(Y1[tr==1] ~ X1 + X2, data = Ux[tr==1,]) 
lm0 <- lm(Y0[tr==0] ~ X1 + X2, data = Ux[tr==0,])

m1 <- predict(lm1,df)
m0 <- predict(lm0,df)

v.kl <- var(m1 - m0) + sigma(lm1)^2*sum(ebKL$w1^2)/N + sigma(lm0)^2*sum(ebKL$w0^2)/N # Eq. ... in the paper

# obtain a confidence interval för ATE
ci <- ATE.kl + c(-1.96*sqrt(v.kl)/sqrt(N),1.96*sqrt(v.kl)/sqrt(N))
ci