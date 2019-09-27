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
  
