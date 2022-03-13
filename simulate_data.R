library(ebal)

# Construct data-frame with 2nd order terms + interactions of covariates (the balancing vector U26 in the paper)
# For the real diabetes data set, real = T is used

cov.int <- function(data, real=F){
  
  xl <- colnames(data)[substr(colnames(data),1,1) == "X"]
  X <- data[,xl]
  
  p <- ncol(X)
  N <- nrow(X)
  
  X2 <- matrixmaker(as.matrix(X))
  id <- cumsum(c(p+1,2:p))
  id.rm <- apply(X,2,function(x) length(table(x)) == 2)
  
  if(any(id.rm)){
    if(!real){
      data.frame(X2[,-id[id.rm]],Y=data$Y, Y1 = data$Y1, Y0 = data$Y0, tr= data$tr)
    }else{
      data.frame(X2[,-id[id.rm]],Y=data$Y, tr= data$tr)
    }
  } else {
    if(!real){
      data.frame(X2,Y=data$Y, tr= data$tr)  
    }else{
      data.frame(X2,Y=data$Y, tr= data$tr)  
    }
  }
  
}

## Simulate covariate vector X = (X1,X2,X3,X4,X5,X6) ##
Xsim <- function(n){
  
  # data parameters
  m <- c(0,0,0) 
  Sig <- matrix(0,3,3) 
  diag(Sig) <- c(2,1,1) 
  Sig[2,1] <- Sig[1,2] <- 1  
  Sig[3,1] <- Sig[1,3] <- -1 
  Sig[2,3] <- Sig[3,2] <- -.5  
  
  X123 <- mvrnorm(n=n, mu=m, Sigma=Sig)
  x4 <- runif(n,-3,3) 
  x5 <- rchisq(n,df = 1)
  x6 <- rbinom(n = n,size = 1,prob = .5)
  X <- cbind(X123,x4,x5,x6)  
  colnames(X) <- paste('X',1:6,sep="")
  return(X)
}

### Simulate data. Data-frame with outcome Y, Y(1), Y(0), treatment, conditional mean ###
# or.desig is "a","b", or "c"
# tr.effect is zero (z) or heterogenous (h)

Datasim <- function(N, or.design = "a", tr.effect = "c"){
  
  # simulate confounder vector
  X <- Xsim(N)
  
  # coefficients, treatment variable
  beta <- c(1 , 2, -2, -1, -0.5, 1)
  
  #treatment variable
  eps <-  rnorm(N,0,sd = 30^.5)
  tr <- round(X %*% beta + eps > 0)
  
  # propensity score
  ps <- pnorm(X %*% beta, 0, sd = 30^.5)
  
  # conditional outcome
  
  # constants
  a1.0 <- c( .5, -.7, .6)
  a0.0 <- c(-.7, .6, -.4)
  
  # errors epst = Y(t) - bt, for t = 1,0
  eps1 <- rnorm(N)  
  eps0 <- rnorm(N) 
  
  if(or.design == "a"){
    
    RX <- as.matrix(cbind(1,X))
    
    # outcome coefficients
    
    if(tr.effect == "z" | tr.effect == "c"){
      
      a1 <- c(0, 1, 1, 1, -1, 1 ,1)
      a0 <- c(0, 1, 1, 1, -1, 1 ,1)
      
    } 
    
    if(tr.effect == "h"){
      
      a1 <- c( 0.5, 1.7, 1.5, 1.4, -0.5, 1.3, 2.1)
      a0 <- c(-0.7, 0.7, 0.3, 0.8, -1.6, 0.9, 0.5)
      
    }
  }
  
  if(or.design == "b"){
    
    RX <- as.matrix(cbind(1, X[,1:2],X[,3]*X[,4],X[,5]^0.5))  
    
    if(tr.effect == "z" | tr.effect == "c"){
      
      a1 <- c(0, 1, 1, 0.2, -1)
      a0 <- c(0, 1, 1, 0.2, -1)
      
    }
    
    if(tr.effect == "h"){
      
      a1 <- c(-0.7, 1.7, 0.7, -0.2, -1.5)
      a0 <- c(0.6, 0.6, 0.5, 0.3, -0.5)
      
    }
  }
  
  if(or.design == "c"){
    
    RX <- as.matrix(cbind(1, X[,c(1:2,5)]^2,X[,1]*X[,2],X[,1]*X[,5],X[,2]*X[,5]))  
    
    if(tr.effect == "z" | tr.effect == "c"){
      
      a1 <- c(0,1,1,1,2,2,2)
      a0 <- c(0,1,1,1,2,2,2)
      
    }
    
    if(tr.effect == "h"){
      
      a1 <- c(0.6, 1.4, 1.5, 1.3, 2.5, 2.7, 2.2)
      a0 <- c(-0.4, 0.7, 0.5, 0.6, 1.5, 1.9, 1.6)
      
    }
    
  }
  
  #conditional mean functions
  b1 <- RX %*% a1
  b0 <- RX %*% a0
  
  # outcome
  Y1 <- b1 + eps1
  Y0 <- b0 + eps0
  
  Y <- ifelse(tr, Y1, Y0)
  
  colnames(X) <- paste0("X",1:6)
  
  # output
  data.frame(X, Y1, Y0, Y, tr, unadj = mean(Y1[tr==1])-mean(Y0[tr==0]), ps, b1, b0)
  
}
