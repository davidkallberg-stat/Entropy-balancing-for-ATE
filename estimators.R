# the follwoing packages are needed
library(ebal)
library(ATE)
library(PSW)

# Help functions #
# Linear equation solver that uses the Cholesky decomposition, used for the QR entropy balancing estimator
ch.solve <- function(A,b){backsolve(chol(A),forwardsolve(t(chol(A)), b))}

# To supress printing in large scale simulations #
#quiet <- function(x) { sink(tempfile()); on.exit(sink());invisible(force(x))} 

# Entropy balancing estimators. The function also estimates the corresponding OLS variance estimator and a nonparametric variance estimator #
# Kullback-Leibler divergence #
ebalKL <- function(data){
  
  xl <- colnames(data)[substr(colnames(data),1,1) == "X"]  
  X <- as.matrix(data[,xl])
  
  tr <- data$tr
  N <- nrow(X)
  
  w <- rep(1,N)
  k2 <- NA
  
  conv <- F
  
  c1 <- c2 <- rep(0, K+1)
  
  est.sd.np   <-  NA #nonparametric variance estimator
  
  try({
    
    ate <- ATE(Y = data$Y, X = X, Ti = tr)
    
    # weights
    w <- ifelse(tr, ate$weights.p, ate$weights.q)*N
    
    #estimation of variance
    form.or <- paste("Y",paste(xl, collapse=" + "),sep=" ~ ")
    
    or1 <- lm(form.or, data[tr==1,])
    or0 <- lm(form.or, data[tr==0,])
    
    m1 <- predict(or1, newdata = data)
    m0 <- predict(or0, newdata = data)
    
    k2 <- var(m1 - m0) + sigma(or1)^2*mean(data$tr*w^2) + sigma(or0)^2*mean((1-data$tr)*w^2)
    
    c1 <- -ate$lam.p
    c2 <- -ate$lam.q
    
    conv = T
    
    # estimated covariance matrix
    cm <- ate$vcov
    est.sd.np = sqrt(cm[1,1] + cm[2,2] - 2*cm[1,2])
    
  })
  
  est <- weighted.mean(data$Y, w = data$tr*w) - weighted.mean(data$Y, w = (1-data$tr)*w)
  # alt: est <- ate$est[3]
  
  list(w=w, est = est, est.sd = sqrt(k2/N), est.sd.np  = est.sd.np,  c1=c1, c2=c2, conv=conv)
  
}

# Quadratic Rényi diverence #
ebalQR <- function(data){
  
  xl <- colnames(data)[substr(colnames(data),1,1) == "X"]  
  tr <- data$tr
  N  <- length(tr)
  
  X <- cbind(1,as.matrix(data[xl]))
  
  # sum of complete sample
  b <- colSums(X)
  
  X1 <- X*tr # treated
  X0 <- X*(1-tr) # control
  
  A1 <- t(X1) %*% X1
  A0 <- t(X0) %*% X0
  
  # dual parameters
  
  k1 <- c(mean(tr),0*is.na(xl))
  k0 <- c(mean(1-tr),0*is.na(xl))
  
  conv <- F
  
  try({k1 <- ch.solve(A1,b); k0 <- ch.solve(A0,b); conv <- T})
  
  w1 <- X1 %*% k1
  w0 <- X0 %*% k0
  
  w <- w1 + w0
  
  est <- weighted.mean(data$Y,w=data$tr*w) - weighted.mean(data$Y,w=(1-data$tr)*w)
  
  #estimation of variance
  form.or <- paste("Y",paste(xl, collapse=" + "),sep=" ~ ")
  
  or1 <- lm(form.or, data[tr==1,])
  or0 <- lm(form.or, data[tr==0,])
  
  m1 <- predict(or1,newdata = data)
  m0 <- predict(or0,newdata = data)
  
  k2 <- var(m1 - m0) + sigma(or1)^2*mean(data$tr*w^2) + sigma(or0)^2*mean((1-data$tr)*w^2)
  
  form <- paste("Y",paste(xl, collapse=" + "),sep=" ~ ")
  
  # nonparametric variance estimator
  
    b1 <- lm(form, data, subset = tr==1)
    b0 <- lm(form, data, subset = tr==0)
    
    b1 <- predict(b1, newdata = data)
    b0 <- predict(b0, newdata = data)
  
  # estimate the efficient influence function
  eif <- w*Y*(2*tr - 1) - est - b1*(tr*w - 1) + b0*((1-tr)*w - 1)

  list(w = w, est = est, est.sd = sqrt(k2/N), est.sd.np = sqrt(mean(eif^2)/N), kappa1 = k1, kappa0 = k0, conv = conv)
  
}

# IPW logistic regression with estimate of standard deviation (sandwich matrix)
IPW <- function(data){
  
  xl <- colnames(data)[substr(colnames(data),1,1) == "X"]
  form.tr <- paste("tr",paste(xl, collapse=" + "),sep=" ~ ")
  ipw <- psw(data,form.tr, weight ="ATE",wt=T,out.var = "Y")
  return(list(w=ipw$W, est = ipw$est.wt, est.sd = ipw$std.wt, ps = ipw$ps.hat, conv = ipw$ps.model$converged))
  
}

 

