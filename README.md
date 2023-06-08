# Entropy-balancing-for-ATE

R code implementing the simulations in the paper "Large Sample Properties of Entropy Balancing Estimators of Average Causal Effects".
The following code requires that you load the script simulate_data.R (for generating the covariates, treatment and outcome) and estimators.R, and rely on some packages:

```R
install.packages(c("MASS","parallel","ebal","ATE","PSW"))
library(MASS)
library(parallel)
library(ebal)
library(ATE)
library(PSW) 

source("estimators.R") # for generating covariates, treatment and outcome
source("simulate_data.R") # entropy balancing estimators and corresponding variance estimators

```

Script for generating data and computing estimators

```R
N <- 500 # sample size
Nsim <- 1000 # number of simulated datasets
B <- 10000; Bextra <- B + 2000 # number of resamples for bootstrap variance estimators
nc <- detectCores() # number of cpu cores to be used, used by mclapply (parallell computing)
k <- 1

DataList <- list()

for(des in c("a", "b", "c")){
    for(te in c("z", "h")){
    
    # Generate Nsim datasets for Design des, treatment effect te (z(ero) or (h)eterogeneous)
    Data <- lapply(1:Nsim, function(x){Datasim(N, tr.effect = te, or.design = des)})
    
    # Collect data
    DataList[[k]] <- Data
    
    k <- k + 1 

# Estimates of the average treatment effect (beta)
    labs1 <- c("kl1", "qr1", "kl2", "qr2", "IPW1", "IPW2", "IPW3", "AIPW")
    betahat <- data.frame(matrix(ncol = length(labs1), nrow = Nsim, dimnames = list(NULL, labs1)))

# Variance estimators (parametric, nonparametric, bootstrap)
    labs2 <- c("kl1","qr1","kl2","qr2","kl1_np","qr1_np","kl2_np","qr2_np", "kl1_boot", "qr1_boot", "kl2_boot", "qr2_boot")
    varests <- data.frame(matrix(ncol = length(labs2), nrow = Nsim, dimnames = list(NULL, labs2)))

# For tracking that enough bootstrap samples have a solution
  boot_conv <- betahat[,1:4]

for(i in 1:Nsim){

    data <- Data[[i]]
    N <- nrow(data)

    X <- data[substr(names(data),1,1)=="X"]
    tr <- data$tr
    Y <- data$Y

  # 2nd order + interactions
    data2 <- cov.int(data, real = T)
    X2 <- data2[substr(names(data2),1,1)=="X"]

  # IPW-estimators, probit link
  e1 <- glm(tr~X[,1]+X[,2]+X[,3]+X[,4]+X[,5]+X[,6], family=binomial(link = "probit"))
  e2 <- glm(tr~I(X[,1]^2)+I(X[,2]^2)+X[,3]+I(X[,4]^2)+I(X[,5]^2)+X[,6], family=binomial(link = "probit"))
  e3 <- glm(tr~I(X[,1]*X[,3])+I(X[,2]^2)+X[,4]+X[,5]+X[,6], family=binomial(link = "probit"))

  ipw1 <- weighted.mean(Y,w = tr/e1$fitted.values) - weighted.mean(Y,w = (1-tr)/(1-e1$fitted.values))
  ipw2 <- weighted.mean(Y,w = tr/e2$fitted.values) - weighted.mean(Y,w = (1-tr)/(1-e2$fitted.values))
  ipw3 <- weighted.mean(Y,w = tr/e3$fitted.values) - weighted.mean(Y,w = (1-tr)/(1-e3$fitted.values))
    
  # Augmented IPW estimator, logit link, and linear outcome regression
  form.tr <- paste("tr", paste(xl, collapse=" + "), sep=" ~ ")
  form.or <- paste("Y", paste(xl, collapse=" + "), sep=" ~ ")
  
  aipw_fit <- psw.aug(data, form.tr, "ATE", form.or)  

  # KL estimator
  kl1 <- ebalKL(data)
  kl2 <- ebalKL(data2)

  # QR estimators
  qr1 <- ebalQR(data)
  qr2 <- ebalQR(data2)

  betahat$kl1[i] <- kl1$est
  betahat$qr1[i] <- qr1$est
  betahat$kl2[i] <- kl2$est
  betahat$qr2[i] <- qr2$est

  betahat$IPW1[i] <- ipw1
  betahat$IPW2[i] <- ipw2
  betahat$IPW3[i] <- ipw3
    
  betahat$AIPW[i] <- aipw_fit$est_aug   

  # Variance estimators
  # Parametric (proposed in the paper)
  varests$kl1[i] <- kl1$est.sd^2
  varests$qr1[i] <- qr1$est.sd^2
  varests$kl2[i] <- kl2$est.sd^2
  varests$qr2[i] <- qr2$est.sd^2

  # Non-parametric estimator from Chan
  varests$kl1_np[i] <- kl1$est.sd.np^2
  varests$qr1_np[i] <- qr1$est.sd.np^2
  varests$kl2_np[i] <- kl2$est.sd.np^2
  varests$qr2_np[i] <- qr2$est.sd.np^2

  # bootstrap variances
  qr1_boot <- mclapply(1:Bextra, function(x){ebalQR(data[sample(N,N,T),])}, mc.cores = nc)
  kl1_boot <- mclapply(1:Bextra, function(x){ebalKL(data[sample(N,N,T),])}, mc.cores = nc)
  kl2_boot <- mclapply(1:Bextra, function(x){ebalKL(data2[sample(N,N,T),])}, mc.cores = nc)
  qr2_boot <- mclapply(1:Bextra, function(x){ebalQR(data2[sample(N,N,T),])}, mc.cores = nc)

  boot_conv_kl1 <- unlist(lapply(kl1_boot, '[[', "conv"))
  boot_conv_qr1 <- unlist(lapply(qr1_boot, '[[', "conv"))
  boot_conv_kl2 <- unlist(lapply(kl2_boot, '[[', "conv"))
  boot_conv_qr2 <- unlist(lapply(qr2_boot, '[[', "conv"))

  est_boot_kl1 <- unlist(lapply(kl1_boot, '[[', "est"))[boot_conv_kl1][1:B]
  est_boot_qr1 <- unlist(lapply(qr1_boot, '[[', "est"))[boot_conv_qr1][1:B]
  est_boot_kl2 <- unlist(lapply(kl2_boot, '[[', "est"))[boot_conv_kl2][1:B]
  est_boot_qr2 <- unlist(lapply(qr2_boot, '[[', "est"))[boot_conv_qr2][1:B]
  
  boot_conv$kl1[i] <- sum(!is.na(est_boot_kl1)) 
  boot_conv$kl2[i] <- sum(!is.na(est_boot_kl2)) 
  boot_conv$qr1[i] <- sum(!is.na(est_boot_qr1)) 
  boot_conv$qr2[i] <- sum(!is.na(est_boot_qr2)) 

  varests$kl1_boot[i] <- var(est_boot_kl1, na.rm = T)
  varests$qr1_boot[i] <- var(est_boot_qr1, na.rm = T)
  varests$kl2_boot[i] <- var(est_boot_kl2, na.rm = T)
  varests$qr2_boot[i] <- var(est_boot_qr2, na.rm = T)
  
  print(i)

}

# save result for des, te
Res <-  list(betahat = betahat, varests = varests, boot_conv)
save(Res, file=paste0("Res_test",N, des,te,".RDa"))

  }
}

```
# Tables 

Which design should be analyzed?
```R
des <- "a" # design (a, b, or c)
te <- "z" # treatment effect (zero or heterogeneous)
load(file = paste0("Res_",N, des, te,".RDa")) # Res-file has been saved in loop above
```

True values of the average treatment effect for the different scenarios
```R
ate_tr_df <- c(az = 0, bz = 0, cz = 0, 
               ah  =  1.2 + 0.4*1 + 1.6 *.5, 
               bh = -1.3 - (2/pi)^.5, 
               ch = 1 + 0.7*2+1+0.7*3+1)

ate_tr <- unname(ate_tr_df[paste0(des,te)])
```
Bias-sd-rmse- table
```R

betahat <- Res$betahat
varests <- Res$varests
varests$kl2_np[rm_kl2np] <- NA

res <- betahat - ate_tr
bias <- colMeans(res)
std <- apply(res, 2, sd)
rmse <- colMeans(res^2)^.5

tab <- cbind(bias, std, rmse)
tab

```
# Evaluations of variance estimators

```R
betaV <- betahat[,rep(1:4, times = 3)]
v_emp <- apply(betaV, 2, var, na.rm=T) # empirical (Monte-Carlo) variance
Rn <- colMeans(varests, na.rm=T)/v_emp # ratio of estimated to empirical
```
Approximate standard error of Rn 
```R
d1 <- (colMeans(varests, na.rm=T))^2/(v_emp^2)/N
d2 <- apply(betaV, 2, kurtosis, na.rm=T)
d3 <- 0
d3 <- sapply(1:length(Rn), function(j) cov((betaV[,j] - mean(betaV[,j], na.rm=T))^2, varests[,j],  use = "na.or"))/v_emp/colMeans(betaV^2, na.rm=T)
d4 <- apply(varests, 2, var, na.rm=T)

Rn_sd = sqrt(d1*(d2 - 1 - 2*d3 + d4)) # 
sd_v <- paste0(round(Rn,2),"(", round(Rn_sd,2), ")")
tabV <- cbind(Rn, Rn_sd)
tabV
```

# Coverage probabilites (for 95% confidence intervals)

```R
coverage <- colMeans(abs(res[,c(1:4,1:4,1:4)]) < 1.96*varests^.5, na.rm = T)
N_SD <- apply(betahat[,c(1:4,1:4,1:4)],2,sd)*sqrt(N)
mean_est <- colMeans(varests^.5, na.rm=T)*sqrt(N)
  
tabC <- cbind(coverage, N_SD, mean_est)
tabC
```  
# Asymptotic error table, Corollary 1, Table 2
```R
error <- matrix(NA, nrow = 4, ncol = 10)

for(des in c("a","b","c")){
  for(tr.eff in c("z","h")){
    for(i in 1:10){
      
      d1 <- Datasim(1e6, or.design = des, tr.effect = tr.eff)
      d2 <- Datasim(1e6, or.design = des, tr.effect = tr.eff)
      
      X <- cbind(1,as.matrix(d2[,1:6])) # Balance vector U6
      X2 <- cbind(1,as.matrix(cov.int(d2)[,1:26])) # Balance vector U26
      
      ebKL1 <- ebalKL(d1)
      ebKL2 <- ebalKL(cov.int(d1))
      
      ebQR1 <- ebalQR(d1)
      ebQR2 <- ebalQR(cov.int(d1))
      
      # Errors from Corollary 1
      error[1,i] <- cov(d2$ps*exp(X %*% ebKL1$c1), d2$b1) - cov((1-d2$ps)*exp(X %*% ebKL1$c2), d2$b0)
      error[2,i] <- cov(d2$ps*(X %*% ebQR1$kappa1), d2$b1) - cov((1-d2$ps)*(X %*% ebQR1$kappa0), d2$b0)
      error[3,i] <- cov(d2$ps*exp(X2 %*% ebKL2$c1), d2$b1) - cov((1-d2$ps)*exp(X2 %*% ebKL2$c2), d2$b0)
      error[4,i] <- cov(d2$ps*(X2 %*% ebQR2$kappa1), d2$b1) - cov((1-d2$ps)*(X2 %*% ebQR2$kappa0), d2$b0)
      
      est <- round(rowMeans(error),3)
      est_sd <- round(apply(error,1,sd)/sqrt(10),3)
      
      assign(paste0("error", des, tr.eff), paste0(est,"(",est_sd, ")"))
      print(des)
      
    }
  }
}

error_table <- xtable(rbind(cbind(erroraz,errorbz, errorcz),cbind(errorah,errorbh, errorch)))
error_table

```



