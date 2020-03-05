# Entropy-balancing-for-ATE
R code that implements the method in "Entropy Balancing for Average Causal Effects".
The following code requires that you load the script ebalanceATE.R.
We generate example data:

```R
N <- 1000
X1 <- rnorm(N)
X2 <- rnorm(N)
Ux <- data.frame(X1,X2)
ps <- 1/(1+ exp(1.5*X1-X2))
tr <- rbinom(n = N,1,prob = ps)
Y1 <- 2*X1 - X2 + rnorm(N,0,sd = 1)
Y0 <- X1 - .5*X2 + rnorm(N,0,sd = 1)
Y <- ifelse(tr == 1, Y1, Y0)
```

To estimate the entropy balancing weights, using either the KL-or QR divergence:
```R
ebKL <- ebalanceATE(Ux, tr, div="KL") # med KL-divergence
ebQR <- ebalanceATE(Ux, tr, div="QR") # med QR-divergence
```
The estimates for the average treatment effect are obtained as
```R
# ATE-estimates
ATE.kl <- weighted.mean(Y1[tr==1], w=ebKL$w1) - weighted.mean(Y0[tr==0], w=ebKL$w0)
ATE.qr <- weighted.mean(Y1[tr==1], w=ebQR$w1) - weighted.mean(Y0[tr==0], w=ebQR$w0)
```
To estimate the asymptotic variance from Equation ... in the paper
```R
# variance estimation
lm1 <- lm(Y[tr==1] ~ X1 + X2, data = Ux[tr==1,]) 
lm0 <- lm(Y[tr==0] ~ X1 + X2, data = Ux[tr==0,])

m1 <- predict(lm1, newdata = Ux)
m0 <- predict(lm0, newdata = Ux)
v.kl <- var(m1 - m0) + 
          sigma(lm1)^2*sum(ebKL$w1^2)/N + sigma(lm0)^2*sum(ebKL$w0^2)/N
```
Then we can obtain a confidence interval för ATE:
```R
ci <- ATE.kl + c(-1.96*sqrt(v.kl)/sqrt(N),1.96*sqrt(v.kl)/sqrt(N))

