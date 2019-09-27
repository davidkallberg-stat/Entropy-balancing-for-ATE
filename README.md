# Entropy-balancing-for-ATE
R code that implements the method in "Entropy Balancing for Average Causal Effects"


Example usage:

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
