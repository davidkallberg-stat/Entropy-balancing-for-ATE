####################################################################################
# R code for the NHANES data example on the effect of smoking of blood lead levels #
####################################################################################

library(xtable)
source('estimators.R')
set.seed(2021126)

# load data
Data <- read.csv("dataMale.csv")

# Covariates to be balanced
X.mat <- model.matrix(smoking ~ married + birth.country + edu + race + income + income.mis + cage + army + cfamily.size, data = Data)
df.X <- as.data.frame(X.mat)
names(df.X) <- paste0("X",names(df.X))

# Data frame
df <- data.frame(df.X[,-1], tr = Data$smoking, Y  = Data$lead)

# entropy balancing estimators based on the KL and QR divergences
kl <- ebalKL(df)
qr <- ebalQR(df)

#bootstrap variance estimators
N  <- nrow(df)
B <- 10000 # number of resamples 
boots <- data.frame(replicate(B,sample(N,N,T))) #

kl_boot <- lapply(1:B, function(x){ebalKL(df[boots[,x],])})
qr_boot <- lapply(1:B, function(x){ebalQR(df[boots[,x],])})

estsKL <- unlist(lapply(kl_boot, '[[', "est"))
estsQR <- unlist(lapply(qr_boot, '[[', "est"))

kl$est.sd.boot <- sd(estsKL)
qr$est.sd.boot <- sd(estsQR)

# Create result table
# Point estimates and standard errors
est <- c(rep(kl$est,3), rep(qr$est,3))
sd.est <- c(kl$est.sd,kl$est.sd.np,kl$est.sd.boot,qr$est.sd, qr$est.sd.np, qr$est.sd.boot)

# confidence intervals
ci_lower <- round(est - 1.96*sd.est,2)
ci_upper <- round(est + 1.96*sd.est,2)

ci <- paste0("(",ci_lower,", ",ci_upper,")")

# table
xtable(cbind(round(est,2),round(sd.est,2),ci))

