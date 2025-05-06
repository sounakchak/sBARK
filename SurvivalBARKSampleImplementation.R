## Parametric and Nonparametric Survival Analysis with Bayesian Additive Regression Kernels (sBARK)

## Author: Sounak Chakraborty and Peng Zhao



library(dummies)
library(survival)
library(mlr)
library(BART)
library(survsim)
source("surv.bark.R")
source("aft.tree.R")

library(devtools)
#install_github("cran/bark")
library(bark)

######################################################## Real data #####################################################


## ============================================ data cleaning & preprocessing =========================================

# Kidney catheter data: recurrence times to infection, at the point of insertion of the catheter
# https://cran.r-project.org/web/packages/survival/survival.pdf p43

data(kidney)
kidney <- kidney[,-1]  # drop patient id
kidney.surv <- kidney[,c(1,2)]
kidney.x <- kidney[,c(3,4,5,6)]
kidney.x$sex <- as.factor(kidney.x$sex)
kidney.x.dummy <- mlr::createDummyFeatures(kidney.x[,-c(1,4)])
kidney.x.scale <- scale(kidney.x[,c(1,4)])
kidney.cleaned <- cbind(kidney.surv, kidney.x.dummy, kidney.x.scale)
kidney.uncensored <- kidney.cleaned[kidney.cleaned$status == 1,]
kidney.censored <- kidney.cleaned[kidney.cleaned$status == 0,]


## =========================================== nonparametric survival analysis =========================================

# discrete-time bark
bark.result <- surv.bark(x.train=as.matrix(kidney.cleaned[,-c(1,2)]), times.train=kidney.cleaned[,"time"],
                         delta.train=kidney.cleaned[,"status"], x.test=as.matrix(kidney.cleaned[,-c(1,2)]), times.test=kidney.cleaned[,"time"])

bark.K <- bark.result$K
bark.M <- ncol(bark.result$yhat.test)
bark.N <- length(kidney$status)

bark.surv <- matrix(nrow=bark.M, ncol=bark.K)

for(j in 1:bark.K) {
  h <- seq(j, bark.N*bark.K, by=bark.K)
  bark.surv[ , j] <- apply(bark.result$surv.test[h, ], 2, mean)
}

bark.surv.mu  <- apply(bark.surv, 2, mean)

# discrete-time bart
bart.result <- surv.bart(x.train=as.matrix(kidney.cleaned[,-c(1,2)]), times=kidney.cleaned[,"time"],
                         delta=kidney.cleaned[,"status"], x.test=as.matrix(kidney.cleaned[,-c(1,2)]))

bart.K <- bart.result$K
bart.M <- nrow(bart.result$yhat.test)
bart.N <- length(kidney$status)

bart.surv <- matrix(nrow=bart.M, ncol=bart.K)

for(j in 1:bart.K) {
  h <- seq(j, bart.N*bart.K, by=bart.K)
  bart.surv[ , j] <- apply(bart.result$surv.test[ ,h], 1, mean)
}

bart.surv.mu  <- apply(bart.surv, 2, mean)

# Kaplan-Meier (baseline)
km.result <- survfit(formula=Surv(kidney.cleaned[,"time"], kidney.cleaned[,"status"]) ~ 1, data=as.data.frame(kidney.cleaned[,-c(1,2)]))

# make plot
plot(c(0, bark.result$events), c(1, bark.surv.mu), type='s', col='red', ylim=0:1, ylab='Survival Probability', xlab='Time (days)', main='Kidney Catheter - Time to Infection')
lines(c(0, bart.result$times), c(1, bart.surv.mu), col='green', type='s')
lines(km.result$time, km.result$surv, col='blue', type='s')
legend(x="topright", legend=c('Survival BARK', 'Survival BART', 'Kaplan-Meier'), col=c('red', 'green','blue'), lty=1)


## =========================================== parametric survival analysis =========================================

cv_fold = 5

# shuffle the data 
set.seed(.Random.seed)
indices <- sample(1:nrow(kidney.uncensored))  

# split indices into cv_fold chunks
indices_chunks <- split(indices, ceiling(seq_along(indices)/(nrow(kidney.uncensored)/cv_fold)))  

errors <- matrix(ncol=2,nrow=cv_fold)

for (c in 1:cv_fold) {
  train_set <- kidney.uncensored[unlist(indices_chunks[-c]),]
  train_set <- rbind(train_set,kidney.censored)
  x.train <- as.matrix(train_set[,!names(train_set) %in% c("time","status")])
  test_set <- kidney.uncensored[unlist(indices_chunks[c]),]
  x.test <- as.matrix(test_set[,!names(test_set) %in% c("time","status")])
  
  bark.result <- surv.bark(x.train=x.train, times.train=train_set[,"time"], surv.type = "parametric",
                             delta.train=train_set[,"status"], x.test=x.test, times.test=test_set[,"time"])
  errors[c,1] <- bark.result$test.error
  
  y.train <- train_set[,c("time","status")]
  names(y.train) <- c("Survival", "Status")
  bart.result <- aft.tree(X.train=x.train,Y.train=y.train,X.test=x.test,n.iter=100,ntree=100,burn.in=20)
  bart.yhat.test.mean <- apply(bart.result$predictions$yhat.test, 2, mean)
  errors[c,2] <- sqrt(mean((bart.yhat.test.mean - log(test_set[,"time"]))^2))
  
}

errors.mean <- apply(errors, 2, mean)


######################################################## Simulated data #####################################################


## ================================================ generate simulated data  ================================================

# 80 instances
# maximum follow-up time is 365 days
# time to event follows log-normal distribution
# 7 covariates with each following a standard normal distribution

sim.data <- simple.surv.sim(n=80, foltime=365, dist.ev="lnorm", anc.ev=0.6, beta0.ev=5, anc.cens=1.2, beta0.cens=7, dist.cens="lnorm", beta=c(0.1,0.2,0.3,0.4,0.5,0.6,0.7),
          x=list(c("normal", 0, 1), c("normal", 0, 1), c("normal", 0, 1), c("normal", 0, 1), c("normal", 0, 1), c("normal", 0, 1), c("normal", 0, 1)))

sim.data <- sim.data[,-c(1,3,5)]
names(sim.data)[2] <- "time"
sim.data$time <- as.integer(sim.data$time)

sim.data.uncensored <- sim.data[sim.data$status == 1,]
sim.data.censored <- sim.data[sim.data$status == 0,]

# export the simulated data into csv file
write.csv(sim.data, file = "sim_data1.csv", row.names=FALSE)


## =========================================== nonparametric survival analysis =========================================

# discrete-time bark
sim.bark.result <- surv.bark(x.train=as.matrix(sim.data[,-c(1,2)]), times.train=sim.data[,"time"],
                         delta.train=sim.data[,"status"], x.test=as.matrix(sim.data[,-c(1,2)]), times.test=sim.data[,"time"])

sim.bark.K <- sim.bark.result$K
sim.bark.M <- ncol(sim.bark.result$yhat.test)
sim.bark.N <- length(sim.data$status)

sim.bark.surv <- matrix(nrow=sim.bark.M, ncol=sim.bark.K)

for(j in 1:sim.bark.K) {
  h <- seq(j, sim.bark.N*sim.bark.K, by=sim.bark.K)
  sim.bark.surv[ , j] <- apply(sim.bark.result$surv.test[h, ], 2, mean)
}

sim.bark.surv.mu  <- apply(sim.bark.surv, 2, mean)

# discrete-time bart
sim.bart.result <- surv.bart(x.train=as.matrix(sim.data[,-c(1,2)]), times=sim.data[,"time"],
                         delta=sim.data[,"status"], x.test=as.matrix(sim.data[,-c(1,2)]))

sim.bart.K <- sim.bart.result$K
sim.bart.M <- nrow(sim.bart.result$yhat.test)
sim.bart.N <- length(sim.data$status)

sim.bart.surv <- matrix(nrow=sim.bart.M, ncol=sim.bart.K)

for(j in 1:sim.bart.K) {
  h <- seq(j, sim.bart.N*sim.bart.K, by=sim.bart.K)
  sim.bart.surv[ , j] <- apply(sim.bart.result$surv.test[ ,h], 1, mean)
}

sim.bart.surv.mu  <- apply(sim.bart.surv, 2, mean)

# Kaplan-Meier (baseline)
sim.km.result <- survfit(formula=Surv(sim.data[,"time"], sim.data[,"status"]) ~ 1, data=as.data.frame(sim.data[,-c(1,2)]))

# make plot
plot(c(0, sim.bark.result$events), c(1, sim.bark.surv.mu), type='s', col='red', ylim=0:1, ylab='Survival Probability', xlab='Time (days)', main='Simulated Data')
lines(c(0, sim.bart.result$times), c(1, sim.bart.surv.mu), col='green', type='s')
lines(sim.km.result$time, sim.km.result$surv, col='blue', type='s')
legend(x="topright", legend=c('Survival BARK', 'Survival BART', 'Kaplan-Meier'), col=c('red', 'green','blue'), lty=1)

## =========================================== parametric survival analysis =========================================

sim.cv_fold = 5

# shuffle the data 
sim.indices <- sample(1:nrow(sim.data.uncensored))  

# split indices into cv_fold chunks
sim.indices_chunks <- split(sim.indices, ceiling(seq_along(sim.indices)/(nrow(sim.data.uncensored)/sim.cv_fold)))  

sim.errors <- matrix(ncol=2,nrow=sim.cv_fold)

for (c in 1:sim.cv_fold) {
  train_set <- sim.data.uncensored[unlist(sim.indices_chunks[-c]),]
  train_set <- rbind(train_set,sim.data.censored)
  x.train <- as.matrix(train_set[,!names(train_set) %in% c("time","status")])
  test_set <- sim.data.uncensored[unlist(sim.indices_chunks[c]),]
  x.test <- as.matrix(test_set[,!names(test_set) %in% c("time","status")])
  
  bark.result <- surv.bark(x.train=x.train, times.train=train_set[,"time"], surv.type = "parametric",
                           delta.train=train_set[,"status"], x.test=x.test, times.test=test_set[,"time"])
  sim.errors[c,1] <- bark.result$test.error
  
  y.train <- train_set[,c("time","status")]
  names(y.train) <- c("Survival", "Status")
  bart.result <- aft.tree(X.train=x.train,Y.train=y.train,X.test=x.test,n.iter=100,ntree=100,burn.in=20)
  bart.yhat.test.mean <- apply(bart.result$predictions$yhat.test, 2, mean)
  sim.errors[c,2] <- sqrt(mean((bart.yhat.test.mean - log(test_set[,"time"]))^2))
  
}

sim.errors.mean <- apply(sim.errors, 2, mean)

# =========================================== nonparametric survival analysis =========================================

# discrete-time bark
bark.result <- surv.bark(x.train=as.matrix(kidney.cleaned[,-c(1,2)]), times.train=kidney.cleaned[,"time"],
                         delta.train=kidney.cleaned[,"status"], x.test=as.matrix(kidney.cleaned[,-c(1,2)]), times.test=kidney.cleaned[,"time"])




########################################################
################# hospital readmit data ################
########################################################

source("surv.bark.R")
source("aft.tree.R")
library(rlist)

data <- read.csv(file="readmission_cleaned.csv", header=TRUE, sep=",")

data.uncensored <- data[data$status == 1,]
data.censored <- data[data$status == 0,]  

train_set <- data 

x.train <- as.matrix(data[,!names(data) %in% c("time","status")]) 
x.test <- x.train

cc <- c(2,5,19,41,42,60,71,78,83,84)

x.train1 <- x.train[,cc] 
x.test1 <- x.test[,cc]

y.train <- data[,c("time","status")]
names(y.train) <- c("Survival", "Status")

## =========================================== nonparametric survival analysis =========================================

# discrete-time bark # NOT YET WORKING WITH READMIN DATA

sim.bark.result <- surv.bark(x.train=x.train1) times.train=y.train[,"time"],
                             delta.train=y.train[,"status"], x.test=x.train1, times.test=y.train[,"time"])

bark.result <- surv.bark(x.train=x.train, times.train=train_set[,"time"], surv.type="parametric", delta.train=train_set[,"status"], 
                         x.test=x.test, times.test=test_set[,"time"], keepevery=1, nburn=2, nkeep=10, p.mcmc.iter=100, p.nburn=20)

sim.bark.K <- sim.bark.result$K
sim.bark.M <- ncol(sim.bark.result$yhat.test)
sim.bark.N <- length(sim.data$status)

sim.bark.surv <- matrix(nrow=sim.bark.M, ncol=sim.bark.K)

for(j in 1:sim.bark.K) {
  h <- seq(j, sim.bark.N*sim.bark.K, by=sim.bark.K)
  sim.bark.surv[ , j] <- apply(sim.bark.result$surv.test[h, ], 2, mean)
}

sim.bark.surv.mu  <- apply(sim.bark.surv, 2, mean)

# discrete-time bart
sim.bart.result <- surv.bart(x.train=x.train1, times=y.train[,"time"],
                             delta=y.train[,"status"], x.test=x.train1)

sim.bart.K <- sim.bart.result$K
sim.bart.M <- nrow(sim.bart.result$yhat.test)
sim.bart.N <- length(sim.data$status)

sim.bart.surv <- matrix(nrow=sim.bart.M, ncol=sim.bart.K)

for(j in 1:sim.bart.K) {
  h <- seq(j, sim.bart.N*sim.bart.K, by=sim.bart.K)
  sim.bart.surv[ , j] <- apply(sim.bart.result$surv.test[ ,h], 1, mean)
}

sim.bart.surv.mu  <- apply(sim.bart.surv, 2, mean)

# Kaplan-Meier (baseline)
sim.km.result <- survfit(formula=Surv(y.train[,"time"], y.train[,"status"]) ~ 1, data=as.data.frame(x.train1))

# make plot
#plot(c(0, sim.bark.result$events), c(1, sim.bark.surv.mu), type='s', col='red', ylim=0:1, ylab='Survival Probability', xlab='Time (days)', main='Simulated Data')
plot(c(0, sim.bart.result$times), c(1, sim.bart.surv.mu), type='s', col='red', ylab='Survival Probability', xlab='Time (days)', main='Simulated Data')
lines(c(0, sim.bart.result$times), c(1, sim.bart.surv.mu), col='green', type='s')
lines(sim.km.result$time, sim.km.result$surv, col='blue', type='s')
legend(x="topright", legend=c('Survival BARK', 'Survival BART', 'Kaplan-Meier'), col=c('red', 'green','blue'), lty=1)
