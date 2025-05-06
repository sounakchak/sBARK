
## Parametric and Nonparametric Survival Analysis with Bayesian Additive Regression Kernels (sBARK)

## Author: Sounak Chakraborty and Peng Zhao



surv.bark <- function(
    x.train = matrix(0.0, 0L, 0L),   # Training data covariates, must be a matrix of doubles
    times.train = NULL,              # The time of event or right-censoring for training data
    delta.train = NULL,              # The event indicator: delta = 1: event, delta = 0: right-censoring
    x.test = matrix(0.0, 0L, 0L),    # Test data covariates, must be a matrix of doubles, should have the same structure as x.train
    times.test = NULL,               # The time of event or right-censoring for test data
    K = NULL,                        # K quantiles for time grid
    type = "e",                      # BARK type, e, d, se, or sd, default choice is e
                                         # e: BARK with equal weights
                                         # d: BARK with different weights
                                         # se: BARK with selection and equal weights
                                         # sd: BARK with selection and different weights
    keepevery = 10,                 # Every keepevery draw is kept 
    nburn = 20,                     # Number of MCMC iterations to be treated as burn in
    nkeep = 10,                     # Number of MCMC iterations to be kept for the posterior inference after burn in
    printevery = 10,               # As the MCMC runs, a message is printed every printevery draws
    keeptrain = TRUE,                # Logical, whether to keep results for training samples
    fixed = list(),                  # A list of fixed hyperparameters, using the default values if not sepecified
    tune = list(lstep=0.5, frequL=.2, dpow=1, upow=0, varphistep=.5, phistep=1),   # A list of tuning parameters, not expected to change
    theta = list(),                  # A list of the starting values for the parameter theta, use the defaults if nothing is given
    surv.type = "nonparametric",     # switch the type of survival analysis, default choice is nonparametric 
                                         # "nonparametric": discrete time nonparametric survival analysis
                                         # "parametric": accelerated failure time (AFT) model
    p.mcmc.iter = 100,               # Number of MCMC iterations to estimate omega, sigma, alpha in the parametric model
    p.nburn = 20                     # Number of MCMC iterations to be treated as burn in when estimating omega, sigma, alpha in the parametric model
) 

{
  
  start.time <- Sys.time()

# import libraries and functions
  
  library(bark)     # contains the original BARK algorithm
  library(MCMCglmm) # contains rtnorm(), the truncated normal distribution
  library(survival)
  
# nonparametric survival analysis
  
  if (surv.type == "nonparametric") {
    
    # Convert each (t, delta) pair into binary-coded y.train for all available observed event and censoring time points, see the report for details
    # Put t and the original x.train into the new tx.train 
    # The restructured data tx.train and y.train can be directly used in BARK

    re.struct <- surv.bark.restruct(times.train, delta.train, x.train, times.test, x.test, K=K)
    y.train   <- re.struct$y.train
    x.train   <- re.struct$tx.train
    x.test    <- re.struct$tx.test
    y.test    <- re.struct$y.test
    times.train  <- re.struct$times
    K         <- re.struct$K    # the number of unique times
    
    # fit a bark model with the restructured data
    
    bark.obj <- bark(x.train=x.train, y.train=y.train, x.test=x.test, type = type, keepevery = keepevery, nburn = nburn, 
                     nkeep = nkeep, printevery = printevery, keeptrain = keeptrain, fixed = fixed, tune = tune, theta = theta, classification=TRUE)
    
    # evaluate the model with test data if provided
    
    bark.obj$x.train <- x.train
    bark.obj$x.test <- x.test
    bark.obj$y.train <- y.train
    bark.obj$y.test <- y.test
    bark.obj$K <- K
    bark.obj$events <- times.train
    
    if (length(x.test) > 0) {

      # the number of different settings
      
      H <- nrow(x.test)/K 
      
      # the conditional probility of the event in the interval given no events in preceding intervals
      
      bark.obj$surv.test <- 1 - pnorm(bark.obj$yhat.test)
      for(h in 1:H) {
        for(j in 2:K) {
          l <- K*(h-1)+j
          bark.obj$surv.test[l,] <- bark.obj$surv.test[l-1,] * bark.obj$surv.test[l,]
        }
      }
      #bark.obj$surv.test.mean <- apply(bark.obj$surv.test, 1, mean)
      
      
      # get the misclassification error on test data
      
      bark.obj$prob.test <- pnorm(bark.obj$yhat.test)
      bark.obj$prob.test.mean <- apply(bark.obj$prob.test, 1, mean)
      bark.obj$test.error <- mean(((bark.obj$prob.test.mean > 0.5) != as.logical(y.test)))
    }
    
    return(bark.obj)
    
  # parametric survival analysis  
    
  } else if (surv.type == "parametric") {
    
    censored.id <- which(delta.train==0)

    # define containers for alpha, omega, and sigma
    
    sigma.error <- matrix(NA,p.mcmc.iter,1)
    alpha.posterior <- matrix(0,p.mcmc.iter+1,1)
    omega.train <- matrix(NA,p.mcmc.iter,length(times.train))
    omega <- matrix(NA,p.mcmc.iter,length(censored.id))
    omega.test<-matrix(NA,p.mcmc.iter,length(times.test))
    
    # log time
    log.time <- log(times.train)   
    
    # MCMC
    for (i in 1:p.mcmc.iter) {
      
      bark.obj <- bark(x.train = x.train, y.train = log.time - alpha.posterior[i], x.test = x.test, type = type, keepevery = keepevery, nburn = nburn, 
                       nkeep = nkeep, printevery = printevery, keeptrain = keeptrain, fixed = fixed, tune = tune, theta = theta, classification=FALSE)
      
      sigma.error[i]  <- 1/(bark.obj$theta.phi[nkeep]) 
      omega.train[i,] <- bark.obj$yhat.train.mean 
      omega.test[i,]  <- bark.obj$yhat.test.mean 
      alpha.posterior[i+1] <- mean(abs(log.time - omega.train[i,])) 
      
      # sample log time for censored cases from the truncated normal distribution
      
      for (j in 1:length(censored.id)) { 
        
        N.censored <- rtnorm(n = 1, mean = alpha.posterior[i] + omega.train[i,censored.id[j]], sd = sigma.error[i], lower = log.time[censored.id[j]], upper = Inf)
        log.time[censored.id[j]] <- N.censored                             
        omega[i,j] <- N.censored 
        
      }
      
      print(i)
      
    }
    
    # draw after burning in 
    
    draws <- seq(p.nburn,p.mcmc.iter,1)
    omega.train[,censored.id] <- omega
    yhat.train <- alpha.posterior[draws] + omega.train[draws,]
    yhat.train.mean <- apply(yhat.train, 2, mean)
    yhat.test  <- alpha.posterior[draws] + omega.test[draws,]
    yhat.test.mean <- apply(yhat.test, 2, mean)
    test.error <- sqrt(mean((yhat.test.mean - log(times.test))^2))
    
    end.time <- Sys.time()
    
    list(yhat.train = yhat.train, yhat.test = yhat.test, yhat.train.mean = yhat.train.mean, yhat.test.mean = yhat.test.mean,
         omega.train = omega.train[draws,], omega.test = omega.test[draws,], alpha = alpha.posterior[draws], 
         sigma.error = sigma.error[draws], test.error = test.error, duration = end.time - start.time)
    
  } else {
    
    stop("This program only supports nonparametric and parametric survival analysis.")
  
  }

}
  
  
  
## The following function is used to prepare the data for discrete-time nonparametric survival analysis in BARK
## It is created based on the function surv.pre.bart() in the BART package. The original function does not support 
## restructing y.test and it's implemented now. 

## Convert each (t, delta) pair into binary-coded y.train for all available observed event and censoring time points
## Put t and the original x.train into the new tx.train 
## The restructed data tx.train and y.train can be directly used in BARK

surv.bark.restruct <- function(
  train.times,    ## vector of survival times of training instances. contains both censored and uncensored cases
  train.delta,    ## vector of event indicators of training instances: 1 event, 0 right-censoring
  x.train=NULL,   ## matrix of covariate regressors of training instances
  test.times,     ## vector of survival times of test instances. contains only uncensored cases
  x.test=NULL,    ## matrix of covariate regressors of test instances
  K=NULL          ## if specified, then use K quantiles for time grid
) 
{
  
  N <- length(train.times)
  
  if(N!=length(train.delta))
    stop('The length of times and delta must be identical')
  
  if(length(x.train)>0 && N!=nrow(x.train))
    stop('The length of times and the number of rows in x.train, if any, must be identical')
  
  if(length(K)>0) {
    events <- unique(quantile(train.times, probs=(1:K)/K))
    attr(events, 'names') <- NULL
    
    for(i in 1:N) {
      k <- min(which(train.times[i]<=events))
      train.times[i] <- events[k]
    }
  }
  else {
    events <- unique(sort(train.times))
    ## time grid of events including censoring times
  }
  
  K <- length(events)
  
  if(events[1]<=0)
    stop('Time points exist less than or equal to time zero.')
  
  ## get y.train
  
  y.train <- integer(N)  # at least N long
  
  k <- 1
  
  for(i in 1:N) for(j in 1:K) if(events[j] <= train.times[i]) {
    y.train[k] <- train.delta[i]*(train.times[i] == events[j])
    
    k <- k+1
  }
  
  ## get tx.train
  
  m <- length(y.train)
  
  if(length(x.train)==0) {
    p <- 0
    n <- 1
    
    X.train <- matrix(nrow=m, ncol=1, dimnames=list(NULL, 't'))
  } else {
    p <- ncol(x.train)
    
    if(length(x.test)>0) n <- nrow(x.test)
    
    X.train <- matrix(nrow=m, ncol=p+1)
    
    if(length(dimnames(x.train)[[2]])>0)
      dimnames(X.train)[[2]] <- c('t', dimnames(x.train)[[2]])
    else dimnames(X.train)[[2]] <- c('t', paste0('x', 1:p))
  }
  
  k <- 1
  
  for(i in 1:N) for(j in 1:K) if(events[j] <= train.times[i]) {
    if(p==0) X.train[k, ] <- c(events[j])
    else X.train[k, ] <- c(events[j], x.train[i, ])
    
    k <- k+1
  }
  
  ## get tx.test and y.test
  
  ## the real survial times (without censored cases) of test data will be compared to the unique event times of training 
  ## data to get 0 (test event did not occur) or 1 (test event occurred)
  
  if(p==0 | length(x.test)>0) {
    X.test <- matrix(nrow=K*n, ncol=p+1, dimnames=dimnames(X.train))
    y.test <- integer(K*n)
    for(i in 1:n) for(j in 1:K) {
      if(p==0) {
        X.test[j, ] <- c(events[j])
        if (events[j] >= test.times[i]) {
          y.test[j] <- 1
        }  
      }
      else {
        X.test[(i-1)*K+j, ] <- c(events[j], x.test[i, ])
        if (events[j] >= test.times[i]) {
          y.test[(i-1)*K+j] <- 1
        }  
      }
    }
  }
  else X.test <- matrix(nrow=0, ncol=0)*0
  
  return(list(y.train=y.train, y.test=y.test, tx.train=X.train, tx.test=X.test, times=events, K=K))
}

  
  
  
  
  
  


