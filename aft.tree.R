# The following code is downloaded from the corresponding author's website: http://odin.mdacc.tmc.edu/~vbaladan/Veera_Home_Page/Software.html
# Reference: Vinicius Bonato et al. Bayesian ensemble methods for survival prediction in gene expression data. Bioinformatics. 2011, 27(3), 359-367
# The original code has bugs. Slight modifications have been made to debug. 

aft.tree<-function(X.train,Y.train,X.test,n.iter,sigma.prior.bart=10,ntree,init.alpha=NULL,burn.in=1,thinning=1){

library(BayesTree)
library(msm)
library(rms)  # contains the new psm (parametric survival model). This library was not loaded in the original script
library(MCMCglmm) # contains rtnorm(). This library was not loaded in the original script

# load data
X.train<-as.matrix(X.train)
X.test<-as.matrix(X.test)
censored.id<-which(Y.train$Status==0)
n.samples<-dim(X.train)[1]
n.samples.test<-dim(X.test)[1]
n.censored<-length(censored.id)
ndpost<-100    # The number of posterior draws returned after burn in
nskip<-50      # Number of MCMC iterations to be treated as burn in

# define containers for alpha, omega, and sigma
sigma.error<-matrix(NA,n.iter,1)
alpha.posterior<-matrix(0,n.iter+1,1)
omega.train<-matrix(NA,n.iter,n.samples)

# impute omega for censored values
omega<-matrix(NA,n.iter,n.censored)
omega.test<-matrix(NA,n.iter,n.samples.test)

# log time
log.time<-lower.bound<-as.vector(log(Y.train$Survival))   # survival time

rel.freq.aft<-matrix(NA,n.iter,dim(X.train)[2])

# MCMC
# n.iter is the number of MCMC iterations
for (aft in 1:n.iter){
	bart1<-bart(X.train,log.time-alpha.posterior[aft],x.test=X.test,sigest=sigma.prior.bart,ntree=ntree,ndpost=ndpost,nskip=nskip,verbose=FALSE)
	sigma.error[aft]<-bart1$sigma[ndpost]
	omega.train[aft,]<-bart1$yhat.train.mean
	omega.test[aft,]<-bart1$yhat.test.mean
	alpha.posterior[aft+1]<-mean(abs(log.time-omega.train[aft,]))
	
	# sample log time for censored cases from the truncated normal distribution
	for(j in 1:length(censored.id)){
		N.censored<-rtnorm(1,mean=alpha.posterior[aft]+omega.train[aft,censored.id[j]],sd=sigma.error[aft],lower=lower.bound[censored.id[j]], upper = Inf)
		log.time[censored.id[j]]<-N.censored
		omega[aft,j]<-N.censored
	}

	varcount<-apply(bart1$varcount,2,sum)  # column-wise sum
	                                       # varcount: a matrix with (ndpost/keepevery) rows and nrow(x.train) columns. 
	                                       # Each row is for a draw. For each variable (corresponding to the columns), 
	                           # the total count of the number of times that variable is used in a tree decision rule (over all trees) is given.
	rel.freq.aft[aft,]<-varcount/(sum(varcount))
	print(sprintf("DRAW_AFT.TREE_%i",aft))
}

#post-processing
draws<-seq(burn.in,n.iter,thinning)
omega.train[,censored.id]<-omega
omega.train<-as.data.frame(omega.train)
names(omega.train)<-row.names(X.train)
omega.test<-as.data.frame(omega.test)
names(omega.test)<-row.names(X.test)
rel.freq.aft<-as.data.frame(rel.freq.aft)
colnames(rel.freq.aft)<-colnames(X.train)

#predictions
Y.hat.train.posterior<-alpha.posterior[draws]+omega.train[draws,]
Y.hat.test.posterior<-alpha.posterior[draws]+omega.test[draws,]
Y.hat.train.stats<-apply(Y.hat.train.posterior,2,quantile,probs=c(0.05,0.5,0.95))
Y.hat.test.stats<-apply(Y.hat.test.posterior,2,quantile,probs=c(0.05,0.5,0.95))
Y.hat.train.stats<-as.data.frame(t(Y.hat.train.stats))
Y.hat.test.stats<-as.data.frame(t(Y.hat.test.stats))

list(predictions=list(yhat.train=Y.hat.train.posterior,yhat.test=Y.hat.test.posterior,yhat.train.stats=Y.hat.train.stats, 
                      yhat.test.stats=Y.hat.test.stats),parameters=list(omega.train=omega.train[draws,],omega.test=omega.test[draws,],
                      alpha=alpha.posterior[draws],sigma.error=sigma.error[draws]),variable.use=data.frame(rel.freq.aft[draws,]))
}

