# This set of code is to minimize information criterion
# BIC: lambda = log(n)
# AIC: lambda = 2


library(mvtnorm)
library(glmnet)


minBIC = function(X, y, beta, interVal, inter, tau3, lambda, preweights, thresh, meth, fam){
	
	#
	#       X: design matrix
	#       y: vector of response
	#    beta: initial guess of coefficients
	#interVal: initial value for intercept
	#   inter: logical, to include intercept or not
	#    tau3: a small constant used to approximate 0-norm
	#  lambda: penalty coefficient (could be a seq) (lambda*sig2)
	#    meth: take values on {"IC", "Reg"}
							# IC: information criterion 
							#Reg: 0-norm regression
	#preweights: weights to make 0-norm approx better
	
	# Return: coefficient estimates and error variances.
	# Moreover: returns BIC values for continuous data and binary.
	# Return NA when meth is specified as not IC. 
	
	#
	# Getting dimensions
	n = length(y)
	p = dim(X)[2]
	
	#
	# Initial values
	sig2 = ifelse(fam=="gaussian" & meth=="IC", t(y-cbind(1,X)%*%c(interVal,beta))%*%(y-cbind(1,X)%*%c(interVal,beta))/n,1)
		
	#
	# Saving devices
	BETA = matrix(-99, nrow = 0, ncol = p+1)
	BETA_inter = c()
	SIGMA2 = c()
	BIC = c()
	
	# 
	# Startting While loop
	currThresh = 1
	iter = 0 
	#print("Starting while loop!")
	while(currThresh > thresh){
		
		beta0 = beta
		sig20 = sig2  

        # Update Beta via EM				
		innerIter = 0
		currthresh1 = 1
		while(currthresh1 > thresh){
			beta00<-beta			
            wLasso = preweights*1/(abs(beta00) + tau3) 
			lassoP = as.numeric(sum(wLasso)/(p)*(lambda*sig20/log(1+1/tau3)))/(2*n)   ###*******
			fit= glmnet(X,y,family=fam,alpha=1,lambda=lassoP,penalty.factor=wLasso,standardize=FALSE, intercept=inter)
			# beta2 = fit$beta
			beta = fit$beta
			currthresh1 = max(abs(beta-beta00))
			# print(paste("inner:",currthresh1))
			# print(beta)
		}

		coeff = coef(fit)
		# Update sigma^2 closed form
		sig2 = ifelse(fam == "gaussian" & meth=="IC",as.numeric(t(y-cbind(1,X)%*%coeff)%*%(y-cbind(1,X)%*%coeff)/n), 1)
					
		currThresh1 = max(abs(beta-beta0))
		currThresh2 = abs(sig20 - sig2)
        currThresh  = max(c(currThresh1,currThresh2))
        # print(paste("outer:",currThresh))        
                     
		BETA = rbind(BETA,as.vector(coeff))
		SIGMA2 = c(SIGMA2, sig2)
		intercept = coef(fit)[1]
		
	}
	
	BIC_cont = as.numeric(n*log(2*pi*sig2) + t(y-cbind(1,X)%*%coeff)%*%(y-cbind(1,X)%*%coeff)/sig2 + lambda*sum(coeff!=0))
	BIC_binary = (-2) * (sum(y*(cbind(1,X)%*%coeff)-log(1+exp(cbind(1,X)%*%coeff)))) + log(n)*sum(coeff!=0)
	
	VAL = NA	
    if(fam == "gaussian" & meth =="IC")  VAL = BIC_cont
    if(fam == "binomial") VAL = BIC_binary
    
	return(list("BETA"=BETA, "SIG2"=SIGMA2, "VAL"=VAL))	
}


