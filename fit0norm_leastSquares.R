# This function is to optimize the following objective function
# ||y-X%*%beta||_2^2 + lambda * ||beta||_0
# (0-norm penalized least squares problem)


fit0norm_leastSquares=function(X,y,start, fitInter,lambda,tau_seq){
	
	# X: design matrix (no intercept column)
	# y: response variable
	# start: initial value for beta
	# inter: logic, to fit intercept or not
	# lambda: penalty parameter
	#tau_seq: sequence of tau_vals
	
	n=dim(X)[1]
	p=dim(X)[2]
	
	X = t((t(X)-apply(X,2,mean))/apply(X,2,sd))
	
	weight_est = coef(lm(y~X-1))
	weight = log(1+1/tau_seq[1])/log(1+abs(weight_est)/tau_seq[1])
	run = minBIC(X,y,beta=start,interVal=mean(y), inter= fitInter,tau3=tau_seq[1],lambda=lambda,preweights=weight,thresh=1e-6, meth="reg", fam="gaussian")
	finish = run$BETA[dim(run$BETA)[1],]

	if(length(tau_seq) >= 2){
		for(i in 2:length(tau_seq)){
			
			start = finish[-1]	
			weight_est = finish[-1]
			weight_est[weight_est==0] = 1e-5
			weight = log(1+1/tau_seq[i])/log(1+abs(weight_est)/tau_seq[i])
			
			run = minBIC(X,y,beta=start,interVal=finish[1], inter= fitInter,tau3=tau_seq[i],lambda=lambda,preweights=weight,thresh=1e-6, meth="reg", fam="gaussian")
			finish = run$BETA[dim(run$BETA)[1],]
		
		}  # end of for loop
	}  #end of if condition	

	Val = t(y-cbind(1,X)%*%finish)%*%(y-cbind(1,X)%*%finish) + lambda*sum(finish!=0)
	if(fitInter){coef = finish
		}else{coef = finish[-1]}
	
	return(list("Val"=Val, "beta"=coef))
}

