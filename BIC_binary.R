# This function is to optimize the BIC for data sets with binary response


BIC_binary=function(X,y,start,fitInter,tau_seq){
	
	# X: design matrix (no intercept column)
	# y: response variable
	# start: initial value for beta
	# fitInter: logic, to fit intercept (T) or not (F)
	#tau_seq: sequence of tau_vals
	
	n=dim(X)[1]
	p=dim(X)[2]
	
	X = t((t(X)-apply(X,2,mean))/apply(X,2,sd))
	

	weight_est = coef(glm(y~X-1, family="binomial"))
	weight = log(1+1/tau_seq[1])/log(1+abs(weight_est)/tau_seq[1])
	
	run = minBIC(X,y,beta=start,interVal=1, inter= fitInter,tau3=tau_seq[1],lambda=log(n),preweights=weight,thresh=1e-6, meth="IC", fam="binomial")
	finish = run$BETA[dim(run$BETA)[1],]
	
	if(length(tau_seq)>1){
		for(i in 2:length(tau_seq)){
			# i=i+1
			# print(i)
			start = finish[-1]
			weight_est = finish[-1]
			weight_est[weight_est==0] = 1e-5
			weight = log(1+1/tau_seq[i])/log(1+abs(weight_est)/tau_seq[i])
		
			run = minBIC(X,y,beta=start,interVal=finish[1], inter= fitInter,tau3=tau_seq[i],lambda=log(n),preweights=weight,thresh=1e-6, meth="IC", fam="binomial")
			finish = run$BETA[dim(run$BETA)[1],]
		}
		finish		
	}	
	

	Val = run$VAL
	if(fitInter){coef = finish
		}else{coef = finish[-1]}
	
	return(list("Val"=Val, "beta"=coef))
}
