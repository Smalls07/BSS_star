# This function is to optimize the BIC for data sets with continuous response


BIC_linear=function(X,y,beta0,Intercept,tau_seq){
	
	# X: design matrix (no intercept column)
	# y: response variable
	# start: initial value for beta
	# inter: logic, to fit intercept or not
	#tau_seq: sequence of tau_vals
	
	n=dim(X)[1]
	p=dim(X)[2]
	
	lambda = log(n)
	
	X = t((t(X)-apply(X,2,mean))/apply(X,2,sd))
	
	weight_est = coef(lm(y~X-1))
	weight = log(1+1/tau_seq[1])/log(1+abs(weight_est)/tau_seq[1])
	run = minBIC(X,y,beta=beta0,interVal=mean(y), inter= Intercept,tau3=tau_seq[1],lambda=lambda,preweights=weight,thresh=1e-6, meth="IC", fam="gaussian")
	finish = run$BETA[dim(run$BETA)[1],]

	if(length(tau_seq) >= 2){
		for(i in 2:length(tau_seq)){
			
			start = finish[-1]	
			weight_est = finish[-1]
			weight_est[weight_est==0] = 1e-5
			weight = log(1+1/tau_seq[i])/log(1+abs(weight_est)/tau_seq[i])
			
			run = minBIC(X,y,beta=start,interVal=finish[1], inter= Intercept,tau3=tau_seq[i],lambda=lambda,preweights=weight,thresh=1e-6, meth="IC", fam="gaussian")
			finish = run$BETA[dim(run$BETA)[1],]
		
		}  # end of for loop
	}  #end of if condition	

	Val = run$VAL
	if(Intercept){coef = finish
		}else{coef = finish[-1]}
	
	return(list("Val"=Val, "beta"=coef))
}
