##########R functions for estimation and smoothing state-space harvest-based models
####Citation:Fukasawa, K., Osada, Y., Iijima, H. (in review) Is harvest size a valid indirect measure of abundance for evaluating the population size of game animals using harvest-based estimation?
####2020/07/17

####Auxiliary functions
## function to calculate capture prob.
#[argument]
	#c:capture efficiency, positive-valued scalar
	#E:vector of capture effort, E>=0
#[value]
	#vector of capture probability

capt <- function(c,E){ 1-exp(-c*E) }

##function to calculate log(sum(exp(lp))) avoiding digit loss
#[argument]
	#lp:vector of log probability (lp<=0)
#[value]
	#log(sum(exp(lp)))

sum_lp <- function(lp){ max_lp <- max(lp); log(sum(exp(lp-max_lp)))+max_lp }

## function to normalize exp(lp) 
#[argument]
	#lp:vector of log probability (lp<=0)
#[value]
	#logarithm of normalized exp(lp)

std_lp <- function(lp){ p <- exp(lp-max(lp)); log(p/sum(p)) }


####function to calculate log marginal likelihood
#[argument]
	#pars:vector of parameters, c(log(1+intrinsic growth rate), log(capture efficiency))
	#C:vector of harvest size, integer, C>=0
	#E:vector of capture effort, E>=0
	#Smax:maximal population size after population growth, S[t]
	#scale: multiplier for log marginal likelihood, scalar
#[value]
	#log marginal likelihood*scale

loglf.hbm <- function(pars, C, E, Smax=10000, scale=-1) {
	ny <- length(C)
	r  <- exp(pars[1])
	c  <- exp(pars[2])
	L  <- rep(NA,ny)
	
	# Transition Probability
	TP <- sapply(0:Smax, function(k) dpois(k,r*0:Smax,log=TRUE)) 
	
	# t = 1
	pred <- rep(-log(Smax+1),Smax+1)
	filt <- pred + dbinom(C[1],0:Smax,capt(c,E[1]),log=TRUE)
	L[1] <- sum_lp(filt)
	filt <- c(filt[-(1:C[1])],rep(-Inf,C[1]))
	filt <- std_lp(filt)
	
	# t > 1
	for(t in 2:ny) {
		for(k in 0:Smax) pred[k+1] <- sum_lp(TP[,k+1]+filt)	#prediction
		pred <- std_lp(pred)
		filt <- pred + dbinom(C[t],0:Smax,capt(c,E[t]),log=TRUE)	#filtering
		L[t] <- sum_lp(filt)
		filt <- c(filt[-(1:C[t])],rep(-Inf,C[t]))	#
		filt <- std_lp(filt)
	}
	res <- sum(L) * scale
	cat("pars: ",pars,"; val:",res,"\n")
	return(res)
}

####function to obtain marginal maximum likelihood estimate
#[argument]
	#C: vector of harvest size, integer, C>=0
	#E: vector of capture effort, E>=0
	#initpar: vector of initial parameters, c(log(1+intrinsic growth rate), log(capture efficiency)), default c(0.2, -6)
	#Smax: maximal population size after population growth, S[t], default 10000
	#...: additional parameters passed to nlm()
#[value]
	#list of length 5
		#$summary: summary table of estimate
		#$nlm: output from nlm()
		#$C, $E, $Smax: input data

fit.hbm<-function(C,E,initpar=c(0.2,-6),Smax=10000,...){
	res<-nlm(loglf.hbm,initpar,C=C,E=E,Smax=Smax,hessian=T,...)
	se<-sqrt(diag(solve(res$hessian)))
	ci2.5<-res$estimate+qnorm(0.025,0,1)*se
	ci97.5<-res$estimate+qnorm(0.975,0,1)*se
	summ<-cbind(mle=res$estimate,se=se,ci2.5=ci2.5,ci97.5=ci97.5)
	rownames(summ)<-c("ln(r)","ln(c)")
	out<-list(summary=summ,nlm=res,C=C,E=E,Smax=Smax)
	return(out)
}


####function to obtain smoothed distribution of population size
#[argument]
	#hbmfit:output from fit.hbm()
#[value]
	#list of length 2
		#$muN: vector of mean population size, Nt
		#$smooth: smoothed distribution, list of length 2
			#$N0: smoothed distribution of initial population size, N0
			#$St: smoothed distribution of population size after reproduction, St

smooth.hbm <- function(hbmfit) {
	pars<-hbmfit$nlm$estimate
	C<-hbmfit$C
	E<-hbmfit$E
	Smax<-hbmfit$Smax

	ny <- length(C)
	r  <- exp(pars[1])
	c  <- exp(pars[2])
	L  <- rep(NA,ny)
	
	# Transition Probability
	TP <- sapply(0:Smax, function(k) dpois(k,r*0:Smax,log=TRUE)) 
	
	# t = 1
	pred<-matrix(NA,nrow=Smax+1,ncol=ny)
	filt<-matrix(NA,nrow=Smax+1,ncol=ny)
	pred[,1] <- rep(-log(Smax+1),Smax+1)
	filt[,1] <- pred[,1] + dbinom(C[1],0:Smax,capt(c,E[1]),log=TRUE)
	L[1] <- sum_lp(filt[,1])
	filt[,1] <- std_lp(filt[,1])
	temp <- c(filt[-(1:C[1]),1],rep(-Inf,C[1]))
	
	# t > 1
	for(t in 2:ny) {
		pred[,t]<-npforeach(k=0:Smax,.combine="c")({
				sum_lp(TP[,k+1]+temp)
		})
		#for(k in 0:Smax) pred[k+1,t] <- sum_lp(TP[,k+1]+filt[,t-1])
		pred[,t] <- std_lp(pred[,t])
		filt[,t] <- pred[,t] + dbinom(C[t],0:Smax,capt(c,E[t]),log=TRUE)
		L[t] <- sum_lp(filt[,t])
		filt[,t] <- std_lp(filt[,t])
		temp <- c(filt[-(1:C[t]),t],rep(-Inf,C[t]))
	}

	cat("smooth")
	smooth<-matrix(NA,nrow=Smax+1,ncol=ny)
	smooth[,ny]<-filt[,ny]
	for(t in (ny-1):1){
		cat("t=",t,"\n")
		for(k in 0:Smax) {
			temp<-TP[k+1,]-pred[,t+1]+smooth[,t+1]
			temp<-ifelse(is.nan(temp),-Inf,temp)
			smooth[k+1,t]<-sum_lp(temp)
		}
		smooth[,t]<-ifelse(is.nan(smooth[,t]),-Inf,smooth[,t])
		smooth[,t]<-c(rep(-Inf,C[t]),smooth[1:(Smax+1-C[t]),t])
		smooth[,t]<-filt[,t]+smooth[,t]
		smooth[,t]<-std_lp(smooth[,t])
	}
	n0<-rep(NA,Smax+1)
	for(k in 0:Smax) {
		temp<-TP[k+1,]-pred[,1]+smooth[,1]
		temp<-ifelse(is.nan(temp),-Inf,temp)
		n0[k+1]<-sum_lp(temp)
	}
	n0<-ifelse(is.nan(n0),-Inf,n0)
	n0<-std_lp(n0)
	
	N0estim<-sum(exp(n0)*0:Smax)
	Ntestim<-apply(exp(smooth)*outer(0:Smax,rep(1,ny)),2,sum)-C
	muN<-c(N0estim,Ntestim[-ny])

	return(list(muN=muN,smooth=list(N0=n0,St=smooth)))
}


