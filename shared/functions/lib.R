source_url("https://raw.githubusercontent.com/PerceptionAndCognitionLab/ctx-indiff/public/shared/functions/extendAssign.R")


myF=function(dat)
{
	K=table(dat$sub,dat$cond)
	N=sum(K)
	I=dim(K)[1]
	J=dim(K)[2]
	kstar=K[,1]*K[,2]/(K[,1]+K[,2])
	condMean=tapply(dat$y,list(dat$sub,dat$cond),mean)
	d=condMean[,2]-condMean[,1]
	dstar=d*sqrt(kstar)
	s1=var(dstar)
	ss=sum((K-1)*tapply(dat$y,list(dat$sub,dat$cond),var))
	s2=ss/(N-I*J)  
	return(s1/s2)
}
	
dsinvgamma=function(s,shape,scale) 2*s*dinvgamma(s*s,shape,scale)	


dphnorm=function(x,mu,sigma,log=F){
	a=dnorm(x,mu,sigma,log=T)-pnorm(mu/sigma,log=T)
	if (log) {return(a)} else {return(exp(a))}
	}

rphnorm=function(N,mean,sd) rtnorm(N,mean,sd,lower=0,upper=Inf)

modOne=function(dat,prior,myZero=-1,M=1000)
{
	y=dat$y
	sub=dat$sub
	cond=dat$cond
	if (sum(cond==0)>0) cond[dat$cond==0]=2
	K=table(sub,cond)
	g(I,J)%=%dim(K)
	N=sum(K)

	condMean=tapply(y,list(sub,cond),mean)
	effect=condMean[,2]-condMean[,1]
	meanEffect=mean(effect)
	if (myZero== -1) myZero=meanEffect
	
	g(a,b,muAlpha.m,muAlpha.s2,s2Alpha.a,s2Alpha.b,theta0.mu,theta0.s2)%=%prior
	s2=1:M
	theta0=1:M
	alpha=matrix(nrow=M,ncol=I)
	muAlpha=1:M
	s2Alpha=1:M
	postLNDens=1:M
	
	y.2=y[cond==2]
    sub.2=sub[cond==2]
    SK2=sum(K[,2])
    Ki=apply(K,1,sum)
    
    theta0[1]=meanEffect
	alpha[1,]=condMean[,1]
	s2[1]=.3^2
	muAlpha[1]=mean(condMean[,1])
	s2Alpha[1]=.1^2

	priorLDens=dnorm(myZero,theta0.m,sqrt(theta0.s2),log=T)


	
	for (m in 2:M){
		#alpha
    		vAlpha=1/(Ki/s2[m-1]+1/s2Alpha[m-1])
    		cAlpha=tapply(y-(cond-1)*theta0[m-1],sub,sum)/s2[m-1]+muAlpha[m-1]/s2Alpha[m-1]
	    alpha[m,]=rnorm(I,vAlpha*cAlpha,sqrt(vAlpha))
		#theta0
		vTheta0=1/(SK2/s2[m-1]+1/theta0.s2)
		cTheta0=sum(y.2-alpha[m,sub.2])/s2[m-1]+theta0.m/theta0.s2
		theta0[m]=rnorm(1,vTheta0*cTheta0,sqrt(vTheta0))
		postLNDens[m]=dnorm(myZero,vTheta0*cTheta0,sqrt(vTheta0),log=T)-priorLDens
		#s2
		SSE=sum((y-(alpha[m,sub]+theta0[m]*(cond-1)))^2)
		s2scale=b+SSE/2
		s2[m]=rinvgamma(1,shape=a+N/2,scale=s2scale)
		vMuAlpha=1/((I/s2Alpha[m-1])+1/muAlpha.s2)
		cMuAlpha=(sum(alpha[m,])/s2Alpha[m-1])+(muAlpha.m/muAlpha.s2)
		muAlpha[m]=rnorm(1,vMuAlpha*cMuAlpha,sqrt(vMuAlpha))
		s2scale=s2Alpha.b+sum((alpha[m,]-muAlpha[m])^2)/2
		s2Alpha[m]=rinvgamma(1,shape=s2Alpha.a+I/2,scale=s2scale)}

	out=list(
	's2'=s2,
	'alpha'=alpha,
	'muAlpha'=muAlpha,
	's2Alpha'=s2Alpha,
	'theta0'=theta0,
	'postLNDens'=postLNDens,
	'priorLDens'=priorLDens)	
}


pinvgamma=function(x,a,b,lower.tail=T) pgamma(1/x,a,rate=b,lower.tail=!lower.tail)


modGenLim=function(dat,prior,myZero=-1,M=1000,offset=300,lowerS2Theta=.020^2,upperS2Theta=.1^2){
	y=dat$y
	sub=dat$sub
	cond=dat$cond
	if (sum(cond==0)>0) cond[dat$cond==0]=2
	K=table(sub,cond)
	g(I,J)%=%dim(K)
	N=sum(K)

	condMean=tapply(y,list(sub,cond),mean)
	effect=condMean[,2]-condMean[,1]
	meanEffect=mean(effect)
	if (myZero== -1) myZero=meanEffect
	
	g(a,b,muAlpha.m,muAlpha.s2,s2Alpha.a,s2Alpha.b,muTheta.m,muTheta.s2,s2Theta.a,s2Theta.b)%=%prior

	priorVal=function(w){
		V=matrix(ncol=I,nrow=I,rep(muTheta.s2))+w*diag(I)
		return(exp(dmvnorm(rep(myZero,I),rep(muTheta.m,I),V,log=T)
		+log(dinvgamma(w,s2Theta.a,s2Theta.b))-log(pinvgamma(lowerS2Theta,s2Theta.a,s2Theta.b,lower.tail=F)) -offset))}
	priorValVec=Vectorize(priorVal)
	priorDens=integrate(priorValVec,lower=lowerS2Theta,upper=1)$value
	priorLDens=log(priorDens)


	s2=1:M
	alpha=matrix(nrow=M,ncol=I)
	theta=matrix(nrow=M,ncol=I)
	muAlpha=s2Alpha=1:M
	muTheta=s2Theta=1:M
	postLNDens=1:M
	
	y.2=y[cond==2]
    sub.2=sub[cond==2]
    Ki=apply(K,1,sum)
    
    alpha[1,]=condMean[,1]
	theta[1,]=effect
	s2[1]=.3^2
	muAlpha[1]=mean(alpha[1,])
	s2Alpha[1]=var(alpha[1,])
	muTheta[1]=mean(theta[1,])
	s2Theta[1]=var(theta[1,])
	countPos=0
 	for (m in 2:M){
		#alpha
    		vAlpha=1/(Ki/s2[m-1]+1/s2Alpha[m-1])
    		cAlpha=tapply(y-(cond-1)*theta[m-1,sub],sub,sum)/s2[m-1]+muAlpha[m-1]/s2Alpha[m-1]
    		alpha[m,]=rnorm(I,vAlpha*cAlpha,sqrt(vAlpha))
		#theta
		vTheta=1/(K[,2]/s2[m-1]+1/s2Theta[m-1])
		cTheta=tapply(y.2-alpha[m,sub.2],sub.2,sum)/s2[m-1]+muTheta[m-1]/s2Theta[m-1]
		theta[m,]=rnorm(I,vTheta*cTheta,sqrt(vTheta))
		if (sum(theta[m,]>0)==I) countPos=countPos+1
		postLNDens[m]=sum(dnorm(myZero,vTheta*cTheta,sqrt(vTheta),log=T))-offset-priorLDens
		#s2
		SSE=sum((y-(alpha[m,sub]+theta[m,sub]*(cond-1)))^2)
		s2scale=b+SSE/2
		s2[m]=rinvgamma(1,shape=a+N/2,scale=s2scale)
		#alpha hierarchy
		vMuAlpha=1/((I/s2Alpha[m-1])+1/muAlpha.s2)
		cMuAlpha=(sum(alpha[m,])/s2Alpha[m-1])+(muAlpha.m/muAlpha.s2)
		muAlpha[m]=rnorm(1,vMuAlpha*cMuAlpha,sqrt(vMuAlpha))
		s2scale=s2Alpha.b+sum((alpha[m,]-muAlpha[m])^2)/2
		s2Alpha[m]=rinvgamma(1,shape=s2Alpha.a+I/2,scale=s2scale)
		#theta hierarchyTheta
		vMuTheta=1/((I/s2Theta[m-1])+1/muTheta.s2)
		cMuTheta=(sum(theta[m,])/s2Theta[m-1])+(muTheta.m/muTheta.s2)
		muTheta[m]=rphnorm(1,vMuTheta*cMuTheta,sqrt(vMuTheta))
		s2scale=s2Theta.b+sum((theta[m,]-muTheta[m])^2)/2
		repeat{
			s2Theta[m]=rinvgamma(1,shape=s2Theta.a+I/2,scale=s2scale)
			if(s2Theta[m]>lowerS2Theta) {break}}
		}


	out=list(
	's2'=s2,
	'alpha'=alpha,
	'muAlpha'=muAlpha,
	's2Alpha'=s2Alpha,
	'theta'=theta,
	'muTheta'=muTheta,
	's2Theta'=s2Theta,	
	'postLNDens'=postLNDens,
	'priorLDens'=priorLDens,
	'countPos'=countPos)	
}


modTruncLim=function(dat,prior,myZero=-1,M=1000,muThetaSD=.015,s2ThetaSD=.001, lowerS2Theta = .010^2){
  y=dat$y
  sub=dat$sub
  cond=dat$cond
  if (sum(cond==0)>0) cond[dat$cond==0]=2
  K=table(sub,cond)
  g(I,J)%=%dim(K)
  N=sum(K)
  
  condMean=tapply(y,list(sub,cond),mean)
  effect=condMean[,2]-condMean[,1]
  meanEffect=mean(effect)
  if (myZero== -1) myZero=meanEffect
  
  g(a,b,muAlpha.m,muAlpha.s2,s2Alpha.a,s2Alpha.b,muTheta.m,muTheta.s2,s2Theta.a,s2Theta.b)%=%prior
  s2=1:M
  alpha=matrix(nrow=M,ncol=I)
  theta=matrix(nrow=M,ncol=I)
  muAlpha=s2Alpha=1:M
  muTheta=s2Theta=1:M
  
  y.2=y[cond==2]
  sub.2=sub[cond==2]
  Ki=apply(K,1,sum)
  
  alpha[1,]=condMean[,1]
  theta[1,]=effect
  s2[1]=.3^2
  muAlpha[1]=mean(alpha[1,])
  s2Alpha[1]=var(alpha[1,])
  muTheta[1]=mean(theta[1,])
  s2Theta[1]=var(theta[1,])
  
  #MetHasting
  counter=c(0,0)
  
  for (m in 2:M){
    #alphA
    vAlpha=1/(Ki/s2[m-1]+1/s2Alpha[m-1])
    cAlpha=tapply(y-(cond-1)*theta[m-1,sub],sub,sum)/s2[m-1]+muAlpha[m-1]/s2Alpha[m-1]
    alpha[m,]=rnorm(I,vAlpha*cAlpha,sqrt(vAlpha))
    #theta
    vTheta=1/(K[,2]/s2[m-1]+1/s2Theta[m-1])
    cTheta=tapply(y.2-alpha[m,sub.2],sub.2,sum)/s2[m-1]+muTheta[m-1]/s2Theta[m-1]
    theta[m,]=rtnorm(I,vTheta*cTheta,sqrt(vTheta),0,Inf)
    #postDens[m]=exp(sum(dphnorm(myZero,vTheta*cTheta,sqrt(vTheta),log=T)))
    #s2
    SSE=sum((y-(alpha[m,sub]+theta[m,sub]*(cond-1)))^2)
    s2scale=b+SSE/2
    s2[m]=rinvgamma(1,shape=a+N/2,scale=s2scale)
    #alpha hierarchy
    vMuAlpha=1/((I/s2Alpha[m-1])+1/muAlpha.s2)
    cMuAlpha=(sum(alpha[m,])/s2Alpha[m-1])+(muAlpha.m/muAlpha.s2)
    muAlpha[m]=rnorm(1,vMuAlpha*cMuAlpha,sqrt(vMuAlpha))
    s2scale=s2Alpha.b+sum((alpha[m,]-muAlpha[m])^2)/2
    s2Alpha[m]=rinvgamma(1,shape=s2Alpha.a+I/2,scale=s2scale)
    #theta hierarchyTheta
    muTheta[m]=muTheta[m-1]
    cand=muTheta[m]+rnorm(1,0,muThetaSD)
    if (cand>0){		
      fullCurM=sum(dtnorm(theta[m,],muTheta[m],sqrt(s2Theta[m-1]),0,Inf,log=T)) +
        dphnorm(muTheta[m],muTheta.m,sqrt(muTheta.s2),log=T)
      fullCandM=sum(dtnorm(theta[m,],cand,sqrt(s2Theta[m-1]),0,Inf,log=T)) + 
        dphnorm(cand,muTheta.m,sqrt(muTheta.s2),log=T)
      prob=min(exp(fullCandM-fullCurM),1)
      if (rbinom(1,1,prob)){muTheta[m]=cand;counter[1]=counter[1]+1}}	
    s2Theta[m]=s2Theta[m-1]
    cand=s2Theta[m]+rnorm(1,0,s2ThetaSD)
    if (cand>lowerS2Theta){
      fullCurS=sum(dtnorm(theta[m,],muTheta[m],sqrt(s2Theta[m]),0,Inf,log=T)) +
        log(dinvgamma(s2Theta[m],s2Theta.a,s2Theta.b))
      fullCandS=sum(dtnorm(theta[m,],muTheta[m],sqrt(cand),0,Inf,log=T)) + 
        log(dinvgamma(cand,s2Theta.a,s2Theta.b))
      prob=min(exp(fullCandS-fullCurS),1)
      #p=ifelse(is.na(prob),0,prob)
      if (rbinom(1,1,prob)){s2Theta[m]=cand;counter[2]=counter[2]+1}}}
  
  out=list(
    's2'=s2,
    'alpha'=alpha,
    'muAlpha'=muAlpha,
    's2Alpha'=s2Alpha,
    'theta'=theta,
    'muTheta'=muTheta,
    's2Theta'=s2Theta,	
    'counter'=counter)	
}


#####################

readStroopOb=function()
{	
filename <- curl("https://raw.githubusercontent.com/PerceptionCognitionLab/data0/master/contexteffects/FlankerStroopSimon/LEF_stroop.csv")
stroop <- read.csv2(filename, header=TRUE, dec=".")
stroop$cond <- as.numeric(stroop$congruency)  #congruent -> 1, incongruent -> 2, neutral -> 3
ntrial <- length(stroop[stroop$id == stroop$id[1], 1])
nsub <- length(unique(stroop$id))
stroop$trial <- rep(1:ntrial, nsub)
stroop$sub <- rep(1:nsub, each = ntrial)
stroop$rt <- stroop$RT/1000
stroop <- stroop[stroop$rt > .2 & stroop$rt < 2, ]
stroop <- subset(stroop, accuracy == 1 & cond != 3)
dat=stroop
dat$y=dat$rt
return(dat)
}

readStroopPratte1=function()
{
filename <- curl("https://raw.githubusercontent.com/PerceptionCognitionLab/data0/master/contexteffects/StroopSimonAPP2010/allsi2.dat")
clnames <- c('exp', 'sub', 'blk', 'trial', 'color', 'distract', 'cond', 'resp', 'acc', 'rt', 'errorTotal')
dat <- read.table(filename)
colnames(dat) <- clnames
dat <- dat[dat$rt > .2 & dat$rt < 2, ]
dat <- subset(dat, acc == 1 & cond != 2 & exp == 1)
dat <- dat[!(dat$trial %in% 1:5), ]
dat$y=dat$rt
dat$cond[dat$cond==0]=2
return(dat)
}