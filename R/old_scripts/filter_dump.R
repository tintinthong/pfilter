#--LIBRARY FOR ALL PARTICLE FILTERS--

#--Load libraries--
source("generate_mod.R")
source("resample_mod.R")
source("diagnostic_mod.R")


#--Generic particle filter(prior proposal)--
# set N.thr=1, if want to resample every time 
# resampling is relatively slow when it is hard-coded(may resort to base resample() function)
# this function is for 1D random walk example
# N- No of particles
# x- vector of states
# y- vector of measurements
# x0- starting value
# sigma- standard deviation of process 
# sigma.meas- standard deviation of measurment
# resample.type - type of resampling (multinomial,stratified,systematic,residual,standard)
particle.filter<-function(N,x,y,x0,sigma=1,sigma.meas=1,resample.type="multinomial",N.thr.per=0.5){
  
  #threshold N
  N.thr<-as.integer(N.thr.per*N)
  
  #set up filter
  x0.pf<-rnorm(N,mean=x0,sd=sigma) #suppose this is the correct sigma
  x.pf <- matrix(NA,ncol=N,nrow=tau)
  w.mat<-matrix(NA,ncol=N,nrow=tau)
  
  #set up diagnostics
  #- adds variables for diagnoistics
  diagnostic.init(tau)
  
  #local function to choose resampling method
  FUN.resample<-function(N,x,w){
    switch(resample.type,
           multinomial = mult.resample(N=N,x=x,w=w) ,
           stratified = strat.resample(N=N,x=x,w=w) ,
           systematic = sys.resample(N=N,x=x,w=w),
           residual= res.resample(N=N,x=x,w=w),
           standard=sample(x=x,size=N,replace=TRUE,prob=w)
    )
  }
  
  # Initialize
  x.pf[1, ] <- x0.pf+rnorm(N,sd=sigma) #project particles from x0 to the first time step
  w.tilde<- dnorm(y[1],mean=x.pf[1,],sd=sigma.meas) #weights of x0 is equal
  w<-w.tilde/sum(w.tilde) 
  w.mat[1,]<-w.tilde/sum(w.tilde) 
  
  
  
  #compute diagnostics
  #- mean, effective sample size, confidence intervals 
  m[1]<-sum(w*x.pf[1,]) 
  MSE.k[1]<-(x[1]-m[1])^2
  N.eff[1]<-1/sum(w^2) 
  var.s[1]<-sum(w*(x.pf[1,]-m[1])^2)
  conf[1,]<-quantile(x.pf[1,],probs=c(0.025,0.975))
  
  
  for (t in 2:(tau)) {
    
    
    #project  state forward
    x.pf[t, ] <-x.pf[t-1,]+rnorm(N,sd=sigma)
    
    #compute weights
    w.tilde <- dnorm(y[t], mean=x.pf[t, ],sd=sigma.meas)*w
    w <- w.tilde/sum(w.tilde) 
    w.mat[t,]<-w.tilde/sum(w.tilde) 
    
    
    #compute diagnostics
    N.eff[t]<-1/sum(w^2)
    m[t]<-sum(w*x.pf[t,])
    MSE.k[t]<-(x[t]-m[t])^2
    var.s[t]<-sum(w*(x.pf[t,]-m[t])^2)
    conf[t,]<-quantile(x.pf[t,],probs=c(0.025,0.975))
    
    
    #resample (bootstrap means resample everytime)
    if(N.eff[t]<N.thr){
      x.pf[t,]<-FUN.resample(N=N,x=x.pf[t,],w=w) 
      w<-rep(1/N,N)
      res.point[t]<-1 #add "yes" for resampling point
    }
    
    
    
  }
  
  #return output
  return(list(x.out=x,y.out=y, x.pf.out=x.pf, m.out=m,MSE.k.out=MSE.k, 
              conf.out=conf,N.eff.out=N.eff,w.out=w,res.point.out=res.point,w.mat.out=w.mat,var.s.out=var.s))
  
}



#--Generic particle filter(prior proposal) Smoothed--
# Resample paths, not particles
# N- No of particles
# x- vector of states
# y- vector of measurements
# x0- starting value
# sigma- standard deviation of process 
# sigma.meas- standard deviation of measurment
# resample.type - type of resampling (multinomial,stratified,systematic,residual,standard)
particle.filter.path<-function(N,x,y,x0,sigma=1,sigma.meas=1,resample.type="multinomial",N.thr.per=0.5){
  
  #compute tau
  tau<-length(x)
  
  #threshold N
  N.thr<-as.integer(N.thr.per*N)
  
  #set up filter
  x0.pf<-rnorm(N,mean=x0,sd=sigma)
  x.pf <- matrix(NA,ncol=N,nrow=tau)
  
  #set up diagnostics
  diagnostic.init(tau)
  
  #local function to choose resampling method
  FUN.resample<-function(N,x,w){
    switch(resample.type,
           multinomial = mult.resample(N=N,x=x,w=w) ,
           stratified = strat.resample(N=N,x=x,w=w) ,
           systematic = sys.resample(N=N,x=x,w=w),
           residual= res.resample(N=N,x=x,w=w),
           standard=sample(x=x,size=N,replace=TRUE,prob=w)
    )
  }
  
  # Initialize
  x.pf[1, ] <- x0.pf+rnorm(N,sd=sigma) #project particles from x0 to the first time step
  w.tilde<- dnorm(y[1],mean=x.pf[1,],sd=sigma.meas) #weights of x0 is equal
  w<-w.tilde/sum(w.tilde) 
  
  #compute diagnostics
  #- mean, effective sample size, confidence intervals 
  m[1]<-sum(w*x.pf[1,]) 
  MSE.k[1]<-(x[1]-m[1])^2
  N.eff[1]<-1/sum(w^2) 
  conf[1,]<-quantile(x.pf[1,],probs=c(0.025,0.975))
  
  
  for (t in 2:(tau)) {
    
    
    #project  state forward
    x.pf[t, ] <-x.pf[t-1,]+rnorm(N,sd=sigma)
    
    #compute weights
    w.tilde <- dnorm(y[t], mean=x.pf[t, ],sd=sigma.meas)*w 
    w <- w.tilde/sum(w.tilde) 
    
    #compute diagnostics
    N.eff[t]<-1/sum(w^2)
    m[t]<-sum(w*x.pf[t,])
    MSE.k[t]<-(x[t]-m[t])^2
    conf[t,]<-quantile(x.pf[t,],probs=c(0.025,0.975))
    
    #resample (bootstrap means resample everytime)
    if(N.eff[t]<N.thr){
      s<-FUN.resample(N=N,x=1:N,w=w) 
      x.pf<-x.pf[,s]
      w<-rep(1/N,N)
      res.point[t]<-1 #add "yes" for resampling point
    }
    
    
  }
  
  #return output
  return(list(x.out=x,y.out=y, x.pf.out=x.pf, m.out=m,MSE.k.out=MSE.k, 
              conf.out=conf,N.eff.out=N.eff,w.out=w,res.point.out=res.point))
  
}


#--Generic particle filter(fixed proposal)--
# this function is for 1D random walk example
# fixed filter chooses normal centred at x0
# N- No of particles
# x- vector of states
# y- vector of measurements
# x0- starting value
# sigma- standard deviation of process 
# sigma.meas- standard deviation of measurment
# resample.type - type of resampling (multinomial,stratified,systematic,residual,standard)

fixed.filter<-function(N,x,y,x0,sigma=1,sigma.meas=1,resample.type="multinomial",N.thr.per=0.5){
  
  #threshold N
  N.thr<-as.integer(N.thr.per*N)
  
  #set up filter
  x0.pf<-rnorm(N,mean=x0,sd=sigma)
  x.pf <- matrix(NA,ncol=N,nrow=tau)
  
  #set up diagnostics
  #-adds variables for diagnoistics
  diagnostic.init(tau)
  
  #local function to choose resampling method
  FUN.resample<-function(N,x,w){
    switch(resample.type,
           multinomial = mult.resample(N=N,x=x,w=w) ,
           stratified = strat.resample(N=N,x=x,w=w) ,
           systematic = sys.resample(N=N,x=x,w=w),
           residual= res.resample(N=N,x=x,w=w),
           standard=sample(x=x,size=N,replace=TRUE,prob=w)
    )
  }
  
  # Initialize
  x.pf[1, ] <-rnorm(N,x0) #project particles from x0 to the first time step
  w.tilde<- dnorm(y[1],mean=x.pf[1,],sd=sigma.meas)*dnorm(x.pf[1,],x0.pf,sd=sigma)/dnorm(x.pf[1,],x0.pf) #weights of x0 is equal
  w<-w.tilde/sum(w.tilde) 
  
  
  
  #compute diagnostics
  #- mean, effective sample size, confidence intervals 
  m[1]<-sum(w*x.pf[1,]) 
  MSE.k[1]<-(x[1]-m[1])^2
  N.eff[1]<-1/sum(w^2) 
  conf[1,]<-quantile(x.pf[1,],probs=c(0.025,0.975))
  
  
  for (t in 2:(tau)) {
    
    #project  state forward
    x.pf[t, ] <-rnorm(N,x0,sd=sigma)
    
    #compute weights
    w.tilde <- (dnorm(y[t], mean=x.pf[t, ],sd=sigma.meas)*dnorm(x.pf[t,],x.pf[t-1,],sd=sigma)*w )/dnorm(x.pf[t,],x0.pf)
    w <- w.tilde/sum(w.tilde) 
    
    #compute diagnostics
    N.eff[t]<-1/sum(w^2)
    m[t]<-sum(w*x.pf[t,])
    MSE.k[t]<-(x[t]-m[t])^2
    conf[t,]<-quantile(x.pf[t,],probs=c(0.025,0.975))
    
    #resample (bootstrap means resample everytime)
    if(N.eff[t]<N.thr){
      x.pf[t,]<-FUN.resample(N=N,x=x.pf[t,],w=w) 
      w<-rep(1/N,N)
      res.point[t]<-1 #add "yes" for resampling point
    }
    
    
    
  }
  
  #return output
  return(list(x.out=x,y.out=y, x.pf.out=x.pf, m.out=m,MSE.k.out=MSE.k, 
              conf.out=conf,N.eff.out=N.eff,w.out=w,res.point.out=res.point))
  
}


#--Generic Paricle Filter (Optimum proposal)--
# this function is for 1D random walk example
# N- No of particles
# x- vector of states
# y- vector of measurements
# x0- starting value
# sigma- standard deviation of process 
# sigma.meas- standard deviation of measurment
# resample.type - type of resampling 

optimum.filter<-function(N,x,y,x0,sigma=1,sigma.meas=1,resample.type="multinomial",N.thr.per=0.5){
  
  #threshold N
  N.thr<-as.integer(N.thr.per*N)
  
  #set up filter
  x0.pf<-rnorm(N,mean=x0,sd=sigma)
  x.pf <-matrix(NA,ncol=N,nrow=tau)
  
  #set up diagnostics
  #-adds variables for diagnoistics
  diagnostic.init(tau)
  
  
  #local function to choose resampling method
  FUN.resample<-function(N,x,w){
    switch(resample.type,
           multinomial = mult.resample(N=N,x=x,w=w) ,
           stratified = strat.resample(N=N,x=x,w=w) ,
           systematic = sys.resample(N=N,x=x,w=w),
           residual= res.resample(N=N,x=x,w=w),
           standard=sample(x=x,size=N,replace=TRUE,prob=w)
    )
  }
  
  
  #initialise
  sigma.prop<-1/sqrt(1/sigma^2+1/sigma.meas^2)
  mu.prop<-sigma.prop^2*(x0.pf/sigma^2+y[1]/sigma.meas^2)
  x.pf[1, ] <- mu.prop+rnorm(N,sd=sigma.prop) 
  w.tilde<- dnorm(y[1],mean=x0.pf,sd=sqrt(sigma^2+sigma.meas^2))
  w<-w.tilde/sum(w.tilde)
  
  
  #compute diagnostics
  #- mean, effective sample size, confidence intervals 
  m[1]<-sum(w*x.pf[1,]) 
  MSE.k[1]<-(x[1]-m[1])^2
  N.eff[1]<-1/sum(w^2) 
  conf[1,]<-quantile(x.pf[1,],probs=c(0.025,0.975))
  
  for (t in 2:(tau)) {
    
    #initialise
    sigma.prop<-1/sqrt(1/sigma^2+1/sigma.meas^2)
    mu.prop<-sigma.prop^2*(x.pf[t-1,]/sigma^2+y[t]/sigma.meas^2)
    
    #project state forward
    x.pf[t, ] <- mu.prop+rnorm(N,sd=sigma.prop) 
    
    #compute weights
    w.tilde<- dnorm(y[t],mean=x.pf[t-1,],sd=sqrt(sigma^2+sigma.meas^2))*w
    w<-w.tilde/sum(w.tilde)
    
    
    #compute diagnostics
    m[t]<-sum(w*x.pf[t,])
    MSE.k[t]<-(x[t]-m[t])^2
    conf[t,]<-quantile(x.pf[t,],probs=c(0.025,0.975))
    N.eff[t]<-1/sum(w^2)
    
    #resample (bootstrap means resample everytime)
    if(N.eff[t]<N.thr){
      x.pf[t,]<-FUN.resample(N=N,x=x.pf[t,],w=w)
      w<-rep(1/N,N)
      res.point[t]<-1 #add "yes" for resampling point
      
    }
    
    
    
  }
  
  #return output
  return(list(x.out=x,y.out=y, x.pf.out=x.pf, m.out=m,MSE.k.out=MSE.k, 
              conf.out=conf, N.eff.out=N.eff,w.out=w,res.point.out=res.point))
  
}

#Exact Particle Filter 
# did not resample on the first step
# N- No of particles
# x- vector of states
# y- vector of measurements
# x0- starting value
# sigma- standard deviation of process 
# sigma.meas- standard deviation of measurment
# resample.type - type of resampling 

exact.filter<-function(N,x,y,x0,sigma=1,sigma.meas=1,resample.type="multinomial"){
  
  
  #set up filter
  x0.pf<-rnorm(N,mean=x0,sd=sigma)
  x.pf <-matrix(NA,ncol=N,nrow=tau)
  
  #set up diagnostics
  #-adds variables for diagnoistics
  diagnostic.init(tau)
  
  #local function to choose resampling method
  FUN.resample<-function(N,x,w){
    switch(resample.type,
           multinomial = mult.resample(N=N,x=x,w=w) ,
           stratified = strat.resample(N=N,x=x,w=w) ,
           systematic = sys.resample(N=N,x=x,w=w),
           residual= res.resample(N=N,x=x,w=w),
           standard=sample(x=x,size=N,replace=TRUE,prob=w)
    )
  }
  
  #initialise
  w.tilde<- dnorm(y[1],mean=x0.pf,sd=sqrt(sigma^2+sigma.meas^2))
  w<-w.tilde/sum(w.tilde)
  x0.pf<-FUN.resample(N=N,x=x0.pf,w=w) #changed t to t-1
  sigma.prop<-1/sqrt(1/sigma^2+1/sigma.meas^2) #1/sqrt(1/sigma^2+1/sigma.meas^2)
  mu.prop<-sigma.prop^2*(x0.pf/sigma^2+y[1]/sigma.meas^2)#sigma.prop^2*(x0.pf/sigma+y[1]/sigma.meas)
  x.pf[1, ] <- mu.prop+rnorm(N,sd=sigma.prop) 
  
  
  
  
  #compute diagnostics
  #- mean, effective sample size, confidence intervals 
  m[1]<-sum(w*x.pf[1,]) #compute mean
  MSE.k[1]<-(x[1]-m[1])^2
  N.eff[1]<-1/sum(w^2) #compute effective sample size
  conf[1,]<-quantile(x.pf[1,],probs=c(0.025,0.975))
  
  for (t in 2:(tau)) {
    
    
    #compute weights
    w.tilde<- dnorm(y[t],mean=x.pf[t-1,],sd=sqrt(sigma^2+sigma.meas^2))#*w
    w<-w.tilde/sum(w.tilde)
    
    #resample
    x.pf[t-1,]<-FUN.resample(N=N,x=x.pf[t-1,],w=w) #changed t to t-1
    #w<-rep(1/N,N)
    
    #initialise
    sigma.prop<-1/sqrt(1/sigma^2+1/sigma.meas^2)#1/sqrt(1/sigma^2+1/sigma.meas^2)
    mu.prop<-sigma.prop^2*(x.pf[t-1,]/sigma^2+y[t]/sigma.meas^2) #sigma.prop^2*(x.pf[t-1,]/sigma^2+y[t]/sigma.meas^2)
    
    #project state forward
    x.pf[t, ] <- mu.prop+rnorm(N,sd=sigma.prop) 
    
    #compute diagnostics
    N.eff[t]<-1/sum(w^2)
    m[t]<-sum(w*x.pf[t,])
    MSE.k[t]<-(x[t]-m[t])^2
    conf[t,]<-quantile(x.pf[t,],probs=c(0.025,0.975))
    
  }
  
  #return output
  return(list(x.out=x,y.out=y, x.pf.out=x.pf, m.out=m,MSE.k.out=MSE.k, conf.out=conf,
              N.eff.out=N.eff,w.out=w,res.point.out=res.point))
  
}


#--Generic particle filter(prior proposal) Population example--
# N- No of particles
# x- vector of states
# y- vector of measurements
# x0- starting value
# sigma- standard deviation of process 
# sigma.meas- standard deviation of measurment
# resample.type - type of resampling 

particle.filter.pop<-function(N,x,y,x0,sigma=1,sigma.meas=1,resample.type="multinomial",N.thr.per=0.5){
  
  #threshold N
  N.thr<-as.integer(N.thr.per*N)
  
  #set up filter
  x0.pf<-rnorm(N,mean=x0,sd=sigma)
  x.pf <- matrix(NA,ncol=N,nrow=tau)
  
  #set up diagnostics
  #-adds variables for diagnoistics
  diagnostic.init(tau)
  
  
  #local function to choose resampling method
  FUN.resample<-function(N,x,w){
    switch(resample.type,
           multinomial = mult.resample(N=N,x=x,w=w) ,
           stratified = strat.resample(N=N,x=x,w=w) ,
           systematic = sys.resample(N=N,x=x,w=w),
           residual= res.resample(N=N,x=x,w=w),
           standard=sample(x=x,size=N,replace=TRUE,prob=w)
    )
    
  }
  
  # Initialize
  x.pf[1, ] <- f(x0.pf,1)+rnorm(N,sd=sigma) #project particles from x0 to the first time step
  w.tilde<- dnorm(y[1],mean=g(x.pf[1,]),sd=sigma.meas) #the previous weight would be equal weights 1/N depending on how x0 was chosen
  w<-w.tilde/sum(w.tilde) #normalise weights
  
  
  
  #compute diagnostics
  #- mean, effective sample size, confidence intervals 
  m[1]<-sum(w*x.pf[1,])
  MSE.k[1]<-(x[1]-m[1])^2
  N.eff[1]<-1/sum(w^2) 
  conf[1,]<-quantile(x.pf[1,],probs=c(0.025,0.975))
  
  for (t in 2:(tau)) {
    
    #project  state forward
    x.pf[t, ] <-f(x.pf[t-1,],t)+rnorm(N,sd=sigma)
    
    #compute weights
    w.tilde <- dnorm(y[t], mean=g(x.pf[t, ]),sd=sigma.meas)*w #when add w resampling does not help
    w <- w.tilde/sum(w.tilde) #normalise weights
    
    
    #compute diagnostics
    m[t]<-sum(w*x.pf[t,])
    MSE.k[t]<-(x[t]-m[t])^2
    conf[t,]<-quantile(x.pf[t,],probs=c(0.025,0.975))
    N.eff[t]<-1/sum(w^2)
    
    #resample (bootstrap means resample everytime)
    if(N.eff[t]<N.thr){
      x.pf[t,]<-FUN.resample(N=N,x=x.pf[t,],w=w)
      w<-rep(1/N,N)
      res.point[t]<-1 #add "yes" for resampling point
    }
    
    
    
  }
  
  #return output
  return(list(x.out=x,y.out=y, x.pf.out=x.pf, m.out=m,MSE.k.out=MSE.k, conf.out=conf, 
              N.eff.out=N.eff,w.out=w,res.point.out=res.point))
  
}


#--Generic particle filter(prior proposal) Population example Paths--
# N- No of particles
# x- vector of states
# y- vector of measurements
# x0- starting value
# sigma- standard deviation of process 
# sigma.meas- standard deviation of measurment
# resample.type - type of resampling 

particle.filter.pop.path<-function(N,x,y,x0,sigma=1,sigma.meas=1,resample.type="multinomial",N.thr.per=0.5){
  
  
  tau<-length(x)
  
  #threshold N
  N.thr<-as.integer(N.thr.per*N)
  
  #set up filter
  x0.pf<-rnorm(N,mean=x0,sd=sigma)
  x.pf <- matrix(NA,ncol=N,nrow=tau)
  
  #set up diagnostics
  #-adds variables for diagnoistics
  diagnostic.init(tau)
  
  
  #local function to choose resampling method
  FUN.resample<-function(N,x,w){
    switch(resample.type,
           multinomial = mult.resample(N=N,x=x,w=w) ,
           stratified = strat.resample(N=N,x=x,w=w) ,
           systematic = sys.resample(N=N,x=x,w=w),
           residual= res.resample(N=N,x=x,w=w),
           standard=sample(x=x,size=N,replace=TRUE,prob=w)
    )
    
  }
  
  # Initialize
  x.pf[1, ] <- f(x0.pf,1)+rnorm(N,sd=sigma) #project particles from x0 to the first time step
  w.tilde<- dnorm(y[1],mean=g(x.pf[1,]),sd=sigma.meas) #the previous weight would be equal weights 1/N depending on how x0 was chosen
  w<-w.tilde/sum(w.tilde) #normalise weights
  
  
  
  #compute diagnostics
  #- mean, effective sample size, confidence intervals 
  m[1]<-sum(w*x.pf[1,])
  MSE.k[1]<-(x[1]-m[1])^2
  N.eff[1]<-1/sum(w^2) 
  conf[1,]<-quantile(x.pf[1,],probs=c(0.025,0.975))
  
  for (t in 2:(tau)) {
    
    #project  state forward
    x.pf[t, ] <-f(x.pf[t-1,],t)+rnorm(N,sd=sigma)
    
    #compute weights
    w.tilde <- dnorm(y[t], mean=g(x.pf[t, ]),sd=sigma.meas)*w #when add w resampling does not help
    w <- w.tilde/sum(w.tilde) #normalise weights
    
    
    #compute diagnostics
    m[t]<-sum(w*x.pf[t,])
    MSE.k[t]<-(x[t]-m[t])^2
    conf[t,]<-quantile(x.pf[t,],probs=c(0.025,0.975))
    N.eff[t]<-1/sum(w^2)
    
    #resample (bootstrap means resample everytime)
    if(N.eff[t]<N.thr){
      s<-FUN.resample(N=N,x=1:N,w=w)
      x.pf<-x.pf[,s]
      w<-rep(1/N,N)
      res.point[t]<-1 #add "yes" for resampling point
    }
    
    
    
    
  }
  
  #return output
  return(list(x.out=x,y.out=y, x.pf.out=x.pf, m.out=m,MSE.k.out=MSE.k, conf.out=conf, 
              N.eff.out=N.eff,w.out=w,res.point.out=res.point))
  
}


#--Likelihood particle filter --
# N- No of particles
# x- vector of states
# y- vector of measurements
# x0- starting value
# sigma- standard deviation of process 
# sigma.meas- standard deviation of measurment
# resample.type - type of resampling 

likelihood.filter<-function(N,x,y,x0,sigma=1,sigma.meas=1,resample.type="multinomial",N.thr.per=0.5){
  
  #threshold N
  N.thr<-as.integer(N.thr.per*N)
  
  #set up filter
  x0.pf<-rnorm(N,mean=x0,sd=sigma) #suppose this is the correct sigma
  x.pf <- matrix(NA,ncol=N,nrow=tau)
  
  #set up diagnostics
  #-adds variables for diagnoistics
  diagnostic.init(tau)
  
  #local function to choose resampling method
  FUN.resample<-function(N,x,w){
    switch(resample.type,
           multinomial = mult.resample(N=N,x=x,w=w) ,
           stratified = strat.resample(N=N,x=x,w=w) ,
           systematic = sys.resample(N=N,x=x,w=w),
           residual= res.resample(N=N,x=x,w=w),
           standard=sample(x=x,size=N,replace=TRUE,prob=w)
    )
  }
  
  # Initialize
  x.pf[1,]<-y[1]-rnorm(N,sd=sigma.meas)
  w.tilde<-dnorm(x.pf[1,],mean=x0,sd=sigma)
  w<-w.tilde/sum(w.tilde) 
  
  
  #compute diagnostics
  #- mean, effective sample size, confidence intervals 
  m[1]<-sum(w*x.pf[1,]) 
  MSE.k[1]<-(x[1]-m[1])^2
  N.eff[1]<-1/sum(w^2) 
  conf[1,]<-quantile(x.pf[1,],probs=c(0.025,0.975))
  
  
  for (t in 2:(tau)) {
    
    
    #project  state forward
    x.pf[t, ] <-y[t]-rnorm(N,sd=sigma.meas)
    
    #compute weights
    w.tilde <- dnorm(x.pf[t,], mean=x.pf[t-1,],sd=sigma)*w
    w <- w.tilde/sum(w.tilde) 
    
    #compute diagnostics
    m[t]<-sum(w*x.pf[t,])
    MSE.k[t]<-(x[t]-m[t])^2
    conf[t,]<-quantile(x.pf[t,],probs=c(0.025,0.975))
    N.eff[t]<-1/sum(w^2)
    
    #resample (bootstrap means resample everytime)
    if(N.eff[t]<N.thr){
      x.pf[t,]<-FUN.resample(N=N,x=x.pf[t,],w=w) 
      w<-rep(1/N,N)
      res.point[t]<-1 #add "yes" for resampling point
    }
    
    
    
  }
  
  #return output
  return(list(x.out=x,y.out=y, x.pf.out=x.pf, m.out=m,MSE.k.out=MSE.k, 
              conf.out=conf,N.eff.out=N.eff,w.out=w,res.point.out=res.point))
  
}

#--Auxiliary particle filter--
# N- No of particles
# x- vector of states
# y- vector of measurements
# x0- starting value
# sigma- standard deviation of process 
# sigma.meas- standard deviation of measurment
# resample.type - type of resampling (multinomial,stratified,systematic,residual,standard)

auxiliary.filter<-function(N,x,y,x0,sigma=1,sigma.meas=1,resample.type="multinomial"){
  
  #set up filter
  x0.pf<-rnorm(N,mean=x0,sd=sigma) #suppose this is the correct sigma
  x.pf <- matrix(NA,ncol=N,nrow=tau)
  
  #set up diagnostics
  #-adds variables for diagnoistics
  diagnostic.init(tau)
  
  #local function to choose resampling method
  FUN.resample<-function(N,x,w){
    switch(resample.type,
           multinomial = mult.resample(N=N,x=x,w=w) ,
           stratified = strat.resample(N=N,x=x,w=w) ,
           systematic = sys.resample(N=N,x=x,w=w),
           residual= res.resample(N=N,x=x,w=w),
           standard=sample(x=x,size=N,replace=TRUE,prob=w)
    )
  }
  
  # Initialize
  mu<-x0.pf #if use true mean x0 must be replicated N times
  #compute parent weights
  v<-dnorm(y[1],mu,sd=sigma.meas)
  #normalise parent weights
  v<-v/ sum(v)
  #sample parents
  samp<-sample(1:N,size=N,replace=TRUE,prob=v)
  #sample from prior 
  x.pf[1,]<-rnorm(N,x0.pf[samp],sd=sigma) 
  #set weights
  w.tilde<-dnorm(y[1],x.pf[1,],sd=sigma.meas)/dnorm(y[1],x0.pf[samp],sd=sigma.meas)
  #normalise weights
  w<-w.tilde/sum(w.tilde)
  
  
  #compute diagnostics
  #- mean, effective sample size, confidence intervals 
  m[1]<-sum(w*x.pf[1,]) 
  MSE.k[1]<-(x[1]-m[1])^2
  N.eff[1]<-1/sum(w^2) 
  conf[1,]<-quantile(x.pf[1,],probs=c(0.025,0.975))
  
  
  for (t in 2:(tau)) {
    
    
    # calculate mean of previous statistics
    mu<-x.pf[t-1,]
    v<-w*dnorm(y[t],mu,sd=sigma.meas)
    #normalise parent weights
    v<-v/ sum(v)
    #sample parent trajectories
    samp<-sample(1:N,size=N,replace=TRUE,prob=v)
    #sample from prior 
    x.pf[t,]<-rnorm(N,x.pf[t-1,samp],sd=sigma) 
    #set weights
    w.tilde<-dnorm(y[t],x.pf[t,],sd=sigma.meas)/dnorm(y[t],x.pf[t-1,samp],sd=sigma.meas)
    #normalise weights
    w<-w.tilde/sum(w.tilde)
    
    #compute diagnostics
    m[t]<-sum(w*x.pf[t,])
    N.eff[t]<-1/sum(w^2) 
    MSE.k[t]<-(x[t]-m[t])^2
    conf[t,]<-quantile(x.pf[t,],probs=c(0.025,0.975))
    
  }
  
  #return output
  return(list(x.out=x,y.out=y, x.pf.out=x.pf, m.out=m,MSE.k.out=MSE.k, 
              conf.out=conf,N.eff.out=N.eff,w.out=w,res.point.out=res.point))
  
}


#--Auxiliary particle filter Population Model--
# N- No of particles
# x- vector of states
# y- vector of measurements
# x0- starting value
# sigma- standard deviation of process 
# sigma.meas- standard deviation of measurment
# resample.type - type of resampling (multinomial,stratified,systematic,residual,standard)

auxiliary.filter.pop<-function(N,x,y,x0,sigma=1,sigma.meas=1,resample.type="multinomial"){
  
  #set up filter
  x0.pf<-rnorm(N,mean=x0,sd=sigma) #suppose this is the correct sigma
  x.pf <- matrix(NA,ncol=N,nrow=tau)
  
  #set up diagnostics
  #-adds variables for diagnoistics
  diagnostic.init(tau)
  
  #local function to choose resampling method
  FUN.resample<-function(N,x,w){
    switch(resample.type,
           multinomial = mult.resample(N=N,x=x,w=w) ,
           stratified = strat.resample(N=N,x=x,w=w) ,
           systematic = sys.resample(N=N,x=x,w=w),
           residual= res.resample(N=N,x=x,w=w),
           standard=sample(x=x,size=N,replace=TRUE,prob=w)
    )
  }
  
  # Initialize
  mu<-f(x0.pf,1) #rep(mean(f(x0.pf,1)),N) #if use true mean x0 must be replicated N times
  #compute parent weights
  v<-dnorm(y[1],g(mu),sd=sigma.meas)
  #normalise parent weights
  v<-v/ sum(v)
  #sample parents
  samp<-sample(1:N,size=N,replace=TRUE,prob=v)
  #sample from prior 
  x.pf[1,]<-rnorm(N,f(x0.pf[samp],1),sd=sigma) 
  #set weights
  w.tilde<-dnorm(y[1],g(x.pf[1,]),sd=sigma.meas)/dnorm(y[1],g(x0.pf[samp]),sd=sigma.meas)
  #normalise weights
  w<-w.tilde/sum(w.tilde)
  
  #compute diagnostics
  #- mean, effective sample size, confidence intervals 
  m[1]<-sum(w*x.pf[1,]) 
  MSE.k[1]<-(x[1]-m[1])^2
  N.eff[1]<-1/sum(w^2) 
  conf[1,]<-quantile(x.pf[1,],probs=c(0.025,0.975))
  
  
  for (t in 2:(tau)) {
    
    
    # calculate mean of previous statistics
    mu<-f(x.pf[t-1,],t) #rep(mean(f(x.pf[t-1,],t)),N)
    v<-w*dnorm(y[t],g(mu),sd=sigma.meas)
    #normalise parent weights
    v<-v/ sum(v)
    #sample parent trajectories
    samp<-sample(1:N,size=N,replace=TRUE,prob=v)
    #sample from prior 
    x.pf[t,]<-rnorm(N,f(x.pf[t-1,samp],t),sd=sigma) 
    #set weights
    w.tilde<-dnorm(y[t],g(x.pf[t,]),sd=sigma.meas)/dnorm(y[t],g(x.pf[t-1,samp]),sd=sigma.meas)
    #normalise weights
    w<-w.tilde/sum(w.tilde)
    
    #compute diagnostics
    m[t]<-sum(w*x.pf[t,])
    N.eff[t]<-1/sum(w^2) 
    MSE.k[t]<-(x[t]-m[t])^2
    conf[t,]<-quantile(x.pf[t,],probs=c(0.025,0.975))
    
  }
  
  #return output
  return(list(x.out=x,y.out=y, x.pf.out=x.pf, m.out=m,MSE.k.out=MSE.k, 
              conf.out=conf,N.eff.out=N.eff,w.out=w,res.point.out=res.point))
  
}


#--Regularized particle filter(prior proposal) NOT USED--
#random walk
# we use optimal bin width, given Epanechnikov kernel, assuming filtered distribution is gaussian
# N- No of particles
# x- vector of states
# y- vector of measurements
# x0- starting value
# sigma- standard deviation of process 
# sigma.meas- standard deviation of measurment
# resample.type - type of resampling (multinomial,stratified,systematic,residual,standard)

regularized.filter<-function(N,x,y,x0,sigma=1,sigma.meas=1,resample.type="multinomial",N.thr.per=0.5){
  
  #threshold N
  N.thr<-as.integer(N.thr.per*N)
  
  #set up filter
  x0.pf<-rnorm(N,mean=x0,sd=sigma) #suppose this is the correct sigma
  x.pf <- matrix(NA,ncol=N,nrow=tau)
  
  #set up diagnostics
  #-adds variables for diagnoistics
  diagnostic.init(tau)
  
  # Compute constants
  dim<-1 #dimesnion of state 
  c<-2  #volume of a unit sphere in [dim] dimensions
  A<-((8/c)*(dim+4)*2*sqrt(pi))^(1/(dim+4))
  h<-A*N^(1/(dim+4)) #compute optimal bandwidth
  
  #local function to choose resampling method
  FUN.resample<-function(N,x,w){
    switch(resample.type,
           multinomial = mult.resample(N=N,x=x,w=w) ,
           stratified = strat.resample(N=N,x=x,w=w) ,
           systematic = sys.resample(N=N,x=x,w=w),
           residual= res.resample(N=N,x=x,w=w),
           standard=sample(x=x,size=N,replace=TRUE,prob=w)
    )
  }
  
  # Initialize
  x.pf[1, ] <- x0.pf+rnorm(N,sd=sigma) #project particles from x0 to the first time step
  w.tilde<- dnorm(y[1],mean=x.pf[1,],sd=sigma.meas) #weights of x0 is equal
  w<-w.tilde/sum(w.tilde) 
  
  
  #compute diagnostics
  #- mean, effective sample size, confidence intervals 
  m[1]<-sum(w*x.pf[1,]) 
  MSE.k[1]<-(x[1]-m[1])^2
  N.eff[1]<-1/sum(w^2) 
  conf[1,]<-quantile(x.pf[1,],probs=c(0.025,0.975))
  
  
  for (t in 2:(tau)) {
    
    
    #project  state forward
    x.pf[t, ] <-x.pf[t-1,]+rnorm(N,sd=sigma)
    
    #compute weights
    w.tilde <- dnorm(y[t], mean=x.pf[t, ],sd=sigma.meas)*w 
    w <- w.tilde/sum(w.tilde) 
    
    #compute diagnostics
    m[t]<-sum(w*x.pf[t,])
    MSE.k[t]<-(x[t]-m[t])^2
    conf[t,]<-quantile(x.pf[t,],probs=c(0.025,0.975))
    N.eff[t]<-1/sum(w^2)
    
    #compute empirical covariance
    #d<-sqrt(sum(w*x.pf[t,]^2))
    d<-1
    
    
    
    #resample (bootstrap means resample everytime)
    if(N.eff[t]<N.thr){
      x.pf[t,]<-FUN.resample(N=N,x=x.pf[t,],w=w) 
      w<-rep(1/N,N)
      res.point[t]<-1 #add "yes" for resampling point
      
      #fit
      fit<-density(x.pf[t,])
      x.pf[t,]<- rnorm(N, sample(x.pf[t,], size = N, replace = TRUE), fit$bw)
      
      #sample from Epanechnikov
      #tsi<-apply(matrix(runif(3*N, -1, 1), 3), 2, median)
      
      
      #jitter particles
      #x.pf[t,]<-x.pf[t,]+tsi*d*h
    }
    
    
    
  }
  
  #return output
  return(list(x.out=x,y.out=y, x.pf.out=x.pf, m.out=m,MSE.k.out=MSE.k, 
              conf.out=conf,N.eff.out=N.eff,w.out=w,res.point.out=res.point))
  
}


#--Regularized particle filter Population Model(prior proposal)--
# we use optimal bin width, given Epanechnikov kernel, assuming filtered distribution is gaussian
# N- No of particles
# x- vector of states
# y- vector of measurements
# x0- starting value
# sigma- standard deviation of process 
# sigma.meas- standard deviation of measurment
# resample.type - type of resampling (multinomial,stratified,systematic,residual,standard)

regularized.filter.pop<-function(N,x,y,x0,sigma=1,sigma.meas=1,resample.type="multinomial",N.thr.per=0.5){
  
  #threshold N
  N.thr<-as.integer(N.thr.per*N)
  
  #set up filter
  x0.pf<-rnorm(N,mean=x0,sd=sigma) #suppose this is the correct sigma
  x.pf <- matrix(NA,ncol=N,nrow=tau)
  
  #set up diagnostics
  #-adds variables for diagnoistics
  diagnostic.init(tau)
  
  # Compute constants
  dim<-1 #dimesnion of state 
  c<-2  #volume of a unit sphere in [dim] dimensions
  A<-((8/c)*(dim+4)*2*sqrt(pi))^(1/(dim+4))
  h<-A*N^(1/(dim+4)) #compute optimal bandwidth
  
  #local function to choose resampling method
  FUN.resample<-function(N,x,w){
    switch(resample.type,
           multinomial = mult.resample(N=N,x=x,w=w) ,
           stratified = strat.resample(N=N,x=x,w=w) ,
           systematic = sys.resample(N=N,x=x,w=w),
           residual= res.resample(N=N,x=x,w=w),
           standard=sample(x=x,size=N,replace=TRUE,prob=w)
    )
  }
  
  # Initialize
  x.pf[1, ] <- f(x0.pf,1)+rnorm(N,sd=sigma) #project particles from x0 to the first time step
  w.tilde<- dnorm(y[1],mean=g(x.pf[1,]),sd=sigma.meas) #weights of x0 is equal
  w<-w.tilde/sum(w.tilde) 
  
  
  #compute diagnostics
  #- mean, effective sample size, confidence intervals 
  m[1]<-sum(w*x.pf[1,]) 
  MSE.k[1]<-(x[1]-m[1])^2
  N.eff[1]<-1/sum(w^2) 
  conf[1,]<-quantile(x.pf[1,],probs=c(0.025,0.975))
  
  
  for (t in 2:(tau)) {
    
    
    #project  state forward
    x.pf[t, ] <-f(x.pf[t-1,],t)+rnorm(N,sd=sigma)
    
    #compute weights
    w.tilde <- dnorm(y[t], mean=g(x.pf[t, ]),sd=sigma.meas)*w 
    w <- w.tilde/sum(w.tilde) 
    
    #compute diagnostics
    m[t]<-sum(w*x.pf[t,])
    MSE.k[t]<-(x[t]-m[t])^2
    conf[t,]<-quantile(x.pf[t,],probs=c(0.025,0.975))
    N.eff[t]<-1/sum(w^2)
    
    #compute empirical covariance
    d<-sqrt(sum(w*x.pf[t,]^2))
    
    
    #resample (bootstrap means resample everytime)
    if(N.eff[t]<N.thr){
      x.pf[t,]<-FUN.resample(N=N,x=x.pf[t,],w=w) 
      #samp<-FUN.resample(N=N,x=1:N,w=w) 
      #x.pf[t,]<-x.pf[t,samp]
      w<-rep(1/N,N)
      res.point[t]<-1 #add "yes" for resampling point
      
      #fit
      fit<-density(x.pf[t,])
      x.pf[t,]<- rnorm(N, sample(x.pf[t,], size = N, replace = TRUE), fit$bw)
      
      #sample from Epanechnikov
      #tsi<-apply(matrix(runif(3*N, -1, 1), 3), 2, median)
      
      #jitter particles
      #x.pf[t,]<-x.pf[t,]+tsi*samp*h
      #x.pf[t,]<-x.pf[t,]+tsi*d*h
    }
    
    
    
  }
  
  #return output
  return(list(x.out=x,y.out=y, x.pf.out=x.pf, m.out=m,MSE.k.out=MSE.k, 
              conf.out=conf,N.eff.out=N.eff,w.out=w,res.point.out=res.point))
  
}



#--Regularized particle filter Population Model(prior proposal)--
# we use optimal bin width, given Epanechnikov kernel, assuming filtered distribution is gaussian
# N- No of particles
# x- vector of states
# y- vector of measurements
# x0- starting value
# sigma- standard deviation of process 
# sigma.meas- standard deviation of measurment
# resample.type - type of resampling (multinomial,stratified,systematic,residual,standard)

regularized.filter.pop.path<-function(N,x,y,x0,sigma=1,sigma.meas=1,resample.type="multinomial",N.thr.per=0.5){
  
  #threshold N
  N.thr<-as.integer(N.thr.per*N)
  
  #set up filter
  x0.pf<-rnorm(N,mean=x0,sd=sigma) #suppose this is the correct sigma
  x.pf <- matrix(NA,ncol=N,nrow=tau)
  
  #set up diagnostics
  #-adds variables for diagnoistics
  diagnostic.init(tau)
  
  # Compute constants
  dim<-1 #dimesnion of state 
  c<-2  #volume of a unit sphere in [dim] dimensions
  A<-((8/c)*(dim+4)*2*sqrt(pi))^(1/(dim+4))
  h<-A*N^(1/(dim+4)) #compute optimal bandwidth
  
  #local function to choose resampling method
  FUN.resample<-function(N,x,w){
    switch(resample.type,
           multinomial = mult.resample(N=N,x=x,w=w) ,
           stratified = strat.resample(N=N,x=x,w=w) ,
           systematic = sys.resample(N=N,x=x,w=w),
           residual= res.resample(N=N,x=x,w=w),
           standard=sample(x=x,size=N,replace=TRUE,prob=w)
    )
  }
  
  # Initialize
  x.pf[1, ] <- f(x0.pf,1)+rnorm(N,sd=sigma) #project particles from x0 to the first time step
  w.tilde<- dnorm(y[1],mean=g(x.pf[1,]),sd=sigma.meas) #weights of x0 is equal
  w<-w.tilde/sum(w.tilde) 
  
  
  #compute diagnostics
  #- mean, effective sample size, confidence intervals 
  m[1]<-sum(w*x.pf[1,]) 
  MSE.k[1]<-(x[1]-m[1])^2
  N.eff[1]<-1/sum(w^2) 
  conf[1,]<-quantile(x.pf[1,],probs=c(0.025,0.975))
  
  
  for (t in 2:(tau)) {
    
    
    #project  state forward
    x.pf[t, ] <-f(x.pf[t-1,],t)+rnorm(N,sd=sigma)
    
    #compute weights
    w.tilde <- dnorm(y[t], mean=g(x.pf[t, ]),sd=sigma.meas)*w 
    w <- w.tilde/sum(w.tilde) 
    
    #compute diagnostics
    m[t]<-sum(w*x.pf[t,])
    MSE.k[t]<-(x[t]-m[t])^2
    conf[t,]<-quantile(x.pf[t,],probs=c(0.025,0.975))
    N.eff[t]<-1/sum(w^2)
    
    #compute empirical covariance
    d<-sqrt(sum(w*x.pf[t,]^2))
    
    
    #resample (bootstrap means resample everytime)
    if(N.eff[t]<N.thr){
      #x.pf[t,]<-FUN.resample(N=N,x=x.pf[t,],w=w) 
      samp<-FUN.resample(N=N,x=1:N,w=w) 
      x.pf[t,]<-x.pf[t,samp]
      w<-rep(1/N,N)
      res.point[t]<-1 #add "yes" for resampling point
      
      #fit
      fit<-density(x.pf[t,])
      x.pf[t,]<- rnorm(N, sample(x.pf[t,], size = N, replace = TRUE), fit$bw)
      
      #sample from Epanechnikov
      #tsi<-apply(matrix(runif(3*N, -1, 1), 3), 2, median)
      
      #jitter particles
      #x.pf[t,]<-x.pf[t,]+tsi*samp*h
      #x.pf[t,]<-x.pf[t,]+tsi*d*h
    }
    
    
    
  }
  
  #return output
  return(list(x.out=x,y.out=y, x.pf.out=x.pf, m.out=m,MSE.k.out=MSE.k, 
              conf.out=conf,N.eff.out=N.eff,w.out=w,res.point.out=res.point))
  
}



#--Rejection particle filter(prior proposal)--

# N- No of particles
# x- vector of states
# y- vector of measurements
# x0- starting value
# sigma- standard deviation of process 
# sigma.meas- standard deviation of measurment
# resample.type - type of resampling (multinomial,stratified,systematic,residual,standard)

rejection.filter<-function(N,x,y,x0,sigma=1,sigma.meas=1,Mnorm=1){
  
  #set up filter
  x0.pf<-rnorm(N,mean=x0,sd=sigma) #suppose this is the correct sigma
  x.pf <- matrix(NA,ncol=N,nrow=tau)
  
  #set up diagnostics
  #-adds variables for diagnoistics
  diagnostic.init(tau)
  
  # Initialize
  
  #project particles from x0 to the first time step
  x.pf[1, ] <- x0.pf+rnorm(N,sd=sigma) 
  noise.j<-rnorm(N,sd=sigma) 
  grid<-expand.grid(noise.j, x0.pf) #j-noise, i- coord (t-1 ) particles   #each row is section is first x0.pf with every single noise
  grid.sum<-rowSums(grid)  #each row is sim of x0.pf with every single noise
  c<-(1/N^2)*sum(dnorm(y[1],grid.sum))
  index<-rep(1:N, each=N) #apply on each section on each x0.pf
  ci<-sapply(1:N, function(s) mean(dnorm(y[1],grid.sum[index==s]))) 
  lambda<-ci/N*c #version of weights
  x.old<-sample(x0.pf,  N,replace=TRUE,prob=lambda)
  
  u<-1
  for(i in 1:N){
    a<-0
    while(u>a){
      z<-x.old[i]+rnorm(1,0,sd=sigma) #make a proposal
      a<-( dnorm(y[1],z,sd=sigma.meas)*dnorm(z,x.old[i], sd=sigma)/dnorm(z,x.old[i],sd=sigma) )/Mnorm
      #print(paste("u",u,"x.old",x.old[i],"i",i,"a",a))
      u<-runif(1,0,1)
      if(a>1){
        print("alert. increase M")
      }
    }
    x.pf[1,i]<-z
  }
  
  
  
  #compute diagnostics
  #- mean, effective sample size, confidence intervals 
  m[1]<-mean(x.pf[1,]) 
  MSE.k[1]<-(x[1]-m[1])^2
  conf[1,]<-quantile(x.pf[1,],probs=c(0.025,0.975))
  
  
  for (t in 2:(tau)) {
    
    print(t)
    x.pf[t, ] <- x.pf[t-1,]+rnorm(N,sd=sigma) #project particles from x0 to the first time step
    noise.j<-rnorm(N,sd=sigma) 
    grid<-expand.grid(noise.j, x.pf[t-1,]) #j-noise, i- coord (t-1 ) particles
    grid.sum<-rowSums(grid)
    c<-(1/N^2)*sum(dnorm(y[t],grid.sum))
    index<-rep(1:N, each=N)
    ci<-sapply(1:N, function(s) mean(dnorm(y[t],grid.sum[index==s])))
    lambda<-ci/N*c #version of weights
    x.old<-sample(x.pf[t-1,],N,replace=TRUE,prob=lambda)
    
    u<-1
    for(i in 1:N){
      a<-0
      while(u>a){
        z<-x.old[i]+rnorm(1,0,sd=sigma)
        a<- (dnorm(y[t],z,sd=sigma.meas)*dnorm(z,x.old[i], sd=sigma)  /dnorm(z,x.old[i],sd=sigma) )/Mnorm #x.pf[t-1,i]
        
        #print(paste("t",t,"u",u,"x.old",x.old[i],"i",i,"a",a))
        if(a>1){
          print("alert. increase M")
        }
        u<-runif(1,0,1)
      }
      
      x.pf[t,i]<-z
    }
    
    
    #compute diagnostics
    m[t]<-mean(x.pf[t,])
    MSE.k[t]<-(x[t]-m[t])^2
    conf[t,]<-quantile(x.pf[t,],probs=c(0.025,0.975))
    
    
  }
  
  #return output
  return(list(x.out=x,y.out=y, x.pf.out=x.pf, m.out=m,MSE.k.out=MSE.k, 
              conf.out=conf))
  
}



#--Rejection particle filter Population(prior proposal)--

# N- No of particles
# x- vector of states
# y- vector of measurements
# x0- starting value
# sigma- standard deviation of process 
# sigma.meas- standard deviation of measurment
# resample.type - type of resampling (multinomial,stratified,systematic,residual,standard)

rejection.filter.pop<-function(N,x,y,x0,sigma=1,sigma.meas=1,Mnorm=1){
  
  #set up filter
  x0.pf<-rnorm(N,mean=x0,sd=sigma) #suppose this is the correct sigma
  x.pf <- matrix(NA,ncol=N,nrow=tau)
  
  #set up diagnostics
  #-adds variables for diagnoistics
  diagnostic.init(tau)
  
  # Initialize
  
  #project particles from x0 to the first time step
  x.pf[1, ] <- f(x0.pf,1)+rnorm(N,sd=sigma) 
  noise.j<-rnorm(N,sd=sigma) 
  grid<-expand.grid(noise.j, f(x0.pf,1)) #j-noise, i- coord (t-1 ) particles   #each row is section is first x0.pf with every single noise
  grid.sum<-rowSums(grid)  #each row is sim of x0.pf with every single noise
  c<-(1/N^2)*sum(dnorm(y[1],g(grid.sum)))
  index<-rep(1:N, each=N) #apply on each section on each x0.pf
  ci<-sapply(1:N, function(s) mean(dnorm(y[1],g(grid.sum[index==s]))) )
  lambda<-ci/(N*c) #version of weights
  x.old<-sample(x0.pf,  N,replace=TRUE,prob=lambda)
  
  u<-1
  for(i in 1:N){
    a<-0
    while(log(u)>log(a)){
      z<-f(x.old[i],1)+rnorm(1,0,sd=sigma) #make a proposal
      a<-( dnorm(y[1],g(z),sd=sigma.meas)*dnorm(z,f(x.old[i],1), sd=sigma)/dnorm(z,f(x.old[i],1),sd=sigma) )/Mnorm
      print(paste("u",u,"x.old",x.old[i],"i",i,"a",a))
      u<-runif(1,0,1)
      if(a>1){
        print("alert. increase M")
      }
    }
    x.pf[1,i]<-z
  }
  
  
  
  #compute diagnostics
  #- mean, effective sample size, confidence intervals 
  m[1]<-mean(x.pf[1,]) 
  MSE.k[1]<-(x[1]-m[1])^2
  conf[1,]<-quantile(x.pf[1,],probs=c(0.025,0.975))
  
  
  for (t in 2:(tau)) {
    
    x.pf[t, ] <- f(x.pf[t-1,],t)+rnorm(N,sd=sigma) #project particles from x0 to the first time step
    noise.j<-rnorm(N,sd=sigma) 
    grid<-expand.grid(noise.j, f(x.pf[t-1,],t) )#j-noise, i- coord (t-1 ) particles
    grid.sum<-rowSums(grid)
    c<-(1/N^2)*sum(dnorm(y[t],g(grid.sum)))
    index<-rep(1:N, each=N)
    ci<-sapply(1:N, function(s) mean(dnorm(y[t],g(grid.sum[index==s]))))
    lambda<-ci/(N*c) #version of weights
    x.old<-sample(x.pf[t-1,],N,replace=TRUE,prob=lambda)
    
    u<-1
    for(i in 1:N){
      a<-0
      while(u>a){
        z<-f(x.old[i],t)+rnorm(1,0,sd=sigma)
        a<- (dnorm(y[t],g(z),sd=sigma.meas)*dnorm(z,f(x.old[i],t), sd=sigma)  /dnorm(z,f(x.old[i],t),sd=sigma) )/Mnorm #x.pf[t-1,i]
        
        print(paste("t",t,"u",u,"x.old",x.old[i],"i",i,"a",a))
        if(a>1){
          print("alert. increase M")
        }
        u<-runif(1,0,1)
      }
      
      x.pf[t,i]<-z
    }
    
    
    #compute diagnostics
    m[t]<-mean(x.pf[t,])
    MSE.k[t]<-(x[t]-m[t])^2
    conf[t,]<-quantile(x.pf[t,],probs=c(0.025,0.975))
    
    
  }
  
  #return output
  return(list(x.out=x,y.out=y, x.pf.out=x.pf, m.out=m,MSE.k.out=MSE.k, 
              conf.out=conf))
  
}

#--Rejection particle filter Population(prior proposal)--

# N- No of particles
# x- vector of states
# y- vector of measurements
# x0- starting value
# sigma- standard deviation of process 
# sigma.meas- standard deviation of measurment
# resample.type - type of resampling (multinomial,stratified,systematic,residual,standard)

rejection.filter.pop.log<-function(N,x,y,x0,sigma=1,sigma.meas=1,Mnorm=1){
  
  #set up filter
  x0.pf<-rnorm(N,mean=x0,sd=sigma) #suppose this is the correct sigma
  x.pf <- matrix(NA,ncol=N,nrow=tau)
  
  #set up diagnostics
  #-adds variables for diagnoistics
  diagnostic.init(tau)
  
  # Initialize
  
  #project particles from x0 to the first time step
  x.pf[1, ] <- f(x0.pf,1)+rnorm(N,sd=sigma) 
  noise.j<-rnorm(N,sd=sigma) 
  grid<-expand.grid(noise.j, f(x0.pf,1)) #j-noise, i- coord (t-1 ) particles   #each row is section is first x0.pf with every single noise
  grid.sum<-rowSums(grid)  #each row is sim of x0.pf with every single noise
  c<-(1/N^2)*sum(dnorm(y[1],g(grid.sum)))
  index<-rep(1:N, each=N) #apply on each section on each x0.pf
  ci<-sapply(1:N, function(s) mean(dnorm(y[1],g(grid.sum[index==s]))) )
  lambda<-ci/(N*c) #version of weights
  x.old<-sample(x0.pf,  N,replace=TRUE,prob=lambda)
  
  u<-1
  for(i in 1:N){
    a<-0
    while(u>a){
      z<-f(x.old[i],1)+rnorm(1,0,sd=sigma) #make a proposal
      a<-  log(dnorm(y[1],g(z),sd=sigma.meas)) -log(Mnorm)
      print(paste("u",u,"x.old",x.old[i],"i",i,"a",a))
      u<-log(runif(1,0,1))
      if(a>1){
        print("alert. increase M")
      }
    }
    x.pf[1,i]<-z
  }
  
  #compute diagnostics
  #- mean, effective sample size, confidence intervals 
  m[1]<-mean(x.pf[1,]) 
  MSE.k[1]<-(x[1]-m[1])^2
  conf[1,]<-quantile(x.pf[1,],probs=c(0.025,0.975))
  
  
  for (t in 2:(tau)) {
    
    x.pf[t, ] <- f(x.pf[t-1,],t)+rnorm(N,sd=sigma) #project particles from x0 to the first time step
    noise.j<-rnorm(N,sd=sigma) 
    grid<-expand.grid(noise.j, f(x.pf[t-1,],t) )#j-noise, i- coord (t-1 ) particles
    grid.sum<-rowSums(grid)
    c<-(1/N^2)*sum(dnorm(y[t],g(grid.sum)))
    index<-rep(1:N, each=N)
    ci<-sapply(1:N, function(s) mean(dnorm(y[t],g(grid.sum[index==s]))))
    lambda<-ci/(N*c) #version of weights
    x.old<-sample(x.pf[t-1,],N,replace=TRUE,prob=lambda)
    
    u<-1
    for(i in 1:N){
      a<-0
      while(u>a){
        z<-f(x.old[i],t)+rnorm(1,0,sd=sigma)
        a<- log(dnorm(y[t],g(z),sd=sigma.meas)) -log(Mnorm) #x.pf[t-1,i]
        
        print(paste("t",t,"u",u,"x.old",x.old[i],"i",i,"a",a))
        if(a>1){
          print("alert. increase M")
        }
        u<-log(runif(1,0,1))
      }
      
      x.pf[t,i]<-z
    }
    
    
    #compute diagnostics
    m[t]<-mean(x.pf[t,])
    MSE.k[t]<-(x[t]-m[t])^2
    conf[t,]<-quantile(x.pf[t,],probs=c(0.025,0.975))
    
    
  }
  
  #return output
  return(list(x.out=x,y.out=y, x.pf.out=x.pf, m.out=m,MSE.k.out=MSE.k, 
              conf.out=conf))
  
}




#--Extended filter--
#uses optimal proposal
# N- No of particles
# x- vector of states
# y- vector of measurements
# x0- starting value
# sigma- standard deviation of process 
# sigma.meas- standard deviation of measurment
# resample.type - type of resampling (multinomial,stratified,systematic,residual,standard)

extended.filter<-function(N,x,y,x0,sigma=1,sigma.meas=1,resample.type="multinomial",N.thr.per=0.5){
  
  #threshold N
  N.thr<-as.integer(N.thr.per*N)
  
  #set up filter
  x0.pf<-rnorm(N,mean=x0,sd=sigma)
  x.pf <-matrix(NA,ncol=N,nrow=tau)
  
  #set up diagnostics
  #-adds variables for diagnoistics
  diagnostic.init(tau)
  
  
  #local function to choose resampling method
  FUN.resample<-function(N,x,w){
    switch(resample.type,
           multinomial = mult.resample(N=N,x=x,w=w) ,
           stratified = strat.resample(N=N,x=x,w=w) ,
           systematic = sys.resample(N=N,x=x,w=w),
           residual= res.resample(N=N,x=x,w=w),
           standard=sample(x=x,size=N,replace=TRUE,prob=w)
    )
  }
  
  
  #initialise
  sigma.prop<-1/sqrt( 1/sigma^2+diff.g(f(x0.pf,1))^2/sigma.meas^2 )
  mu.prop<-sigma.prop^2*(f(x0.pf,1)/sigma^2+diff.g(f(x0.pf,1))/sigma.meas^2)*(y[1]-g(f(x0.pf,1))+diff.g(f(x0.pf,1))*f(x0.pf,1))
  x.pf[1, ] <- f(mu.prop,1)+rnorm(N,sd=sigma.prop) 
  w.tilde<- dnorm(y[1],mean=g.tay(x.pf[1,],x0.pf,1),sd=sqrt(sigma^2+diff.g(f(x0.pf,1))^2*sigma.meas^2)) #does not use linearise
  w<-w.tilde/sum(w.tilde)
  
  
  #compute diagnostics
  #- mean, effective sample size, confidence intervals 
  m[1]<-sum(w*x.pf[1,]) 
  MSE.k[1]<-(x[1]-m[1])^2
  N.eff[1]<-1/sum(w^2) 
  conf[1,]<-quantile(x.pf[1,],probs=c(0.025,0.975))
  
  #second step--
  
  for (t in 2:(tau)) {
    
    #initialise
    #sigma.prop<-1/sqrt(1/sigma^2+1/sigma.meas^2)
    #mu.prop<-sigma.prop^2*(x.pf[t-1,]/sigma^2+y[t]/sigma.meas^2)
    sigma.prop<-1/sqrt( 1/sigma^2+diff.g(f(x.pf[t-1,],t))^2/sigma.meas^2 )
    mu.prop<-sigma.prop^2*(f(x.pf[t-1,],t)/sigma^2+diff.g(f(x.pf[t-1,],t))/sigma.meas^2)*(y[t]-g(f(x.pf[t-1,],t))+diff.g(f(x.pf[t-1,],t))*f(x.pf[t-1,],t))
    
    #project state forward
    x.pf[t, ] <- f(mu.prop,t)+rnorm(N,sd=sigma.prop) 
    x.pf[1, ] <- f(mu.prop,1)+rnorm(N,sd=sigma.prop) 
    
    #compute weights
    #w.tilde<- dnorm(y[t],mean=g.tay(x.pf[t-1,],x.pf[t-2],t),sd=sqrt(sigma^2+sigma.meas^2))*w
    w.tilde<- dnorm(y[t],mean=g.tay(x.pf[t,],x.pf[t-1,],t),sd=sqrt(sigma^2+diff.g(f(x.pf[t-1,],1))^2*sigma.meas^2)) #does not use linearise
    
    w<-w.tilde/sum(w.tilde)
    
    
    #effective sample size
    N.eff[t]<-1/sum(w^2)
    
    #resample (bootstrap means resample everytime)
    if(N.eff[t]<N.thr){
      x.pf[t,]<-FUN.resample(N=N,x=x.pf[t,],w=w)
      w<-rep(1/N,N)
      res.point[t]<-1 #add "yes" for resampling point
      
    }
    
    #compute diagnostics
    m[t]<-sum(w*x.pf[t,])
    MSE.k[t]<-(x[t]-m[t])^2
    conf[t,]<-quantile(x.pf[t,],probs=c(0.025,0.975))
    
  }
  
  #return output
  return(list(x.out=x,y.out=y, x.pf.out=x.pf, m.out=m,MSE.k.out=MSE.k, 
              conf.out=conf, N.eff.out=N.eff,w.out=w,res.point.out=res.point))
  
}

