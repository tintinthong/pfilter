#load library
source("filter")

#set model variables
x0<-0
sigma<-1
sigma.meas<-1
tau<-100

#generate states and measurements
set.seed(123)
x<-rand.walk.1D(tau=tau,x0=x0,sigma=sigma)
y<-rand.y.1D(x,sigma.y=sigma.meas)

#particle filter function 
particle.filter<-function(N,x,y,x0,sigma=1,sigma.meas=1,resample.type="multinomial",N.thr.per=0.5){
  #load functions
  source("resample.R")
  source("random-walk.R")
  
  #threshold N
  N.thr<-as.integer(N.thr.per*N)
  
  #set up filter
  x0.pf<-rnorm(N,mean=x0,sd=sigma)
  x.pf <- matrix(NA,ncol=N,nrow=tau)
  
  #set up diagnostics
  #- adds global variables for diagnostic(except for MSE)
  diagnostic.init(tau)
  count<-0 #count number of resamples
  tt<-1 #counter for time between resamples
  tt.between<-c() #list to compute time steps between resamples
  res.point<-rep(0,tau) #resampling points (1's for resampled and 0's  for non-resample)
  
  #local function to choose resampling method
  FUN.resample<-function(N,x,w){
    switch(resample.type,
           multinomial = mult.resample(N=N,x=x,w=w) ,
           stratified = strat.resample(x=x,w=w) ,
           systematic = sys.resample(x=x,w=w),
           residual= res.resample(N=N,x=x,w=w))
  }
  
  # Initialize
  x.pf[1, ] <- x0.pf+rnorm(N,sd=sigma) #project particles from x0 to the first time step
  w.tilde<- dnorm(y[1],mean=x.pf[1,],sd=sigma.meas) #the previous weight would be equal weights 1/N depending on how x0 was chosen
  w<-w.tilde/sum(w.tilde) #normalise weights
  x.pf[1,]<-FUN.resample(N=N,x=x.pf[1,],w=w) #resample
  w<-rep(1/N,N) #set weights to equal
  
    #compute diagnostics
    #- mean, effective sample size, confidence intervals 
    m[1]<-sum(w*x.pf[1,]) #compute mean
    MSE.k[1]<-(x[1]-m[1])^2
    #MSE[1,k]<-(x[1]- mean(x.pf[1,]))^2
    N.eff[1]<-1/sum(w^2) #compute effective sample size
    conf[1,]<-quantile(x.pf[1,],probs=c(0.025,0.975)) #compute confidence intervals
    
    for (t in 2:(tau)) {
      
      #project  state forward
      x.pf[t, ] <-x.pf[t-1,]+rnorm(N,sd=sigma)
      
      #compute weights
      w.tilde <- dnorm(y[t], mean=x.pf[t, ],sd=sigma.meas)*w #when add w resampling does not help
      w <- w.tilde/sum(w.tilde) #normalise weights
      
      
      #effective sample size
      N.eff[t]<-1/sum(w^2)
      
      #resample (bootstrap means resample everytime)
      if(N.eff[t]<N.thr){
        x.pf[t,]<-FUN.resample(N=N,x=x.pf[t,],w=w)
        w<-rep(1/N,N)
        count<-count+1 #add to resample count
        res.point[t]<-1 #add "yes" for resampling point
        tt.between<-c(tt.between,tt) #add number of steps between two resamples
        tt<-0 #reset resample counts
      }
      
      #compute diagnostics
      m[t]<-sum(w*x.pf[t,])
      MSE.k[t]<-(x[t]-m[t])^2
      conf[t,]<-quantile(x.pf[t,],probs=c(0.025,0.975))
      
    }
    
    #return output
    return(list(x.out=x,y.out=y, x.pf.out=x.pf, m.out=m,MSE.k.out=MSE.k, conf.out=conf, N.eff.out=N.eff,w.out=w,count.out=count,res.point.out=res.point,tt.between.out=tt.between))
    
}


#compute MSE monte carlo
M<-1
tau<-100 #remember tau because have to incorporate in line below
MSE<-matrix(NA,tau,M)
for(k in 1:M){
  obj.mult<-particle.filter(N=500,x=x,y=y,x0=x0,sigma=sigma,sigma.meas=sigma.meas,resample.type="multinomial",N.thr.per=1)
  MSE[,k]<-obj.mult$MSE.k.out
}
MSE<-colSums(MSE)/M


#plot Effective sample size
plot(obj.mult$N.eff.out)

#plot MSE
plot(MSE)

#extract objects
x.pf<-obj.mult$x.pf.out
conf<-obj.mult$conf.out
m<-obj.mult$m.out

#check plot
plot(x)
points(y,col="red")
for(i in 1:N){
  lines(x.pf[,i])
}


#plot confidence intervals
library(ggplot2)
df <- data.frame(t =1:tau,mu=x,L=conf[,1],U=conf[,2])
ggplot(df, aes(x = t, y = mu)) +
  geom_point(size = 1) +
  geom_errorbar(aes(ymax = U, ymin = L))









