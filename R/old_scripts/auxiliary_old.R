#auxiliary particle filter
#-sampling from q(xt,i|z1:t)
#-use importance density as prior
#process model N(xt,sd=1)
#measurement model N(yt,xt,sd=1)

N<-500
tau<-100
N.thr<-10

set.seed(111)
x<-cumsum(rnorm(tau)) #true states
y<-x+rnorm(tau) #measurements

x.pf<-matrix(rep(NA,N*(tau+1)),nrow=tau+1)
w.pf<-matrix(rep(NA,N*(tau+1)),nrow=tau+1)
v.pf<-matrix(rep(NA,N*(tau+1)),nrow=tau+1)
w.var<-rep(NA,tau+1)

x.pf[1,]<-rnorm(N)
w.pf[1,]<-1/N

for(t in 2:(tau+1)){
  
  #make f model(mean of projection model)
  mu<-mean(x.pf[t-1,])
  #create parent weights
  v.pf[t,]<-w.pf[t-1,]*dnorm(y[t-1],mu)
  #normalise parent weights
  v.pf[t,]<-v.pf[t,]/ sum(v.pf[t,])
  
  #sample parent trajectories
  samp<-sample(1:N,size=N,replace=TRUE,prob=v.pf[t,])
  
  x.pf[t,]<-rnorm(N,x.pf[t-1,samp]) #prior
  
  w.tilde<-dnorm(y[t-1],x.pf[t,])/dnorm(y[t-1],mu)
  #(dnorm(y[t-1],x.pf[t,])*dnorm(x.pf[t,],x.pf[t-1,])*w.pf[t-1,samp])/(dnorm(x.pf[t,],x.pf[t-1,samp])*v.pf[t,samp])
  #v.pf[t,samp] is questionable
  
  #normalise weights
  w.pf[t,]<-w.tilde/sum(w.tilde)
  
  #compute variance of importance weights
  w.var[t]<-var(w.pf[t,])

  
}


#diagnostic plots
plot(x)
points(y,col="red")

plot(x)
lines(m,col="red")

plot(MSE)

plot(N.eff)

