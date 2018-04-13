#--LIBRARY TO GENERATE STATES AND MEASUREMENTS--
# tau is exclusive of x0
# Sigma is the covariannce so it's diagonals are variance, sigma is standard deviation
# X is a matrix, x is a vector, xi is scalar, x0 is scalar, X0 is vector

#load library
# - only used for high-dimensional case
library(mvtnorm)


#-----------------------------------Random Walk Example----------------------------

#--Generate random walk--
#tau- no of states
#x0-known starting value 
#sigma- fixed standard deviation of process
#OUTPUT: vector of states of length tau
rand.walk.1D<-function(tau=100,x0=10,sigma=1){
  x<-cumsum(c(x0,rnorm(tau,sd=sigma))) 
  return(x[2:(tau+1)])
}

#--Generate measurements from given states--
#x- vector of states
#sigma.meas- standard deviation
#OUTPUT:vector of measurements 
rand.y.1D<-function(x,sigma.meas=1){
  y<-rnorm(length(x),mean=x,sd=sigma.meas)
  return(y)
}
  
#-----------------------------------Non-linear Example----------------------------

#--non-linear function--
# f in process model
#xj- state from previous time
#t- current time
#OUTPUT:deterministic value of function f
f<-function(xj,t){
  return(xj/2+25*xj/(1+xj^2) +8*cos(1.2*t))
}

#--non-linear function--
# g in measurement model
#xi- state from previous time
#OUTPUT:deterministic value of function g
g<-function(xi){
  return(xi^2/20)
} 


#--generate population process--
#tau- no of states 
#x0-known starting value
#sigma- standard deviation of process
#OUTPUT:vector of tau states
process.pop<-function(tau=100,x0=10,sigma=sqrt(10)){
  x<-rep(NA,tau)
  x[1]<-f(x0,1)+rnorm(1,mean=0,sd=sigma)
  for(i in 2:(tau+1)){
    x[i]<-f(x[i-1],i)+rnorm(1,mean=0,sd=sigma)
  }
  return(x[1:tau])
}
#REMEMBER
#\\f is not passed in as a function, but process.1D searches for it in it's parent environment(global)
#\\DANGEROUS but convenient


#--generate  measurements from population process--
#x- vector of states 
#sigma.meas- standard deviation of measurement model
#OUTPUT:vector of measurements 
meas.pop<-function(x,sigma.meas=1){
  tau<-length(x)
  y<- x^2/20+ rnorm(tau,mean=0,sd=sigma.meas)
}



#-------------------------------Functions which are Not Used-------------------------
#May be useful for further research


#--Generate multi-dimensionlal random walk--
#tau- no of states 
#X0- vector of starting values
#Sigma- variance matrix
#OUTPUT: matrix of states 
rand.walk<-function(tau=100,X0=c(0,10),Sigma=diag(c(1,1))){
  d<-length(X0) 
  X<-rmvnorm(tau, mean = rep(0,d), sigma = Sigma) 
  X<-rbind(X0,X)
  X<-apply(X, 2, cumsum)
  return(X[2:(tau+1),])
}


#--generate data(using Normal Process)--
#X- matrix of states(rows for time+ col for dimension)
#Sigma.meas- variance matrix 
#OUTPUT: matrix of measurements
rand.y<-function(X,Sigma.meas=diag(c(1,1))){
  rowX<-nrow(X)
  colX<-ncol(X)
  Y<-matrix(NA,nrow=rowX,ncol=colX) 
  for( i in 1:rowX){
    Y[i,]<-rmvnorm(1,mean=X[i,],sigma= Sigma.meas)
  }
  return(Y)
}


#--Linearise measurement for population measurement --
# used for extended particle filter
# taylor expanded measurement model around f(xj)
#x0- known starting value
#x- vector of states
#sigma.meas- standard deviation of measurement model
#OUTPUTE: vector of measurements
extended.meas.1D<-function(x0,x,sigma.meas){
  tau<-length(x)
  y<-rep(NA,1)
  y[1]<-g(f(x0,1))+diff.g(f(x0,1))*(x0-f(x0,1))+rnorm(1,sigma.meas)
  for(i in 2:tau){
    y[i]<-g(f(x[i-1],i))+diff.g(f(x[i-1],i))*(x[i]-f(x[i-1],i))+rnorm(1,sigma.meas)
  }
  return(y)
}  
#REMEMBER
#\\f,g,diff.g is not passed in as a function, but extended.meas.1D searches for it in it's parent environment(global)
#\\DANGEROUS but convenient  


#-- first derivative of g wrt to xi--
# will be used for taylor expansion in extended pf
#xi- state at current time
#OUTPUT: differential value at xi
diff.g<-function(xi){
  return(xi/10)
}

#--taylor expansion of measurement model--
#xi- state at current time
#xj -state at previous time
g.tay<-function(xi,xj,i){
  return(
    g(f(xj,i))+diff.g(f(xj,i))*(xi-f(xj,i))
  )
}

