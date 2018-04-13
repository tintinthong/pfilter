#--LIBRARY FOR MULTINOMIAL RESAMPLING--
# Resampling performed directly on particles(Not resampling indexes)
# This library used for comparisons of the schemes.


#--multinomial resampling--
#N- no of particles to resample
#x -vector of particles to resample
#w- normalised weights
#OUTPUT:N resampled particles
mult.resample<- function(N,x,w){
  
  if(N==0){      #residual resampling uses mult.resamples 
    x.out<-NULL  #occasionally, N input is 0
    return(x.out)
  }
  
  x.out<-rep(NA,N)
  u.list<-runif(N,0,1) #list of UNORDERED uniform random deviates
  Q<-cumsum(w) #create cdf
  
  #loop through the cdf for uniform deviate
  for(j in 1:N){
    i<-1
    test<-Q[i] #create cdf
    
    #sample particles corresponding to section of CDF deviate falls into
    while(u.list[j]>test){
      i<-i+1
      test<-Q[i] 
    }
    x.out[j]<-x[i]
  }
  return(x.out)
}


#stratified resampling
#N- no of particles to resample
#x -vector of particles to resample
#w - normalised weights
#OUTPUT:N resampled particles
strat.resample<-function(N,x,w){

  x.out<-rep(NA,N)
  u.list<-(seq(1,N)-1+runif(N,0,1))/N #generate ORDERED uniform random deviates
  Q<-cumsum(w)#create cdf
  
  #loop through the cdf for uniform deviate
  for(j in 1:N){
    i<-1
    test<-Q[i]
    #sample particles corresponding to section of CDF deviate falls into
    while(u.list[j]>test){
      i<-i+1
      test<-Q[i]
    }
    x.out[j]<-x[i]
  }
  return(x.out)
}

#systematic resampling
#N- no of particles to resample
#x -vector of particles to resample
#w - normalised weights
#OUTPUT:N resampled particles
sys.resample<-function(N,x,w){

  x.out<-rep(NA,N)
  u.list<-(seq(1,N)-1+runif(1,0,1))/N #generate EVENLY SPREAD ORDERED uniform deviates
  Q<-cumsum(w)#create cdf
  
  #loop through the cdf for uniform deviate
  for(j in 1:N){
    i<-1
    test<-Q[i]
    #sample particles corresponding to section of CDF deviate falls into
    while(u.list[j]>test){
      i<-i+1
      test<-Q[i]
    }
    x.out[j]<-x[i]
  }
  return(x.out)
}

#--residual resampling--
#N- no of particles to resample
#x -vector of particles to resample
#w - normalised weights
#OUTPUT:N resampled particles
res.resample<- function(N,x,w){
  
  x.out<-rep(NA,N)
  
  #include expected number of particles
  ni<-as.integer(N*w)
  ni.pos<-which(ni!=0) #find positions of non-zero ni
  x.out<-rep(x[ni.pos],ni[ni.pos])
  m<-N-sum(ni)
  
  
  #re-compute weights
  w<-N*w-ni
  w<-w/sum(w) 
  
  x.out<-c(x.out,mult.resample(m,x ,w))
  return(x.out)
  
}




