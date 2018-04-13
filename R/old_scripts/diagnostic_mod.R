#--LIBRARY FOR DIAGNOSTICS--
# initializes variables that involve diagnostics

#initial diagnostics
# create variables for diagnostics in calling environment(in side particle filter)
#tau - number of time step
diagnostic.init<-function(tau){
  #mean estimate 
  assign("m",rep(NA,10),pos=sys.frame(-1))
  #effective sample size
  assign("N.eff",rep(NA,10),pos=sys.frame(-1))
  #squared difference between state and mean estimate
  assign("MSE.k",rep(NA,10),pos=sys.frame(-1)) 
  #credbile intervals
  assign("conf",matrix(NA,ncol=2,nrow=tau),pos=sys.frame(-1)) 
  #resampling points (1's for resampled and 0's  for non-resample)
  assign("res.point",rep(0,tau),pos=sys.frame(-1)) 
  #create vector for weighted variance
  assign("var.s",rep(NA,tau),pos=sys.frame(-1)) 
}

#--count time between resamples and no of resamples(Not Used)--
# assume first step is a resample
#OUTPUT: list of 1) vector of time between resamples 2) total no of resamples
count.resamples<-function(res.point){
  M<-length(res.point)
  tt<-c() #time between resamples
  count<-1 #count for NOT resampling
  for(i in 2:M){
    if(res.point[i]==0){
      count<-count+1
    }else{
      if(count!=0){
        tt<-c(tt,count)
      }
      count<-0
    }
    if(i==M){
      if(res.point[i]==0){
        tt<-c(tt,count)
      }
    }
  }
  return(list(tt.out=tt,resamples.out=sum(res.point)))
}


