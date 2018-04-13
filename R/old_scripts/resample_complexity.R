#--COMPARE COMPLEXITY--
# generate random weights and resample


#set vector of times
time.mult<-rep(NA,10)
time.strat<-rep(NA,10)
time.sys<-rep(NA,10)
time.res<-rep(NA,10)

#loop through all times
for(i in 1:10){
  N<-i*500
  w<-runif(N,1,100)
  w<-w/sum(w)
  x<-1:N
  
  start<-proc.time()
  mult.resample(N,x,w)
  time<-proc.time()-start
  time.mult[i]<-time[3]
  
  start<-proc.time()
  strat.resample(N,x,w)
  time<-proc.time()-start
  time.strat[i]<-time[3]
  
  start<-proc.time()
  sys.resample(N,x,w)
  time<-proc.time()-start
  time.sys[i]<-time[3]
  
  start<-proc.time()
  res.resample(N,x,w)
  time<-proc.time()-start
  time.res[i]<-time[3]
}


#plot comparison
plot(500*seq(1,10),time.mult,type="l",main="Complexity Resampling Scheme",ylab="Elapsed Time(s)",xlab="No of Particles, N")
lines(500*seq(1,10),time.strat,col="red")
lines(500*seq(1,10),time.sys,col="blue")
lines(500*seq(1,10),time.res,col="green")
legend(1000, 30, legend=c("Multinomial", "Stratified","Systematic","Residual"),
       col=c( "black","red","blue","green"),lty=c(1,1,1,1), cex=0.8)




