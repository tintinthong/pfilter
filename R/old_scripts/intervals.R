#--COMPUTE CONFIDENCE AND CREDIBLE INTERVALS--

#compare particle filter and particle smooth 
#-- particle smooths get too overconfident and so neglect its 
#--using credible intervals

obj.prior<-particle.filter(N=N,x=x,y=y.mat[,1],x0=5,sigma=sigma,sigma.meas=sigma.meas,resample.type="standard",N.thr.per=1)
x.pf<-obj.prior$x.pf.out
conf.mat<-apply(x.pf,1,quantile,probs=c(0.025,0.975))
conf.mat.main<-matrix(NA,tau,2)
count<-0
for(i in 1:tau){
  confu<-conf.mat[2,i]
  confl<-conf.mat[1,i]
  if(x[i]<confl || x[i]>confu){
    count<-count+1 #this is not really the right way to count
  }
  conf.mat.main[i,]<-c(confl,confu)
  
}


#----------------------------------------------------------------------------------------------------------------------------
#frequentist confidence interval

count<-0
for(j in 1:100){
  obj.prior<-particle.filter(N=N,x=x,y=y.mat[,1],x0=5,sigma=sigma,sigma.meas=sigma.meas,resample.type="standard",N.thr.per=1)
  x.pf<-obj.prior$x.pf.out
  means.prior<-rowMeans(x.pf)
  var.prior<-apply(x.pf,1,var)
  conf.mat<-matrix(NA,tau,2)
  
  for(i in 1:tau){
    
    #confu<-means.prior[i]+1.96*sqrt(var.prior[i]/N)
    #confl<-means.prior[i]-1.96*sqrt(var.prior[i]/N)
    confu<-obj.prior$m.out[i]+1.96*sqrt(obj.prior$var.s.out[i]/N) #this cant be used because its after resampling
    confl<-obj.prior$m.out[i]-1.96*sqrt(obj.prior$var.s.out[i]/N)
    
    conf.mat[i,]<-c(confl,confu)
    
  }
  if(x[50]<conf.mat[50,1] || x[50]>conf.mat[50,2]){
    count<-count+1
  }
  
}

plot(x)
lines(conf.mat[,1],col="red")
lines(conf.mat[,2],col="red")
count
