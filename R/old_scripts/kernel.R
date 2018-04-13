# Given a kernel density estimate, this function 
# carries out a (very quick and dirty) numerical 
# integration, and then fits a spline to get a function 
# which can be used to look up cumulative probabilities. 
smoothed.df <- 
  function(d) 
  { 
    F <- cumsum(d$y) 
    F <- F/F[length(F)] 
    splinefun(d$x, F) 
  } 
# Generate a bimodal test distribution 
# Estimate the desnsity and distribution function 
x <- rnorm(1000) + ifelse(runif(1000) > .5, -3, 3) 
d <- density(x) 
F <- smoothed.df(d) # F returns cumulative probs 
# Plot the true and estimated distribution function 
curve(0.5 * dnorm(x, -3) + 0.5 * dnorm(x, 3), col="red") 
lines(d) 
# Plot the true and estimated distribution function 
F <- smoothed.d
curve(F(x), add=TRUE) 


den.epan<-function(u,h,x0){
  return((1/h)*(3/4*(1-((u-x0)/h)^2)))
}

den.epan<-function(u,x0,h){
  return((1/h)*(3/4*(1-((u-x0)/h)^2)))
}

den.epan((-1+10)*0.1,10,0.1)



#inverse cdf of epnachnikov kernel
z<-0
2*sin(13*asin(2*z-1))




den.epan(-1+10,1,10)

u<-runif(100/0.1,-1,1)
plot(sort(u),den.epan(sort(u),0.1,0))

cdf.epan<-function(u,h){
  return((1/h)*(3/4)*((u/h)- (u/h)^3/3))
}