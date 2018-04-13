#' @title Plot object with class pframe_1d
#'
#' @description Plot method for object with class pframe_1d. The method prints the trend of each
#' particle trajectory.
#'
#' @param obj Object of pframe_1d class
#' @param clean Print a prettified version of plot
#'
#' @return Graphical plot of obj
#'
#' @author Justin Thong \email{justinthong93@gmail.com}
#'
#' @examples
#' plot(particle(pframe_1d()))
#' set.seed(1)
#' x<-c(1,2,3,4,4,4,5,4,5,4,4,4,4,4,4) #true means of state model
#' plot(particle(pframe_1d(N=100,smooth=FALSE,y=rnorm(length(x),mean=x))))
#'
#' @import ggplot2
#' @import reshape2

#'@rdname plot
#' @export

plot.pframe_1d <- function(obj,clean=FALSE,...) {


  #1 Plot is particle paths

  X<-obj$X

  if(nrow(X)>10 && clean){
    X<-X[1:10,]
  }

  df<-as.data.frame(X)
  df$time<-1:nrow(X)
  df_melted= reshape2::melt(df,id.vars='time')
  ggplot2::ggplot(df_melted, aes(x = time, y = value, color = variable, group = variable))+
    geom_line()+theme(legend.position="none")+
    xlab("Time (t)")+
    ylab("x-values")+
    labs(
      title = "Particle Paths",
      subtitle = paste("resample=", attr(obj,"resample"), "smooth=", attr(obj,"smooth"))
      )

  #2 Plot particle paths sequentially


  #3 Plot metrics
    #plot bootstrap confidence intervals
    #plot mean

}

#df = data.frame(cat = LETTERS[1:6], VAR1 = runif(6), VAR2 = runif(6), VAR3 = runif(6), VAR4 = runif(6))
#df_melted = melt(df, id.vars = 'cat')
#ggplot(df_melted, aes(x = variable, y = value)) + geom_line(aes(color = cat, group = cat))












