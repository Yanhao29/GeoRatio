#' Calculating proportion of the pairwise L2-distance that is less than epsilon
#' @param epsilon Epsilon value.
#' @param distance Distance object.
#' @return Proportion of the pairwise L2-distance that is less than epsilon.
C2 = function(epsilon,distance){
  n=length(distance)
  return(sum(distance<=epsilon)/n)
}


#' Correlation dimension estimation
#'
#' Takes in pair-wise distance of data and number of epsilons to estimate intrinsic dimension of data.
#' Function returns, log(epsi) (stored in x), log(C2) (stored in y) values and the numerical derivative of log(C2) w.r.t log(epsi) (stored in deri). The estimated dimension is the level of the plateau (stable part) of the deriative of y w.r.t. x.
#' @param num Number of epsilons.
#' @param distance Distance object (upper triangle) from dist() fn.
#' @return a list of log(epsi) (x), lob(C2) (y) and derivatives of y.
#'
#' @examples
#' ############################## example I: Spiral data
#' ## load data
#' data(spiral)
#'
#' ## intrinsic dimension of the data
#' trueDim = 1
#'
#' ## number of epsilons
#' num = 50
#'
#' ## pair-wise Euclidean distance
#' dis = dist(spiral,method = "euclidean", diag = FALSE, upper = TRUE)
#'
#' ## correlation dimension estimation
#' est = correDimEst(num,dis)
#'
#' ## plot correlation dimension estimation
#' par(mfrow=c(1,2))
#' plot(spiral,xlab='',ylab='')
#' # plot(est$x,est$y,xlab="log(epsilon)",ylab="log(C2(epsilon))",main="log-log plot")
#' plot(est$x[2:(length(est$x)-1)],est$deri,type="l",ylim=c(0,10),xlab="log(epsilon)",ylab="est. correlation dimension",
#'      main="Estimated correlation dimension")
#' abline(trueDim,0,lty=3,col="gray60",lwd=2)
#'
#'
#' ############################## example II: Open box data
#' library(rgl) ## for open3d(), 3d plot
#'
#' ## load data
#' data(OpenBox)
#'
#' ## intrinsic dimension of the data
#' trueDim = 2
#'
#' ## number of epsilons
#' num = 50
#'
#' ## pair-wise Euclidean distance
#' dis = dist(OpenBox,method = "euclidean", diag = FALSE, upper = TRUE)
#'
#' ## correlation dimension estimation
#' est = correDimEst(num,dis)
#'
#' ## plot correlation dimension estimation
#' open3d()
#' plot3d(OpenBox,col=col.f,axes=FALSE,box=FALSE,xlab = "", ylab = "", zlab = "",xlim=c(0,2),ylim=c(0,2),zlim=c(0,2),cex.lab=1)
#'
#' par(mfrow=c(1,1))
#' # plot(est$x,est$y,xlab="log(epsilon)",ylab="log(C2(epsilon))",main="log-log plot")
#' plot(est$x[2:(length(est$x)-1)],est$deri,type="l",ylim=c(0,10),xlab="log(epsilon)",ylab="est. correlation dimension",
#'      main="Estimated correlation dimension")
#' abline(trueDim,0,lty=3,col="gray60",lwd=2)
#'
#' @export
correDimEst = function(num,distance){

  ## minimum of log(epsilon)
  min = log(distance[order(distance)][2])

  ## maximum of log(epsilon)
  max = log(distance[order(-distance)][2])

  ## create log equally spaced sequence of epsilons
  grid = seq(min,max,length=num)
  epsilon = exp(grid)

  y = c(rep(0,num))
  for(i in 1:num){
    y[i] = log(C2(epsilon[i],distance))
  }

  deri = c(rep(0,num-2))
  delta = grid[2]-grid[1]
  for(i in 1:(num-2)){
    deri[i] = (y[i+2]-y[i])/(2*delta)
  }
  list(y=y,x=grid,deri=deri)            #####store log(epsi) (x), log(C2) (y), and derivative of y
}
