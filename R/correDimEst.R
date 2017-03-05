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
#' Function storing the logC2(epsi) values and the numerical derivative of logC2 w.r.t log(epsi).
#' @param num Number of epsilons.
#' @param distance Distance object (upper triangle) from dist() fn.
#' @return a list of log(epsilons) (x), y and derivatives of y.
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
