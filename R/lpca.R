
#' Perform local PCA on each cluster of the data
#'
#' Takes in indicator of cluter membership, data, and intrinsic dimension
#' @param indicatr cluster membership or number of clusters for pam() clustering
#'  if a single positive interger is provided
#' @param X Data.
#' @param d intrinsic dimension
#' @return a list of representation of data (X.rep), mean normalized reconstruction error (mean_error),
#' normalized reconstruction error for all data (all_error), cluster membership (cluster_id),
#' mean normalized reconstruction error in each cluster (each_error), cluster size (cluster_size),
#' variance explained by each PC in each cluster (variance_proportion),
#' d/number of PCs needed to explaine more than d of the variance in each cluster (num_ev).
#'
#' @examples
#'
#' ############################## example I: Open box
#' ## package for 3d plot
#' library(rgl)
#' ## package for pam() kmeans clustering
#' library(cluster)
#'
#' ## load data
#' data(OpenBox)
#'
#' ## intrinsic dimension of the data
#' trueDim = 2
#'
#' ## number of clusters
#' K = 6
#'
#' indi = pam(OpenBox,K)$clustering
#'
#' temp = lpca(indi,OpenBox,trueDim)
#' OpenBox_rep = temp[[1]]
#' error_rep = temp[[2]]
#'
#' open3d()
#' plot3d(OpenBox,col=indi,xlim=c(0,2),ylim=c(0,2),zlim=c(0,2))
#' open3d()
#' plot3d(OpenBox_rep,col=indi,xlim=c(0,2),ylim=c(0,2),zlim=c(0,2))
#'
#'
#'
#' ############################## example II: Swiss roll
#' ## package for 3d plot
#' library(rgl)
#' ## package for pam() kmeans clustering
#' library(cluster)
#'
#' ## load data
#' data(SwissRoll)
#'
#' ## intrinsic dimension of the data
#' trueDim = 2
#'
#' ## number of clusters
#' K = 8
#'
#' indi = pam(SwissRoll,K)$clustering
#'
#' temp = lpca(indi,SwissRoll,trueDim)
#' SwissRoll_rep = temp[[1]]
#' error_rep = temp[[2]]
#'
#' open3d()
#' plot3d(SwissRoll,col=indi,xlim=c(0,2),ylim=c(0,2),zlim=c(0,2))
#' open3d()
#' plot3d(SwissRoll_rep,col=indi,xlim=c(0,2),ylim=c(0,2),zlim=c(0,2))
#'
#'
#' ############################## example III: M-shape
#' ## package for pam() kmeans clustering
#' library(cluster)
#'
#' ## load data
#' data(M_shape)
#'
#' ## intrinsic dimension of the data
#' trueDim = 1
#'
#' ## number of clusters
#' K = 4
#'
#' indi = pam(M_shape,K)$clustering
#'
#' temp = lpca(indi,M_shape,trueDim)
#' M_shape_rep = temp[[1]]
#' error_rep = temp[[2]]
#'
#' indi_true = rep(1:4,each=nrow(M_shape)/4)
#' temp_true = lpca(indi_true,M_shape,trueDim)
#' M_shape_rep_true = temp_true[[1]]
#'
#' par(mfrow=c(1,3))
#' plot(M_shape,col=indi)
#' plot(M_shape_rep,col=indi)
#' plot(M_shape_rep_true,col=indi_true)
#' @export
lpca=function(indicator,X,d){

  ## if length of indicator is one it is the num of clusters o.w. it is the indicator of the cluster membership
  ## if d <= 1 it is chosen as the # of e.vecs to explain var.prop=d of total variance, o.w. its the number of pc used in each cluster
  if(length(indicator)==1){
    k = indicator  ## k is the number of clusters
    indi = pam(X,k)$clustering
  } else {
    k = length(unique(indicator))
    indi = indicator
  }
  nc = ncol(X)  ##  number of grids

  ev.all = matrix(0,k,nc)  ##  store eigenvalues for the PCA for each cluster

  rep.error = 0

  X.rep = matrix(0,nrow(X),nc)  ##  The representation of the data using LPCA

  num.e = rep(0,k)		##  number of e-functions needed to explain 0.95 of variance in each cluster

  cluster.size=numeric(k)
  for(i in 1:k){         #####LPCA for each cluster

      current = matrix(X[indi==i,],ncol=nc)
      j = nrow(current)
      cluster.size[i]=j

      if(j>1){        #####if more than one curve are in this cluster, then do PCA

        current.mean = apply(current,2,mean)
        current.center=current-matrix(current.mean,nrow=nrow(current), ncol=ncol(current), byrow=TRUE)
        current.cov= t(current.center)%*%current.center/j
        temp=eigen(current.cov,symmetric=TRUE)
        current.eigenvec = temp$vectors   ## eigenfunctions
        current.eigenval = temp$values  ## eigenvalues
        current.pcscore = current.center%*%current.eigenvec

        ## percent of cumulative variance explained by eigenvalues in this cluster
        temp.var=cumsum(current.eigenval)/sum(current.eigenval)

        if(d<1){
          num.eigen = which(temp.var>=d)[1]	#####number of eigen vecs to use in the representation
        } else {
          num.eigen = d
        }
        num.e[i] = num.eigen

        ##projected
        X.rep[indi==i,] = matrix(current.mean,nrow=nrow(current), ncol=ncol(current), byrow=TRUE)+matrix(current.pcscore[,1:num.eigen], ncol=num.eigen)%*%t(matrix(current.eigenvec[,1:num.eigen], ncol=num.eigen))
        ev.all[i,1:length(current.eigenval)] = current.eigenval

      }else{##only one curve

        X.rep[indi==i,] =current
        num.e[i] = 1

      }
  }##end for

  ## average nomarlized reconstruction error in each cluster
  rep.e.each = rep(0,k)

  ## nomarlized reconstruction error for all data points
  error.all = apply((X-X.rep)^2,1,mean)/apply(X,1,var)

  ## average nomarlized reconstruction error for all data points
  rep.error = mean(error.all)

  for(i in 1:k){
    if(sum(indi==i)>1){
      rep.e.each[i] = mean(error.all[indi==i])
    } else {
      rep.e.each[i] = error.all[i]
    }
  }

  ev.all[cluster.size==1,1]=1   ## for cluster with only one curve, the curve itself explains all
  total= apply(matrix(ev.all,ncol=ncol(ev.all)),1,sum)
  variance.prop = ev.all/total

  # output = list(variance=variance.prop,rep.error=rep.error,error.all=error.all,cluster.size=cluster.size, X.rep=X.rep, num.e=num.e, K=mean(num.e),clustering=indi, rep.e.each=rep.e.each) ##, indi.theta=indi.theta, indi.alpha=indi.alpha)
  output = list(X.rep=X.rep,mean_error=rep.error,all_error=error.all,cluster_id=indi,each_error=rep.e.each,cluster_size=cluster.size,variance_proportion=variance.prop,num_ev=num.e)
  class(output) = "lpca"
  output
}
