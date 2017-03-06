#' Local principal component analysis (LPCA).
#' Extends k-means to k-planes.
#' In k-means centers are singlton points, but in k-planes centers are the
#' d-dimensional hyperplanes where d is the intrinsic dimension of the data.
#'
#' Takes in data, number of clusters, intrinsic dimension,
#' maximum iteration number and threshold for stopping
#' @param X Data.
#' @param k Number of clusters.
#' @param d Intrinsic dimension.
#' @param iter.max Maximum number of iterations.
#' @param thresh Threshold for stopping, percentage change of normalized reconstruction error.
#' @return a list of cluster membership (indicator), representaion (X.rep), average normalized
#' reconstruction error (error).
#'
#' @examples
#' ############################## example I: Open box
#' ## package for 3d plot
#' library(rgl)
#'
#' ## load data
#' data(OpenBox)
#'
#' ## intrinsic dimension of the data
#' trueDim = 2

#' ## threshold for kplanes stopping (relative reconstruction error change ratio)
#' thresh = 0

#' ## number of clusters
#' K = 6

#' set.seed(0)
#' temp = kplanes(OpenBox, k=K, d=trueDim, thresh=thresh)
#' OpenBox_rep = temp[[2]]
#' cluster_id = temp[[1]]
#' error_rep = temp[[3]]

#' open3d()
#' plot3d(OpenBox,col=cluster_id,xlim=c(0,2),ylim=c(0,2),zlim=c(0,2))
#'
#'
#' ############################## example II: Swiss roll
#' ## package for 3d plot
#' library(rgl)
#'
#' ## load data
#' data(SwissRoll)
#'
#' ## intrinsic dimension of the data
#' trueDim = 2
#'
#' ## threshold for kplanes stopping (relative reconstruction error change ratio)
#' thresh = 0
#'
#' ## number of clusters
#' K = 8
#'
#' set.seed(0)
#' temp = kplanes(SwissRoll, k=K, d=trueDim, iter.max=1000, thresh=thresh)
#' SwissRoll_rep = temp[[2]]
#' cluster_id = temp[[1]]
#' error_rep = temp[[3]]
#'
#' open3d()
#' plot3d(SwissRoll,col=cluster_id,xlim=c(0,2),ylim=c(0,2),zlim=c(0,2))
#'
#' @export
kplanes = function(X,k,d,iter.max=50,thresh=0){
  ## extended kmeans with the centroid not being a point any more but a d-dimenstonal linear space
  n = nrow(X)
  N = ncol(X)
  step = 1

  ## initialize k clusters
  indicator = kmeans(X,k)$cluster

  ## store cluster membership in all iterations
  indi.all = matrix(0,iter.max+1,n)

  ## representation of all data points in all iterations
  rep.all.all = array(0,dim=c(iter.max+1,n,N))

  indi.all[1,]=indicator

  ## record normalized reconstruction error of each point to its
  ## representation in all clusters using the charateristic hyperplane
  ## in that cluster for all data points. This is also the projection distance
  ## of each data point to they cluster-wise hyperplane.
  error.all = matrix(0,k,n)

  ## initialize error and relative change in normalized reconstruction error
  error.now = 1e+5
  error.old = 1e-5
  e.ratio = abs((error.now-error.old)/error.old)

  ## store representation of each cluster in current iteration
  rep.all = array(0,dim=c(k,n,N))

  while(step<iter.max&e.ratio>thresh){
    print(step)

    evec.all = array(0,dim=c(k,N,N))
    eval.all = matrix(0,k,N)
    pc.all = array(0,dim=c(k,n,N))
    #data = X[which(indicator==i),]
    num.e = rep(0,k)		#####number of e-functions needed to explain 0.95 of variae in each cluster

    ## store sample size of each cluster
    cluster.size=numeric(k)

    for(i in 1:k){         #####LPCA for each cluster
      if(sum(indicator==i)>0){
        current = matrix(X[indicator==i,],ncol=N)
        j = nrow(current)
        cluster.size[i]=j

        if(j>1){        #####if more than one curve are in this cluster, then do PCA
          current.mean = apply(current,2,mean)
          current.center=current-matrix(current.mean,nrow=nrow(current), ncol=ncol(current), byrow=TRUE)
          current.cov= t(current.center)%*%current.center/j
          temp=eigen(current.cov,symmetric=TRUE)
          current.eigenvec = temp$vectors      ## eigenfunctions
          current.eigenval = temp$values  ## eigenvalues
          current.pcscore = current.center%*%current.eigenvec

          ## record the range pcscores (not used here, might be useful later to construct
          ## convex connected cluster (needs more investigation))
          upperB = apply(current.pcscore,2,max)
          lowerB = apply(current.pcscore,2,min)

          temp.var=cumsum(current.eigenval)/sum(current.eigenval)

          evec.all[i,,]=current.eigenvec
          eval.all[i,]=current.eigenval

          ## compute the pcscores of all data points using the eigenvectors from this cluster
          all.center = X-matrix(current.mean,nrow=nrow(X), ncol=ncol(X), byrow=TRUE)
          all.pcscore = all.center%*%current.eigenvec

          ## used to bound pc scores in the existing segment
          # for(t in 1:n){
          #   for(c in 1:N){
          #     if(all.pcscore[t,c]>upperB[c]){
          #       all.pcscore[t,c]=upperB[c]
          #     } else if(all.pcscore[t,c]<lowerB[c]){
          #       all.pcscore[t,c]=lowerB[c]
          #     }
          #   }
          # }

          pc.all[i,,]=all.pcscore
          num.eigen = d
          num.e[i] = which(cumsum(current.eigenval)/sum(current.eigenval)>0.95)[1]
          X.rep.now = matrix(current.mean,nrow=nrow(X), ncol=ncol(X), byrow=TRUE)+matrix(all.pcscore[,1:num.eigen], ncol=num.eigen)%*%t(matrix(current.eigenvec[,1:num.eigen], ncol=num.eigen))
          rep.all[i,,]=X.rep.now
          #index.out = which(apply(((all.pcscore>upperB)+(all.pcscore<lowerB)),1,max)>0)
          #print(dim(index.out))
          #print(dim(all.pcscore))
          #print(length(index.out))
          #print(index.out)
          #print(apply(((all.pcscore>upperB)+(all.pcscore<lowerB)),1,max)>0)
          #print(index.out)
          error.all[i,] = apply((X-X.rep.now)^2,1,mean)/apply(X,1,var)
          #error.all[i,index.out] = 1e+10
          #ev.all[i,1:length(current.eigenval)] = current.eigenval
        }else{##only one curve
          X.rep.now = X
          rep.all[i,,]=X.rep.now
          error.all[i,] = apply((X-X.rep.now)^2,1,mean)/apply(X,1,var)
          num.e[i] = 1
        }
      }##end if
    }##end for

    X.rep = matrix(0,n,N)
    for(s in 1:n){
      ## updating cluster membership
      indicator[s] = which(error.all[,s]==min(error.all[,s]))[1]

      ## updating representation
      X.rep[s,] = rep.all[indicator[s],s,]
    }

    print(error.old)

    error.vector = apply((X-X.rep)^2,1,mean)/apply(X,1,var)
    error.now = mean(error.vector)
    e.ratio = abs(error.now-error.old)/error.old
    error.old = error.now

    step=step+1
    indi.all[step,]=indicator
    rep.all.all[step,,]=X.rep
    print(error.now)
    print(e.ratio)
  }

  # output = list(indicator=indicator,X.rep=X.rep,error=error.now,indi.all=indi.all,error.all=error.all,rep.all.all=rep.all.all)
  output = list(indicator=indicator,X.rep=X.rep,error=error.now)
  output
}
