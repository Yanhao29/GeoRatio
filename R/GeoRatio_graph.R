#' Building GeoRatio Graph
#'
#' Takes in data, minimum nbhd size, maximum nbhd size, intrinsic dimension, variance threshold
#' and construct a graph for estimating GeoRatio measure
#' @param X Data.
#' @param k_max maximum neighborhood size
#' @param k_min minimum neighborhood size
#' @param d intrinsic dimension
#' @param thresh threshold of proportion of variance explained by d leading eigenvalues in each neighborhood
#' @param gmode graph edge choosing
#' @param distance if data passed in is a distance matrix
#' @return a list of GeoRatio dissimilarity, geodesic distance, nbhd adjacency matrix, weighted graph, nbhd size
#' @examples
#' ############################## example I: V-shape
#' ## package for pam() kmeans clustering
#' library(cluster)
#' library(igraph) ## used by GeoRatio_graph
#'
#' ## load data
#' data(V_shape)
#'
#' ## intrinsic dimension of the data
#' trueDim = 1
#'
#' ############## correlation dimension estimation
#' ## number of epsilons
#' num = 50
#'
#' ## pair-wise Euclidean distance
#' dis = dist(V_shape,method = "euclidean", diag = FALSE, upper = TRUE)
#'
#' ## correlation dimension estimation
#' est = correDimEst(num,dis)
#'
#' ## plot correlation dimension estimation
#' par(mfrow=c(1,1))
#' plot(est$x[2:(length(est$x)-1)],est$deri,type="l",ylim=c(0,10),xlab="log(epsilon)",ylab="est. correlation dimension",
#'      main="Estimated correlation dimension")
#' abline(trueDim,0,lty=3,col="gray60",lwd=2)
#'
#' ##### decide k_M (plateau log(epsi) from -4.8 to -2.8)
#' epsi_s = -4.5
#' epsi_e = -2.5
#' epsi_seq = seq(epsi_s,epsi_e,length=5)
#'
#' dis_L2 = as.matrix(dist(V_shape))
#'
#' for(i in 1:length(epsi_seq)){
#'   epsi = epsi_seq[i]
#'   dist_thre = exp(epsi)
#'   epsi_nbhd = apply(dis_L2,1,function(x) sum(x<dist_thre))
#'   avg_epsi_nbhd = mean(epsi_nbhd)
#'   print(epsi)
#'   print(avg_epsi_nbhd)
#' }
#'
#' ## parameters for GeoRatio_graph()
#' K = 2
#' k_M = 15
#' k_m = 3
#' pca_thre = 0.9
#'
#' set.seed(100)
#'
#' V_GR = GeoRatio_graph(V_shape,k_max=k_M,k_min=k_m,trueDim,thresh=pca_thre,gmode="max")
#'
#' dis_L2 = dist(V_shape)
#' dis_G = as.dist(V_GR$dis_G)
#' dis_GR = as.dist(V_GR$dis_GR)
#'
#' hc_L2 = hclust(dis_L2,method="average")
#' hc_G = hclust(dis_G,method="average")
#' hc_GR = hclust(dis_GR,method="average")
#'
#' indi_KML2 = pam(dis_L2,K,diss=T)$clustering
#' indi_KMG = pam(dis_G,K,diss=T)$clustering
#' indi_KMGR = pam(dis_GR,K,diss=T)$clustering
#' indi_HL2 = cutree(hc_L2,k=K)
#' indi_HG = cutree(hc_G,k=K)
#' indi_HGR = cutree(hc_GR,k=K)
#'
#' lpca_V = kplanes(V_shape,K,trueDim,iter.max=100,thresh=1e-5)
#' indi_lpca = lpca_V$indicator
#'
#' par(mfrow=c(3,3))
#' plot(V_shape,col=indi_HGR,main="HC-GR")
#' plot(V_shape,col=indi_HG,main="HC-G")
#' plot(V_shape,col=indi_HL2,main="HC-L2")
#' plot(V_shape,col=indi_KMGR,main="KMC-GR")
#' plot(V_shape,col=indi_KMG,main="KMC-G")
#' plot(V_shape,col=indi_KML2,main="KMC-L2")
#' plot(V_shape,col=indi_lpca,main="LPCA")
#'
#' ############################## example II: Swiss roll
#' ## package for 3d plot
#' library(rgl)
#' ## package for pam() kmeans clustering
#' library(cluster)
#'
#' library(igraph) ## used by GeoRatio_graph
#'
#' ## load data
#' data(SwissRoll)
#'
#' ## intrinsic dimension of the data
#' trueDim = 2
#'
#' ############## correlation dimension estimation
#' ## number of epsilons
#' num = 50
#'
#' ## pair-wise Euclidean distance
#' dis = dist(SwissRoll,method = "euclidean", diag = FALSE, upper = TRUE)
#'
#' ## correlation dimension estimation
#' est = correDimEst(num,dis)
#'
#' ## plot correlation dimension estimation
#' par(mfrow=c(1,1))
#' plot(est$x[2:(length(est$x)-1)],est$deri,type="l",ylim=c(0,10),xlab="log(epsilon)",ylab="est. correlation dimension",
#'      main="Estimated correlation dimension")
#' abline(trueDim,0,lty=3,col="gray60",lwd=2)
#'
#' ##### decide k_M (plateau log(epsi) from -4.8 to -2.8)
#' epsi_s = -1.2
#' epsi_e = 1.8
#' epsi_seq = seq(epsi_s,epsi_e,length=5)
#'
#' dis_L2 = as.matrix(dist(SwissRoll))
#'
#' for(i in 1:length(epsi_seq)){
#'   epsi = epsi_seq[i]
#'   dist_thre = exp(epsi)
#'   epsi_nbhd = apply(dis_L2,1,function(x) sum(x<dist_thre))
#'   avg_epsi_nbhd = mean(epsi_nbhd)
#'   print(epsi)
#'   print(avg_epsi_nbhd)
#' }
#'
#' ## number of clusters
#' K = 8
#' k_M = 15
#' k_m = 4
#' pca_thre = 0.9
#'
#' set.seed(100)
#'
#' SwissRoll_GR = GeoRatio_graph(SwissRoll,k_max=k_M,k_min=k_m,trueDim,thresh=pca_thre,gmode="max")
#'
#' dis_L2 = dist(SwissRoll)
#' dis_G = as.dist(SwissRoll_GR$dis_G)
#' dis_GR = as.dist(SwissRoll_GR$dis_GR)
#'
#' hc_L2 = hclust(dis_L2,method="average")
#' hc_G = hclust(dis_G,method="average")
#' hc_GR = hclust(dis_GR,method="average")
#'
#' indi_KML2 = pam(SwissRoll,K)$clustering
#' indi_KMG = pam(dis_G,K,diss=T)$clustering
#' indi_KMGR = pam(dis_GR,K,diss=T)$clustering
#' indi_HL2 = cutree(hc_L2,k=K)
#' indi_HG = cutree(hc_G,k=K)
#' indi_HGR = cutree(hc_GR,k=K)
#'
#' lpca_V = kplanes(SwissRoll,K,trueDim,iter.max=100,thresh=1e-5)
#' indi_lpca = lpca_V$indicator
#'
#' open3d()
#' plot3d(SwissRoll,col=indi_HGR,main="HC-GR")
#' open3d()
#' plot3d(SwissRoll,col=indi_HG,main="HC-G")
#' open3d()
#' plot3d(SwissRoll,col=indi_HL2,main="HC-L2")
#' open3d()
#' plot3d(SwissRoll,col=indi_KMGR,main="KMC-GR")
#' open3d()
#' plot3d(SwissRoll,col=indi_KMG,main="KMC-G")
#' open3d()
#' plot3d(SwissRoll,col=indi_KML2,main="KMC-L2")
#' open3d()
#' plot3d(SwissRoll,col=indi_lpca,main="LPCA")
#'
#' ############################## example III: Open box
#' ## package for 3d plot
#' library(rgl)
#' ## package for pam() kmeans clustering
#' library(cluster)
#' library(igraph) ## used by GeoRatio_graph
#'
#' ## load data
#' data(OpenBox)
#'
#' ## intrinsic dimension of the data
#' trueDim = 2
#'
#' ############## correlation dimension estimation
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
#' par(mfrow=c(1,1))
#' plot(est$x[2:(length(est$x)-1)],est$deri,type="l",ylim=c(0,10),xlab="log(epsilon)",ylab="est. correlation dimension",
#'      main="Estimated correlation dimension")
#' abline(trueDim,0,lty=3,col="gray60",lwd=2)
#'
#' ##### decide k_M (plateau log(epsi) from -4.8 to -2.8)
#' epsi_s = -3.5
#' epsi_e = 0
#' epsi_seq = seq(epsi_s,epsi_e,length=5)
#'
#' dis_L2 = as.matrix(dist(OpenBox))
#'
#' for(i in 1:length(epsi_seq)){
#'   epsi = epsi_seq[i]
#'   dist_thre = exp(epsi)
#'   epsi_nbhd = apply(dis_L2,1,function(x) sum(x<dist_thre))
#'   avg_epsi_nbhd = mean(epsi_nbhd)
#'   print(epsi)
#'   print(avg_epsi_nbhd)
#' }
#'
#' ## number of clusters
#' K = 6
#' k_M = 20
#' k_m = 4
#' pca_thre = 0.9
#'
#' set.seed(100)
#'
#' OpenBox_GR = GeoRatio_graph(OpenBox,k_max=k_M,k_min=k_m,trueDim,thresh=pca_thre,gmode="max")
#'
#' dis_L2 = dist(OpenBox)
#' dis_G = as.dist(OpenBox_GR$dis_G)
#' dis_GR = as.dist(OpenBox_GR$dis_GR)
#'
#' hc_L2 = hclust(dis_L2,method="average")
#' hc_G = hclust(dis_G,method="average")
#' hc_GR = hclust(dis_GR,method="average")
#'
#' indi_KML2 = pam(dis_L2,K,diss=T)$clustering
#' indi_KMG = pam(dis_G,K,diss=T)$clustering
#' indi_KMGR = pam(dis_GR,K,diss=T)$clustering
#' indi_HL2 = cutree(hc_L2,k=K)
#' indi_HG = cutree(hc_G,k=K)
#' indi_HGR = cutree(hc_GR,k=K)
#'
#' lpca_V = kplanes(OpenBox,K,trueDim,iter.max=100,thresh=1e-5)
#' indi_lpca = lpca_V$indicator
#'
#' open3d()
#' plot3d(OpenBox,col=indi_HGR,main="HC-GR",xlim=c(0,2),ylim=c(0,2),zlim=c(0,2))
#' open3d()
#' plot3d(OpenBox,col=indi_HG,main="HC-G",xlim=c(0,2),ylim=c(0,2),zlim=c(0,2))
#' open3d()
#' plot3d(OpenBox,col=indi_HL2,main="HC-L2",xlim=c(0,2),ylim=c(0,2),zlim=c(0,2))
#' open3d()
#' plot3d(OpenBox,col=indi_KMGR,main="KMC-GR",xlim=c(0,2),ylim=c(0,2),zlim=c(0,2))
#' open3d()
#' plot3d(OpenBox,col=indi_KMG,main="KMC-G",xlim=c(0,2),ylim=c(0,2),zlim=c(0,2))
#' open3d()
#' plot3d(OpenBox,col=indi_KML2,main="KMC-L2",xlim=c(0,2),ylim=c(0,2),zlim=c(0,2))
#' open3d()
#' plot3d(OpenBox,col=indi_lpca,main="LPCA",xlim=c(0,2),ylim=c(0,2),zlim=c(0,2))
#'
#' @export
#'
GeoRatio_graph = function(X,k_max=k_M,k_min=k_m,d,thresh=0.9,gmode="max",distance=NaN){
  ## X: data matrix; k: initial number of nbhd points; d: dim of embedding space; epsi: regulary term added to avoid singularity.
  ## thresh: local pca thresh

  library("igraph")

  n=nrow(X)
  N=ncol(X)

  ## compute Euclidean distance if distance is not provided
  if(is.nan(distance)==TRUE){
    dis2 = as.matrix(dist(X))
  } else {
    dis2 = distance
  }
  #   tol = epsi
  #   X2 = apply(X^2,1,sum)
  #   dis2 = ((matrix(rep(X2,n),nrow=n,ncol=n,byrow=T)+matrix(rep(X2,n),nrow=n,ncol=n,byrow=F)-2*X%*%t(X))/500)
  #   dis2 = sqrt(dis2)
  #   dis2 = as.matrix(dist(X))
  #   diag(dis2)=0
  #   mm = apply(X,1,mean)  ## data mean
  #   ss = apply(X,1,var)  ## data variance
  #   ssm = matrix(rep(mm  ,N),n,N,byrow=F)

  #### get the minimum spanning tree of the data
  gall = graph.adjacency(dis2,weight=TRUE,mode=gmode)
  mst = minimum.spanning.tree(gall)

  nbhd = t(apply(dis2,2,order))
  size = rep(0,n)  ## size of nbhd of each i with i-th data included
  for(i in 1:n){
    if(k_min>d+1){ # used to be d
      current = k_min
    }else{
      current = d+2 # used to be d+1
    }
    X.temp = X[nbhd[i,1:current],]
    prop.c = 1e-6
    prop.old = 0

    ## adaptively grow nbhd size
    while(current<=k_max&(prop.c>thresh|prop.c>prop.old)){
      prop.old = prop.c
      temp.mean = apply(X.temp,2,mean)
      n.temp = nrow(X.temp)
      temp.svd = svd((X.temp-as.vector(rep(1,n.temp))%*%t(as.vector(temp.mean)))/sqrt(n.temp))
      temp.eigenval = temp.svd$d^2/N  ## eigenvalues
      prop.c = cumsum(temp.eigenval)[d]/sum(temp.eigenval)
      current = current+1
      X.temp = X[nbhd[i,1:current],]
    }
    size[i] = current-2
  }

  #   X.rep = matrix(0,n,N)
  #   M = diag(1,n)
  #   W = matrix(0,n,n)
  #   L2.error = rep(0,n)
  #   S2 = rep(0,n)
  #   R2 = rep(0,n)
  #   for(i in 1:n){
  #     k = size[i]-1
  #     nbhd.c = nbhd[i,2:(k+1)]
  #     z = X[nbhd.c,]-matrix(rep(X[i,]),nrow=k,ncol=N,byrow=T)
  #     C = z%*%t(z)
  #
  #     if(method=="trace"){
  #       C = C+diag(1,k)*tol*sum(diag(C))
  #       #print("trace")
  #     }else if(method=="diag"){
  #       C = C + tol*diag(diag(C))
  #       #print("diag")
  #     }else if(method=="trace_k"){
  #       C = C + diag(1,k)*tol*sum(diag(C))*k
  #       #print("diag_k")
  #     }
  #     w.c = solve(C,c(rep(1,k)))
  #     w.c = w.c/sum(w.c)
  #
  #     W[i,nbhd.c]=w.c
  #     #W.adp[i,1:k]=w.c
  #     if(length(w.c)==1){
  #       X.rep[i,] = w.c*X[nbhd.c,]
  #     }else{
  #       X.rep[i,] = w.c%*%X[nbhd.c,]
  #     }
  #     L2.error[i] = sqrt(mean((X[i,]-X.rep[i,])^2))
  #     R2[i] = 1-L2.error[i]^2/ss[i]
  #     S2[i] = mean((X.rep[i,]-ssm)^2)
  #     j = nbhd.c
  #     M[i,j] = M[i,j]-w.c
  #     M[j,i] = M[j,i]-w.c
  #     M[j,j] = M[j,j]+w.c%*%t(w.c)
  #   }
  #   eig.adp = eigen(M)
  #   Y = eig.adp$vectors[,(n-d):(n-1)]
  #   eigenvalue.adp = eig.adp$values

  #### Construct a graph using nbhd defined by nbhd_adp
  adjm = matrix(0,n,n)
  for(i in 1:n){
    adjm[i,nbhd[i,2:size[i]]] = dis2[i,nbhd[i,2:size[i]]]
  }

  mst_adjm = get.adjacency(mst)
  adjm_adj = adjm>0
  adjm = adjm+(mst_adjm>adjm_adj)*mst_adjm*dis2 ##### union minimum spanning tree
  graph_adp <- graph.adjacency(adjm, weighted=TRUE,mode=gmode)
  dist_adp = shortest.paths(graph_adp)
  dist_gr = log(dist_adp/as.matrix(dis2))
  diag(dist_gr) = 0

  output = list(dis_GR = dist_gr, dis_G=dist_adp, adjm=adjm, graph_adp=graph_adp, nbhd_size=size)

  class(output) = "GeoRatio_graph"
  output
}


