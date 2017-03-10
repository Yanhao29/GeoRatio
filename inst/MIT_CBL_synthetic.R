##################################
## read in data
##################################
library(pixmap)  ## read in png file
library(cluster)  ## pam clustering
library(clusteval)  ## silhouette information
library(igraph)  ## used by GeoRatio package
library(GeoRatio)
library(xtable)  ## generate latex table

path = getwd()
subject = 1 ## {1,2,...,10}

filenames <- list.files(paste0(path,'/inst/extdata/training-synthetic/'), pattern=paste0("000",toString(subject-1),".*\\.pgm"), full.names=TRUE)
ldf <- lapply(filenames, read.pnm)

## sample size
n = length(ldf)

## number of columns (pixels of each image)
size = dim(ldf[[1]]@grey)
N = size[1]*size[1]

## matrix of images (each row is a vectorized image)
X = matrix(0,n,N)
for(i in 1:n){
  X[i,] = as.vector(ldf[[i]]@grey)
}

#######################################
## global pca (nominal dimension is high)
#######################################
ptm = proc.time()  ## record SVD time for data
X.mean = apply(X,2,mean)
X.svd = svd((X-as.vector(rep(1,n))%*%t(as.vector(X.mean)))/sqrt(n))
X.eigenvec = X.svd$v*sqrt(N)      ## eigenfunctions
X.pcscore.norm = X.svd$u*sqrt(n)       ## normalized PC scores (unit variance in each coordinate)
X.eigenval = X.svd$d^2/N  ## eigenvalues
X.pcscore = X.pcscore.norm %*% diag(X.svd$d/sqrt(N))  ## PC scores
svd_time = proc.time() - ptm

X.eigenval.cumsum = cumsum(X.eigenval)/sum(X.eigenval)

## number of PCs needed to explain at least 0.95 of data variance
v.num95 = which(X.eigenval.cumsum>0.95)[1]

## number of PCs needed to explain at least 0.999 of data variance
v.num999 = which(X.eigenval.cumsum>0.999)[1]

plot(X.eigenval,type="b",pch=1,xlab="index",ylab="eigenvalue")
abline(v=v.num95,lty=3,col="gray60")

#######################################
## correlation dimension estimation
## (intinsic dimension is 3)
#######################################
data = X.pcscore[,1:v.num999]   ## use PCs up to v.num999 (could use all PCs)
num = 50
dis_L2 = dist(data, method = "euclidean", diag = FALSE, upper = TRUE)
est_dim = correDimEst(num,dis_L2)

plot(est_dim$x[2:(length(est_dim$x)-1)],est_dim$deri,type="l"
     ,xlab="logarithm of epsilon",ylab="intrinsic dimension")
abline(3,0,lty=3,col="gray60")

#######################################
## GeoRatio distance estimation
#######################################

##### decide k_M using result correlation dimension estimation
## epsilon range of the plateau at hight 3
epsi_s = -4
epsi_e = -2.5
epsi_seq = seq(epsi_s,epsi_e,length=5)

dis_L2 = as.matrix(dis_L2)
for(i in 1:length(epsi_seq)){
  epsi = epsi_seq[i]
  dist_thre = exp(epsi)
  epsi_nbhd = apply(dis_L2,1,function(x) sum(x<dist_thre))
  avg_epsi_nbhd = mean(epsi_nbhd)
  print(epsi)
  print(avg_epsi_nbhd)
}

##### GeoRatio graph building
k_m = 5
k_M = 10 ## maximum # of nbhd size
pca_thre = 0.9
trueDim = 3

test = GeoRatio_graph(data,k_max=k_M,k_min=k_m,trueDim,thresh=pca_thre,gmode="max")  ### k around 10, epsi about(>) 1e-3 good results

is.connected(test$graph_adp)

## construct GeoRatio distance
dis_G = test$dis_G
dis_GR = log(dis_G/dis_L2)
diag(dis_GR) = 0

#######################################
## clustering using GeoRatio
## look at reconstruction error,
## silhouette information,
## dimension reduction efficiency,
## clustering performance (face poses i.e. K_c=9)
#######################################
## increase number of clusters from 1 to K
K = 15

## 9 face poses
K_c = 9

##### k means on Euclidean distance
sil_KMCL2 = rep(0,K)
error_pc_KMCL2 = rep(0,K)
for(i in 1:K){
  print(i)
  pamm = pam(data,i)
  if(i==1) {sil_KMCL2[i]=0}
  else
  {sil_KMCL2[i] = pamm$silinfo$avg.width}
  indi = pamm$clustering
  now = lpca(indi,data,d=trueDim)
  error_pc_KMCL2[i]=now$mean_error
}

##### hierarchical on Euclidean distance
h_L2 = hclust(as.dist(dis_L2),method="average")

sil_HCL2 = rep(0,K)
error_pc_HCL2 = rep(0,K)
for(i in 1:K){
  print(i)
  indi = cutree(h_L2,i)
  if(i==1) {sil_HCL2[i]=0}
  else
  {sil_HCL2[i] = summary(silhouette(indi,dist=dis_L2))$avg.width}

  now = lpca(indi,data,d=trueDim)
  error_pc_HCL2[i]=now$mean_error
}

##### k means on geodesic distance
sil_KMCG = rep(0,K)
error_pc_KMCG = rep(0,K)
for(i in 1:K){
  print(i)
  pamm = pam(dis_G,i,diss = TRUE)
  if(i==1) {sil_KMCG[i]=0}
  else
  {sil_KMCG[i] = pamm$silinfo$avg.width}
  indi = pamm$clustering
  now = lpca(indi,data,d=trueDim)
  error_pc_KMCG[i]=now$mean_error
}

##### hierarchical on geodesic distance
h_G = hclust(as.dist(dis_G),method="average")

sil_HCG = rep(0,K)
error_pc_HCG = rep(0,K)
for(i in 1:K){
  print(i)
  indi = cutree(h_G,i)
  if(i==1) {sil_HCG[i]=0}
  else
  {sil_HCG[i] = summary(silhouette(indi,dist=dis_G))$avg.width}

  now = lpca(indi,data,d=trueDim)
  error_pc_HCG[i]=now$mean_error
}

##### k means on GeoRatio distance
sil_KMCGR = rep(0,K)
error_pc_KMCGR = rep(0,K)
for(i in 1:K){
  print(i)
  pamm = pam(dis_GR,i,diss = TRUE)
  if(i==1) {sil_KMCGR[i]=0}
  else
  {sil_KMCGR[i] = pamm$silinfo$avg.width}
  indi = pamm$clustering
  now = lpca(indi,data,d=trueDim)
  error_pc_KMCGR[i]=now$mean_error
}

##### hierarchical on GeoRatio distance
h_GR = hclust(as.dist(dis_GR),method="average")

sil_HCGR = rep(0,K)
error_pc_HCGR = rep(0,K)
for(i in 1:K){
  print(i)
  indi = cutree(h_GR,i)
  if(i==1) {sil_HCGR[i]=0}
  else
  {sil_HCGR[i] = summary(silhouette(indi,dist=dis_GR))$avg.width}

  now = lpca(indi,data,d=trueDim)
  error_pc_HCGR[i]=now$mean_error
}

##### local pca
sil_LPCA = rep(0,K)
error_pc_LPCA = rep(0,K)
for(i in 1:K){
  print(i)
  temp = kplanes(data,i,trueDim,iter.max=100,thresh=1e-4)
  indi = temp$indicator
  if(i==1) {sil_LPCA[i]=0}
  else
  {sil_LPCA[i] = summary(silhouette(indi,as.dist(dis_L2)))$avg.width}
  error_pc_LPCA[i] = temp$error
  print("%%%%%")
}

par(oma=c(0,0,0,0),mar=c(5,5,1,1),mfrow=c(1,2))
# limi = max(rbind(sil.h.GR,sil.h.G,sil.h.L2,sil.km.GR,sil.km.G,sil.km.L2))
limi_M = max(rbind(sil_HCGR,sil_HCG,sil_HCL2))
limi_m = min(rbind(sil_HCGR,sil_HCG,sil_HCL2))
plot(1:K,sil_HCGR,type="l",pch=1,lty=1,col=2,lwd=2,ylim=c(limi_m,limi_M),xlab="number of clusters",ylab="sil. information")
points(1:K,sil_HCG,type="l",pch=1,lty=1,col=3)
points(1:K,sil_HCL2,type="l",pch=1,lty=1,col=1)
# points(1:K,sil.lpca,type="l",pch=1,lty=2,col=4)
# points(1:K,sil.km.GR,type="b",pch=2,lty=1,col=2)
# points(1:K,sil.km.G,type="b",pch=2,lty=1,col=3)
# points(1:K,sil.km.L2,type="b",pch=2,lty=1,col=1)
abline(v=K_c,col='gray60',lty=3, lwd = 2)
# legend(6,3,c("HC-GeoRatio","HC-geodesic","HC-L2","KMC-GeoRatio","KMC-geodesic","KMC-L2"),pch=c(1,1,1,2,2,2),col=c(2,3,1,2,3,1),lwd=c(2,1,1,1,1,1))
# legend(6,3,c("HC-GeoRatio","HC-geodesic","HC-L2","LPCA"),col=c(2,3,1,4),lwd=c(2,1,1,1))

limi_M = max(rbind(error_pc_HCGR,error_pc_HCG,error_pc_HCL2,error_pc_LPCA))
limi_m = min(rbind(error_pc_HCGR,error_pc_HCG,error_pc_HCL2,error_pc_LPCA))
plot(1:K,error_pc_HCGR,type="l",pch=1,lty=1,col=2,lwd=2,ylim=c(limi_m,limi_M),xlab="number of clusters",ylab="error")
points(1:K,error_pc_HCG,type="l",pch=1,lty=1,col=3)
points(1:K,error_pc_HCL2,type="l",pch=1,lty=1,col=1)
# points(1:K,rep.e.km.GR,type="b",pch=2,lty=1,col=2)
# points(1:K,rep.e.km.G,type="b",pch=2,lty=1,col=3)
# points(1:K,rep.e.km.L2,type="b",pch=2,lty=1,col=1)
points(1:K,error_pc_LPCA,type="l",pch=1,lty=2,col=4)
abline(v=K_c,col='gray60',lty=3, lwd = 2)
# legend(6,1,c("HC-GeoRatio","HC-geodesic","HC-L2","KMC-GeoRatio","KMC-geodesic","KMC-L2"),pch=c(1,1,1,2,2,2),col=c(2,3,1,2,3,1),lwd=c(2,1,1,1,1,1))
legend(10,limi_M,c("HC-GR","HC-G","HC-L2","LPCA"),col=c(2,3,1,4),lwd=c(2,1,1,1),lty=c(1,1,1,2))

########################################
#####  clustering performance at K_c = 9
########################################

##### k-means clustering on eucleadian distance
indi_KMCL2 = pam(data,K_c)$clustering
now_KMCL2 = lpca(indi_KMCL2,data,d=trueDim)

##### k-means clustering on geodesic distance
pamm = pam(dis_G,K_c,diss=T)
indi_KMCG = pamm$clustering
now_KMCG = lpca(indi_KMCG,data,d=trueDim)

##### k-means clustering on GeoRatio distance
pamm = pam(dis_GR,K_c,diss=T)
indi_KMCGR = pamm$clustering
now_KMCGR = lpca(indi_KMCGR,data,d=trueDim)

##### hierachical clustering on eucleadian distance
indi_HCL2 = cutree(h_L2,k=K_c)
now_HCL2 = lpca(indi_HCL2,data,d=trueDim)

##### hierachical clustering on geodesic distance
indi_HCG = cutree(h_G,k=K_c)
now_HCG = lpca(indi_HCG,data,d=trueDim)

##### hirarchical clustering on GeoRatio distance
indi_HCGR = cutree(h_GR,k=K_c)
now_HCGR = lpca(indi_HCGR,data,d=trueDim)

##### local pca
LPCA_t = kplanes(data,K_c,trueDim,iter.max=100,thresh=1e-4)
indi_LPCA = LPCA_t$indicator

## true cluster membership
col.f = rep(1:9,each=36)

HCGR_cluster = matrix(0,K_c,K_c) # row cluster membership, column true membership
HCG_cluster = matrix(0,K_c,K_c)
HCL2_cluster = matrix(0,K_c,K_c)
KMCGR_cluster = matrix(0,K_c,K_c)
KMCG_cluster = matrix(0,K_c,K_c)
KMCL2_cluster = matrix(0,K_c,K_c)
LPCA_cluster = matrix(0,K_c,K_c)

for(i in 1:K_c){
  for(j in 1:K_c){
    HCGR_cluster[i,j] = length(intersect(which(indi_HCGR==i),(((j-1)*36+1):(j*36))))
    HCG_cluster[i,j] = length(intersect(which(indi_HCG==i),(((j-1)*36+1):(j*36))))
    HCL2_cluster[i,j] = length(intersect(which(indi_HCL2==i),(((j-1)*36+1):(j*36))))
    KMCGR_cluster[i,j] = length(intersect(which(indi_KMCGR==i),(((j-1)*36+1):(j*36))))
    KMCG_cluster[i,j] = length(intersect(which(indi_KMCG==i),(((j-1)*36+1):(j*36))))
    KMCL2_cluster[i,j] = length(intersect(which(indi_KMCL2==i),(((j-1)*36+1):(j*36))))
    LPCA_cluster[i,j] = length(intersect(which(indi_LPCA==i),(((j-1)*36+1):(j*36))))
  }
}

HCGR_assign = apply(HCGR_cluster, 1, function(x) which.max(x))
HCG_assign = apply(HCG_cluster, 1, function(x) which.max(x))
HCL2_assign = apply(HCL2_cluster, 1, function(x) which.max(x))
KMCGR_assign = apply(KMCGR_cluster, 1, function(x) which.max(x))
KMCG_assign = apply(KMCG_cluster, 1, function(x) which.max(x))
KMCL2_assign = apply(KMCL2_cluster, 1, function(x) which.max(x))
LPCA_assign = apply(LPCA_cluster, 1, function(x) which.max(x))

indi_HCGR_a = rep(0,n)
indi_HCG_a = rep(0,n)
indi_HCL2_a = rep(0,n)
indi_KMCGR_a = rep(0,n)
indi_KMCG_a = rep(0,n)
indi_KMCL2_a = rep(0,n)
indi_LPCA_a = rep(0,n)

for(i in 1:K_c){
  indi_HCGR_a[indi_HCGR==i] = HCGR_assign[i]
  indi_HCG_a[indi_HCG==i] = HCG_assign[i]
  indi_HCL2_a[indi_HCL2==i] = HCL2_assign[i]
  indi_KMCGR_a[indi_KMCGR==i] = KMCGR_assign[i]
  indi_KMCG_a[indi_KMCG==i] = KMCG_assign[i]
  indi_KMCL2_a[indi_KMCL2==i] = KMCL2_assign[i]
  indi_LPCA_a[indi_LPCA==i] = LPCA_assign[i]
}

rd_HCGR = cluster_similarity(indi_HCGR_a,col.f,similarity="rand")
rd_HCG = cluster_similarity(indi_HCG_a,col.f,similarity="rand")
rd_HCL2 = cluster_similarity(indi_HCL2_a,col.f,similarity="rand")
rd_KMCGR = cluster_similarity(indi_KMCGR_a,col.f,similarity="rand")
rd_KMCG = cluster_similarity(indi_KMCG_a,col.f,similarity="rand")
rd_KMCL2 = cluster_similarity(indi_KMCL2_a,col.f,similarity="rand")
rd_LPCA = cluster_similarity(indi_LPCA_a,col.f,similarity="rand")

HCGR_mis = mean((indi_HCGR_a-col.f)!=0)
HCG_mis = mean((indi_HCG_a-col.f)!=0)
HCL2_mis = mean((indi_HCL2_a-col.f)!=0)
KMCGR_mis = mean((indi_KMCGR_a-col.f)!=0)
KMCG_mis = mean((indi_KMCG_a-col.f)!=0)
KMCL2_mis = mean((indi_KMCL2_a-col.f)!=0)
LPCA_mis = mean((indi_LPCA_a-col.f)!=0)

dist_true = dist(col.f, method='euclidean')>0
dist_HCGR = dist(indi_HCGR, method='euclidean')>0
dist_HCG = dist(indi_HCG, method='euclidean')>0
dist_HCL2 = dist(indi_HCL2, method='euclidean')>0
dist_KMCGR = dist(indi_KMCGR, method='euclidean')>0
dist_KMCG = dist(indi_KMCG, method='euclidean')>0
dist_KMCL2 = dist(indi_KMCL2, method='euclidean')>0
dist_LPCA = dist(indi_LPCA, method='euclidean')>0

HCGR_tp = sum((dist_HCGR==dist_true)&(dist_true==0))/sum(dist_true==0)
HCGR_fp = sum(dist_HCGR<dist_true)/sum(dist_HCGR==0)
HCGR_tn = sum((dist_HCGR==dist_true)&(dist_true==1))/sum(dist_true==1)
HCGR_fn = sum(dist_HCGR>dist_true)/sum(dist_HCGR==1)

HCG_tp = sum((dist_HCG==dist_true)&(dist_true==0))/sum(dist_true==0)
HCG_fp = sum(dist_HCG<dist_true)/sum(dist_HCG==0)
HCG_tn = sum((dist_HCG==dist_true)&(dist_true==1))/sum(dist_true==1)
HCG_fn = sum(dist_HCG>dist_true)/sum(dist_HCG==1)

HCL2_tp = sum((dist_HCL2==dist_true)&(dist_true==0))/sum(dist_true==0)
HCL2_fp = sum(dist_HCL2<dist_true)/sum(dist_HCL2==0)
HCL2_tn = sum((dist_HCL2==dist_true)&(dist_true==1))/sum(dist_true==1)
HCL2_fn = sum(dist_HCL2>dist_true)/sum(dist_HCL2==1)

KMCGR_tp = sum((dist_KMCGR==dist_true)&(dist_true==0))/sum(dist_true==0)
KMCGR_fp = sum(dist_KMCGR<dist_true)/sum(dist_KMCGR==0)
KMCGR_tn = sum((dist_KMCGR==dist_true)&(dist_true==1))/sum(dist_true==1)
KMCGR_fn = sum(dist_KMCGR>dist_true)/sum(dist_KMCGR==1)

KMCG_tp = sum((dist_KMCG==dist_true)&(dist_true==0))/sum(dist_true==0)
KMCG_fp = sum(dist_KMCG<dist_true)/sum(dist_KMCG==0)
KMCG_tn = sum((dist_KMCG==dist_true)&(dist_true==1))/sum(dist_true==1)
KMCG_fn = sum(dist_KMCG>dist_true)/sum(dist_KMCG==1)

KMCL2_tp = sum((dist_KMCL2==dist_true)&(dist_true==0))/sum(dist_true==0)
KMCL2_fp = sum(dist_KMCL2<dist_true)/sum(dist_KMCL2==0)
KMCL2_tn = sum((dist_KMCL2==dist_true)&(dist_true==1))/sum(dist_true==1)
KMCL2_fn = sum(dist_HCL2>dist_true)/sum(dist_KMCL2==1)

LPCA_tp = sum((dist_LPCA==dist_true)&(dist_true==0))/sum(dist_true==0)
LPCA_fp = sum(dist_LPCA<dist_true)/sum(dist_LPCA==0)
LPCA_tn = sum((dist_LPCA==dist_true)&(dist_true==1))/sum(dist_true==1)
LPCA_fn = sum(dist_LPCA>dist_true)/sum(dist_LPCA==1)

#######################################
## cluster eigen numbers
#######################################

##### HCGR
thresh_t = 0.95
ev_HCGR = NULL
evcumu_HCGR = NULL
ef_HCGR = NULL
pc_HCGR = NULL
mean_HCGR = NULL
eigennum_HCGR = rep(0,K_c)

for(i in 1:K_c){
  X_t = X[indi_HCGR==i,]
  n_t = nrow(X_t)
  mean_HCGR[[i]] = apply(X_t,2,mean)
  X_t.svd = svd((X_t-as.vector(rep(1,n_t))%*%t(as.vector(mean_HCGR[[i]])))/sqrt(n_t))

  ef_HCGR[[i]] = X_t.svd$v*sqrt(N)      ## eigenfunctions
  X_t.pcscore.norm = X_t.svd$u*sqrt(n_t)       ## normalized PC scores (unit variance in each coordinate)
  ev_HCGR[[i]] = X_t.svd$d^2/N  ## eigenvalues
  pc_HCGR[[i]] = X_t.pcscore.norm %*% diag(X_t.svd$d/sqrt(N))  ## PC scores

  evcumu_HCGR[[i]] = cumsum(ev_HCGR[[i]])/sum(ev_HCGR[[i]])
  eigennum_HCGR[i] = which(evcumu_HCGR[[i]]>thresh_t)[1]
}

##### HCG
thresh_t = 0.95
ev_HCG = NULL
evcumu_HCG = NULL
ef_HCG = NULL
pc_HCG = NULL
mean_HCG = NULL
eigennum_HCG = rep(0,K_c)

for(i in 1:K_c){
  X_t = X[indi_HCG==i,]
  n_t = nrow(X_t)
  mean_HCG[[i]] = apply(X_t,2,mean)
  X_t.svd = svd((X_t-as.vector(rep(1,n_t))%*%t(as.vector(mean_HCG[[i]])))/sqrt(n_t))

  ef_HCG[[i]] = X_t.svd$v*sqrt(N)      ## eigenfunctions
  X_t.pcscore.norm = X_t.svd$u*sqrt(n_t)       ## normalized PC scores (unit variance in each coordinate)
  ev_HCG[[i]] = X_t.svd$d^2/N  ## eigenvalues
  pc_HCG[[i]] = X_t.pcscore.norm %*% diag(X_t.svd$d/sqrt(N))  ## PC scores

  evcumu_HCG[[i]] = cumsum(ev_HCG[[i]])/sum(ev_HCG[[i]])
  eigennum_HCG[i] = which(evcumu_HCG[[i]]>thresh_t)[1]
}

##### HCL2
thresh_t = 0.95
ev_HCL2 = NULL
evcumu_HCL2 = NULL
ef_HCL2 = NULL
pc_HCL2 = NULL
mean_HCL2 = NULL
eigennum_HCL2 = rep(0,K_c)

for(i in 1:K_c){
  X_t = X[indi_HCL2==i,]
  n_t = nrow(X_t)
  mean_HCL2[[i]] = apply(X_t,2,mean)
  X_t.svd = svd((X_t-as.vector(rep(1,n_t))%*%t(as.vector(mean_HCL2[[i]])))/sqrt(n_t))

  ef_HCL2[[i]] = X_t.svd$v*sqrt(N)      ## eigenfunctions
  X_t.pcscore.norm = X_t.svd$u*sqrt(n_t)       ## normalized PC scores (unit variance in each coordinate)
  ev_HCL2[[i]] = X_t.svd$d^2/N  ## eigenvalues
  pc_HCL2[[i]] = X_t.pcscore.norm %*% diag(X_t.svd$d/sqrt(N))  ## PC scores

  evcumu_HCL2[[i]] = cumsum(ev_HCL2[[i]])/sum(ev_HCL2[[i]])
  eigennum_HCL2[i] = which(evcumu_HCL2[[i]]>thresh_t)[1]
}

##### KMCGR
thresh_t = 0.95
ev_KMCGR = NULL
evcumu_KMCGR = NULL
ef_KMCGR = NULL
pc_KMCGR = NULL
mean_KMCGR = NULL
eigennum_KMCGR = rep(0,K_c)

for(i in 1:K_c){
  X_t = X[indi_KMCGR==i,]
  n_t = nrow(X_t)
  mean_KMCGR[[i]] = apply(X_t,2,mean)
  X_t.svd = svd((X_t-as.vector(rep(1,n_t))%*%t(as.vector(mean_KMCGR[[i]])))/sqrt(n_t))

  ef_KMCGR[[i]] = X_t.svd$v*sqrt(N)      ## eigenfunctions
  X_t.pcscore.norm = X_t.svd$u*sqrt(n_t)       ## normalized PC scores (unit variance in each coordinate)
  ev_KMCGR[[i]] = X_t.svd$d^2/N  ## eigenvalues
  pc_KMCGR[[i]] = X_t.pcscore.norm %*% diag(X_t.svd$d/sqrt(N))  ## PC scores

  evcumu_KMCGR[[i]] = cumsum(ev_KMCGR[[i]])/sum(ev_KMCGR[[i]])
  eigennum_KMCGR[i] = which(evcumu_KMCGR[[i]]>thresh_t)[1]
}

##### KMCG
thresh_t = 0.95
ev_KMCG = NULL
evcumu_KMCG = NULL
ef_KMCG = NULL
pc_KMCG = NULL
mean_KMCG = NULL
eigennum_KMCG = rep(0,K_c)

for(i in 1:K_c){
  X_t = X[indi_KMCG==i,]
  n_t = nrow(X_t)
  mean_KMCG[[i]] = apply(X_t,2,mean)
  X_t.svd = svd((X_t-as.vector(rep(1,n_t))%*%t(as.vector(mean_KMCG[[i]])))/sqrt(n_t))

  ef_KMCG[[i]] = X_t.svd$v*sqrt(N)      ## eigenfunctions
  X_t.pcscore.norm = X_t.svd$u*sqrt(n_t)       ## normalized PC scores (unit variance in each coordinate)
  ev_KMCG[[i]] = X_t.svd$d^2/N  ## eigenvalues
  pc_KMCG[[i]] = X_t.pcscore.norm %*% diag(X_t.svd$d/sqrt(N))  ## PC scores

  evcumu_KMCG[[i]] = cumsum(ev_KMCG[[i]])/sum(ev_KMCG[[i]])
  eigennum_KMCG[i] = which(evcumu_KMCG[[i]]>thresh_t)[1]
}

##### KMCL2
thresh_t = 0.95
ev_KMCL2 = NULL
evcumu_KMCL2 = NULL
ef_KMCL2 = NULL
pc_KMCL2 = NULL
mean_KMCL2 = NULL
eigennum_KMCL2 = rep(0,K_c)

for(i in 1:K_c){
  X_t = X[indi_KMCL2==i,]
  n_t = nrow(X_t)
  mean_KMCL2[[i]] = apply(X_t,2,mean)
  X_t.svd = svd((X_t-as.vector(rep(1,n_t))%*%t(as.vector(mean_KMCL2[[i]])))/sqrt(n_t))

  ef_KMCL2[[i]] = X_t.svd$v*sqrt(N)      ## eigenfunctions
  X_t.pcscore.norm = X_t.svd$u*sqrt(n_t)       ## normalized PC scores (unit variance in each coordinate)
  ev_KMCL2[[i]] = X_t.svd$d^2/N  ## eigenvalues
  pc_KMCL2[[i]] = X_t.pcscore.norm %*% diag(X_t.svd$d/sqrt(N))  ## PC scores

  evcumu_KMCL2[[i]] = cumsum(ev_KMCL2[[i]])/sum(ev_KMCL2[[i]])
  eigennum_KMCL2[i] = which(evcumu_KMCL2[[i]]>thresh_t)[1]
}

##### LPCA
thresh_t = 0.95
ev_LPCA = NULL
evcumu_LPCA = NULL
ef_LPCA = NULL
pc_LPCA = NULL
mean_LPCA = NULL
eigennum_LPCA = rep(0,K_c)

for(i in 1:K_c){
  X_t = X[indi_LPCA==i,]
  n_t = nrow(X_t)
  mean_LPCA[[i]] = apply(X_t,2,mean)
  X_t.svd = svd((X_t-as.vector(rep(1,n_t))%*%t(as.vector(mean_LPCA[[i]])))/sqrt(n_t))

  ef_LPCA[[i]] = X_t.svd$v*sqrt(N)      ## eigenfunctions
  X_t.pcscore.norm = X_t.svd$u*sqrt(n_t)       ## normalized PC scores (unit variance in each coordinate)
  ev_LPCA[[i]] = X_t.svd$d^2/N  ## eigenvalues
  pc_LPCA[[i]] = X_t.pcscore.norm %*% diag(X_t.svd$d/sqrt(N))  ## PC scores

  evcumu_LPCA[[i]] = cumsum(ev_LPCA[[i]])/sum(ev_LPCA[[i]])
  eigennum_LPCA[i] = which(evcumu_LPCA[[i]]>thresh_t)[1]
}

e_rep_HCGR = rep(0,n)
e_rep_HCG = rep(0,n)
e_rep_HCL2 = rep(0,n)
e_rep_KMCGR = rep(0,n)
e_rep_KMCG = rep(0,n)
e_rep_KMCL2 = rep(0,n)
e_rep_LPCA = rep(0,n)
e_rep_gpca = rep(0,n)

rep_HCGR = matrix(0,n,N)
rep_HCG = matrix(0,n,N)
rep_HCL2 = matrix(0,n,N)
rep_KMCGR = matrix(0,n,N)
rep_KMCG = matrix(0,n,N)
rep_KMCL2 = matrix(0,n,N)
rep_LPCA = matrix(0,n,N)
rep_gpca = matrix(0,n,N)

num_ev = trueDim

for(idx in 1:n){

  ex = as.vector(X.pcscore[idx,1:v.num999]%*%(t(X.eigenvec)[1:v.num999,]) + X.mean)

  rep_gpca[idx,] = as.vector(X.pcscore[idx,1:num_ev]%*%(t(X.eigenvec)[1:num_ev,]) + X.mean)
  e_rep_gpca[idx] = mean((ex-rep_gpca[idx,])^2)/var(ex)

  cluster = indi_HCGR[idx]
  idx_t = which(c(1:n)[indi_HCGR==cluster]==idx)
  rep_HCGR[idx,] = as.vector(pc_HCGR[[cluster]][idx_t,1:num_ev]%*%(t(ef_HCGR[[cluster]])[1:num_ev,]) + mean_HCGR[[cluster]])
  e_rep_HCGR[idx] = mean((ex-rep_HCGR[idx,])^2)/var(ex)

  cluster = indi_HCG[idx]
  idx_t = which(c(1:n)[indi_HCG==cluster]==idx)
  rep_HCG[idx,] = as.vector(pc_HCG[[cluster]][idx_t,1:num_ev]%*%(t(ef_HCG[[cluster]])[1:num_ev,]) + mean_HCG[[cluster]])
  e_rep_HCG[idx] = mean((ex-rep_HCG[idx,])^2)/var(ex)

  cluster = indi_HCL2[idx]
  idx_t = which(c(1:n)[indi_HCL2==cluster]==idx)
  rep_HCL2[idx,] = as.vector(pc_HCL2[[cluster]][idx_t,1:num_ev]%*%(t(ef_HCL2[[cluster]])[1:num_ev,]) + mean_HCL2[[cluster]])
  e_rep_HCL2[idx] = mean((ex-rep_HCL2[idx,])^2)/var(ex)

  cluster = indi_KMCGR[idx]
  idx_t = which(c(1:n)[indi_KMCGR==cluster]==idx)
  rep_KMCGR[idx,] = as.vector(pc_KMCGR[[cluster]][idx_t,1:num_ev]%*%(t(ef_KMCGR[[cluster]])[1:num_ev,]) + mean_KMCGR[[cluster]])
  e_rep_KMCGR[idx] = mean((ex-rep_KMCGR[idx,])^2)/var(ex)

  cluster = indi_KMCG[idx]
  idx_t = which(c(1:n)[indi_KMCG==cluster]==idx)
  rep_KMCG[idx,] = as.vector(pc_KMCG[[cluster]][idx_t,1:num_ev]%*%(t(ef_KMCG[[cluster]])[1:num_ev,]) + mean_KMCG[[cluster]])
  e_rep_KMCG[idx] = mean((ex-rep_KMCG[idx,])^2)/var(ex)

  cluster = indi_KMCL2[idx]
  idx_t = which(c(1:n)[indi_KMCL2==cluster]==idx)
  rep_KMCL2[idx,] = as.vector(pc_KMCL2[[cluster]][idx_t,1:num_ev]%*%(t(ef_KMCL2[[cluster]])[1:num_ev,]) + mean_KMCL2[[cluster]])
  e_rep_KMCL2[idx] = mean((ex-rep_KMCL2[idx,])^2)/var(ex)

  cluster = indi_LPCA[idx]
  idx_t = which(c(1:n)[indi_LPCA==cluster]==idx)
  rep_LPCA[idx,] = as.vector(pc_LPCA[[cluster]][idx_t,1:num_ev]%*%(t(ef_LPCA[[cluster]])[1:num_ev,]) + mean_LPCA[[cluster]])
  e_rep_LPCA[idx] = mean((ex-rep_LPCA[idx,])^2)/var(ex)

  if(idx%%10==0){
    print(idx)
  }
}
mean_e_HCGR = mean(e_rep_HCGR)
mean_e_HCG = mean(e_rep_HCG)
mean_e_HCL2 = mean(e_rep_HCL2)
mean_e_KMCGR = mean(e_rep_KMCGR)
mean_e_KMCG = mean(e_rep_KMCG)
mean_e_KMCL2 = mean(e_rep_KMCL2)
mean_e_LPCA = mean(e_rep_LPCA)
mean_e_gpca = mean(e_rep_gpca)

mean_e_HCGR
mean_e_HCG
mean_e_HCL2
mean_e_KMCGR
mean_e_KMCG
mean_e_KMCL2
mean_e_LPCA
mean_e_gpca


#################################
## All tables
#################################
xtable(cbind(rd_HCGR,rd_HCG,rd_HCL2,rd_KMCGR,rd_KMCG,rd_KMCL2,rd_LPCA))
xtable(rbind(error_pc_HCGR,error_pc_KMCGR,error_pc_HCG,error_pc_KMCG,error_pc_HCL2,error_pc_KMCL2,error_pc_LPCA)*1000)
xtable(rbind(sil_HCGR,sil_KMCGR,sil_HCG,sil_KMCG,sil_HCL2,sil_KMCL2))
xtable(cbind(rbind(HCGR_tp,HCGR_fp,HCGR_tn,HCGR_fn),rbind(HCG_tp,HCG_fp,HCG_tn,HCG_fn),
             rbind(HCL2_tp,HCL2_fp,HCL2_tn,HCL2_fn),rbind(KMCGR_tp,KMCGR_fp,KMCGR_tn,KMCGR_fn),
             rbind(KMCG_tp,KMCG_fp,KMCG_tn,KMCG_fn),rbind(KMCL2_tp,KMCL2_fp,KMCL2_tn,KMCL2_fn),
             rbind(LPCA_tp,LPCA_fp,LPCA_tn,LPCA_fn)))
xtable(cbind(HCGR_mis,HCG_mis,HCL2_mis,KMCGR_mis,KMCG_mis,KMCL2_mis,LPCA_mis))

xtable(rbind(cbind(mean(eigennum_HCGR),mean(eigennum_HCG),mean(eigennum_HCL2),mean(eigennum_KMCGR),mean(eigennum_KMCG),mean(eigennum_KMCL2),mean(eigennum_LPCA))
             ,cbind(sd(eigennum_HCGR),sd(eigennum_HCG),sd(eigennum_HCL2),sd(eigennum_KMCGR),sd(eigennum_KMCG),sd(eigennum_KMCL2),sd(eigennum_LPCA))))



#################################
## All plots
#################################

## sample image
idx_temp = seq(from=1,to=n,n/K_c)
par(mfrow=c(2,5),oma=c(1,1,1,1),mar=c(1,1,1,1))
for(i in 1:K_c){
  image(t(apply(matrix(X[idx_temp[i],],nrow=200),2,rev)),col = grey(seq(0, 1, length = 256)),xaxt= "n", yaxt= "n")
}

## eigenfaces from global pca
par(mfrow=c(2,5),oma=c(1,1,1,1),mar=c(1,1,4,1))
for(j in 1:10){
  image(t(apply(matrix(X.eigenvec[,j],nrow=200),2,rev)),main=paste0(toString(j),"-PC: ",round(X.eigenval.cumsum[j],2)),col = grey(seq(0, 1, length = 256)),xaxt="n", yaxt="n",cex.main=3)
}

## eigenfaces from HCGR in each cluster
par(mfrow=c(3,K_c),mar=c(0,0,4,0),oma=c(1,4,4,1))
for(i in 1:3){
  for(j in 1:K_c){
    image(t(apply(matrix(ef_HCGR[[j]][,i],nrow=200),2,rev)),main=round(evcumu_HCGR[[j]][i],2),col = grey(seq(0, 1, length = 256)),xaxt="n", yaxt="n",cex.main=3)
  }
}
mtext("3rd PC",  side=2, at=0.12,line=1,cex=2.5, col="black",outer=TRUE)
mtext("2ed PC",  side=2, at=0.46,line=1,cex=2.5, col="black",outer=TRUE)
mtext("1st PC",  side=2, at=0.80,line=1,cex=2.5, col="black",outer=TRUE)
mtext("cluster 1",  side=3, at=0.06,line=1,cex=2.5, col="black",outer=TRUE)
mtext("cluster 2",  side=3, at=0.17,line=1,cex=2.5, col="black",outer=TRUE)
mtext("cluster 3",  side=3, at=0.28,line=1,cex=2.5, col="black",outer=TRUE)
mtext("cluster 4",  side=3, at=0.39,line=1,cex=2.5, col="black",outer=TRUE)
mtext("cluster 5",  side=3, at=0.50,line=1,cex=2.5, col="black",outer=TRUE)
mtext("cluster 6",  side=3, at=0.61,line=1,cex=2.5, col="black",outer=TRUE)
mtext("cluster 7",  side=3, at=0.72,line=1,cex=2.5, col="black",outer=TRUE)
mtext("cluster 8",  side=3, at=0.83,line=1,cex=2.5, col="black",outer=TRUE)
mtext("cluster 9",  side=3, at=0.94,line=1,cex=2.5, col="black",outer=TRUE)


## eigenfaces from HCG in each cluster
par(mfrow=c(3,K_c),mar=c(0,0,4,0),oma=c(1,4,4,1))
for(i in 1:3){
  for(j in 1:K_c){
    image(t(apply(matrix(ef_HCG[[j]][,i],nrow=200),2,rev)),main=round(evcumu_HCG[[j]][i],2),col = grey(seq(0, 1, length = 256)),xaxt="n", yaxt="n",cex.main=3)
  }
}
mtext("3rd PC",  side=2, at=0.12,line=1,cex=2.5, col="black",outer=TRUE)
mtext("2ed PC",  side=2, at=0.46,line=1,cex=2.5, col="black",outer=TRUE)
mtext("1st PC",  side=2, at=0.80,line=1,cex=2.5, col="black",outer=TRUE)
mtext("cluster 1",  side=3, at=0.06,line=1,cex=2.5, col="black",outer=TRUE)
mtext("cluster 2",  side=3, at=0.17,line=1,cex=2.5, col="black",outer=TRUE)
mtext("cluster 3",  side=3, at=0.28,line=1,cex=2.5, col="black",outer=TRUE)
mtext("cluster 4",  side=3, at=0.39,line=1,cex=2.5, col="black",outer=TRUE)
mtext("cluster 5",  side=3, at=0.50,line=1,cex=2.5, col="black",outer=TRUE)
mtext("cluster 6",  side=3, at=0.61,line=1,cex=2.5, col="black",outer=TRUE)
mtext("cluster 7",  side=3, at=0.72,line=1,cex=2.5, col="black",outer=TRUE)
mtext("cluster 8",  side=3, at=0.83,line=1,cex=2.5, col="black",outer=TRUE)
mtext("cluster 9",  side=3, at=0.94,line=1,cex=2.5, col="black",outer=TRUE)

## eigenfaces from HCL2 in each cluster
par(mfrow=c(3,K_c),mar=c(0,0,4,0),oma=c(1,4,4,1))
for(i in 1:3){
  for(j in 1:K_c){
    image(t(apply(matrix(ef_HCL2[[j]][,i],nrow=200),2,rev)),main=round(evcumu_HCL2[[j]][i],2),col = grey(seq(0, 1, length = 256)),xaxt="n", yaxt="n",cex.main=3)
  }
}
mtext("3rd PC",  side=2, at=0.12,line=1,cex=2.5, col="black",outer=TRUE)
mtext("2ed PC",  side=2, at=0.46,line=1,cex=2.5, col="black",outer=TRUE)
mtext("1st PC",  side=2, at=0.80,line=1,cex=2.5, col="black",outer=TRUE)
mtext("cluster 1",  side=3, at=0.06,line=1,cex=2.5, col="black",outer=TRUE)
mtext("cluster 2",  side=3, at=0.17,line=1,cex=2.5, col="black",outer=TRUE)
mtext("cluster 3",  side=3, at=0.28,line=1,cex=2.5, col="black",outer=TRUE)
mtext("cluster 4",  side=3, at=0.39,line=1,cex=2.5, col="black",outer=TRUE)
mtext("cluster 5",  side=3, at=0.50,line=1,cex=2.5, col="black",outer=TRUE)
mtext("cluster 6",  side=3, at=0.61,line=1,cex=2.5, col="black",outer=TRUE)
mtext("cluster 7",  side=3, at=0.72,line=1,cex=2.5, col="black",outer=TRUE)
mtext("cluster 8",  side=3, at=0.83,line=1,cex=2.5, col="black",outer=TRUE)
mtext("cluster 9",  side=3, at=0.94,line=1,cex=2.5, col="black",outer=TRUE)

## eigenfaces from LPCA in each cluster
par(mfrow=c(3,K_c),mar=c(0,0,4,0),oma=c(1,4,4,1))
for(i in 1:3){
  for(j in 1:K_c){
    image(t(apply(matrix(ef_LPCA[[j]][,i],nrow=200),2,rev)),main=round(evcumu_LPCA[[j]][i],2),col = grey(seq(0, 1, length = 256)),xaxt="n", yaxt="n",cex.main=3)
  }
}
mtext("3rd PC",  side=2, at=0.12,line=1,cex=2.5, col="black",outer=TRUE)
mtext("2ed PC",  side=2, at=0.46,line=1,cex=2.5, col="black",outer=TRUE)
mtext("1st PC",  side=2, at=0.80,line=1,cex=2.5, col="black",outer=TRUE)
mtext("cluster 1",  side=3, at=0.06,line=1,cex=2.5, col="black",outer=TRUE)
mtext("cluster 2",  side=3, at=0.17,line=1,cex=2.5, col="black",outer=TRUE)
mtext("cluster 3",  side=3, at=0.28,line=1,cex=2.5, col="black",outer=TRUE)
mtext("cluster 4",  side=3, at=0.39,line=1,cex=2.5, col="black",outer=TRUE)
mtext("cluster 5",  side=3, at=0.50,line=1,cex=2.5, col="black",outer=TRUE)
mtext("cluster 6",  side=3, at=0.61,line=1,cex=2.5, col="black",outer=TRUE)
mtext("cluster 7",  side=3, at=0.72,line=1,cex=2.5, col="black",outer=TRUE)
mtext("cluster 8",  side=3, at=0.83,line=1,cex=2.5, col="black",outer=TRUE)
mtext("cluster 9",  side=3, at=0.94,line=1,cex=2.5, col="black",outer=TRUE)

# #################################################
# individual example plot
# #################################################
idx_l = which((e_rep_HCGR-e_rep_HCL2) == min(e_rep_HCGR-e_rep_HCL2))
idx_l = rbind(idx_l,which((e_rep_HCGR-e_rep_LPCA) == min(e_rep_HCGR-e_rep_LPCA)))
idx_l = rbind(idx_l,which((e_rep_HCGR-e_rep_HCG) == min(e_rep_HCGR-e_rep_HCG)))
idx_l = rbind(idx_l,which((e_rep_HCG-e_rep_HCGR) == min(e_rep_HCG-e_rep_HCGR)))

num_ev = trueDim

for(i in 1:length(idx_l)){
  idx = idx_l[i]
  par(mfrow=c(1,6),oma=c(1,1,1,1),mar=c(1,1,4,1))

  ex = as.vector(X.pcscore[idx,1:v.num999]%*%(t(X.eigenvec)[1:v.num999,]) + X.mean)
  image(t(apply(matrix(X[idx,],nrow=200),2,rev)),main='Original',col = grey(seq(0, 1, length = 256)),xaxt="n", yaxt="n",cex.main=4)
  mean((ex-X[idx,])^2)/var(ex)

  ex_gpca = as.vector(X.pcscore[idx,1:num_ev]%*%(t(X.eigenvec)[1:num_ev,]) + X.mean)
  image(t(apply(matrix(ex_gpca,nrow=200),2,rev)),main='GPCA',col = grey(seq(0, 1, length = 256)),xaxt="n", yaxt="n",cex.main=4)
  e_gpca = mean((ex_gpca-ex)^2)/var(ex)

  cluster = indi_HCGR[idx]
  idx_t = which(c(1:n)[indi_HCGR==cluster]==idx)
  ex_HCGR = as.vector(pc_HCGR[[cluster]][idx_t,1:num_ev]%*%(t(ef_HCGR[[cluster]])[1:num_ev,]) + mean_HCGR[[cluster]])
  e_HCGR = mean((ex-ex_HCGR)^2)/var(ex)
  image(t(apply(matrix(ex_HCGR,nrow=200),2,rev)),main="HC-GR",col = grey(seq(0, 1, length = 256)),xaxt="n", yaxt="n",cex.main=4)

  cluster = indi_HCG[idx]
  idx_t = which(c(1:n)[indi_HCG==cluster]==idx)
  ex_HCG = as.vector(pc_HCG[[cluster]][idx_t,1:num_ev]%*%(t(ef_HCG[[cluster]])[1:num_ev,]) + mean_HCG[[cluster]])
  e_HCG = mean((ex-ex_HCG)^2)/var(ex)
  image(t(apply(matrix(ex_HCG,nrow=200),2,rev)),main="HC-G",col = grey(seq(0, 1, length = 256)),xaxt="n", yaxt="n",cex.main=4)

  cluster = indi_HCL2[idx]
  idx_t = which(c(1:n)[indi_HCL2==cluster]==idx)
  ex_HCL2 = as.vector(pc_HCL2[[cluster]][idx_t,1:num_ev]%*%(t(ef_HCL2[[cluster]])[1:num_ev,]) + mean_HCL2[[cluster]])
  e_HCL2 = mean((ex-ex_HCL2)^2)/var(ex)
  image(t(apply(matrix(ex_HCL2,nrow=200),2,rev)),main="HC-L2",col = grey(seq(0, 1, length = 256)),xaxt="n", yaxt="n",cex.main=4)

  cluster = indi_LPCA[idx]
  idx_t = which(c(1:n)[indi_LPCA==cluster]==idx)
  ex_LPCA = as.vector(pc_LPCA[[cluster]][idx_t,1:num_ev]%*%(t(ef_LPCA[[cluster]])[1:num_ev,]) + mean_LPCA[[cluster]])
  e_LPCA = mean((ex-ex_LPCA)^2)/var(ex)
  image(t(apply(matrix(ex_LPCA,nrow=200),2,rev)),main="LPCA",col = grey(seq(0, 1, length = 256)),xaxt="n", yaxt="n",cex.main=4)

  xtable(cbind(e_gpca,e_HCGR,e_HCG,e_HCL2,e_LPCA))
}
