############################## example I: V-shape
rm(list = ls())
## package for pam() kmeans clustering
library(cluster)
library(igraph) ## used by GeoRatio_graph
library(GeoRatio)

# ## load data
# data(V_shape)

## SD of noise
nsd = 0.003

## number of data points
n=200
set.seed(100)

x = runif(n,min=-0.5,max=0.5)
i.f = which(x>0)
y = rep(0,n)
y[i.f] = 0.2*(x[i.f])
x[i.f] = -x[i.f]

V_shape_noiseless = cbind(-x,y)

x = x+rnorm(n,mean=0,sd=nsd)
y = y+rnorm(n,mean=0,sd=nsd)
x = -x

V_shape = cbind(x,y)


## intrinsic dimension of the data
trueDim = 1

############## correlation dimension estimation
## number of epsilons
num = 50

## pair-wise Euclidean distance
dis = dist(V_shape,method = "euclidean", diag = FALSE, upper = TRUE)

## correlation dimension estimation
est = correDimEst(num,dis)

## plot correlation dimension estimation
par(mfrow=c(1,1))
plot(est$x[2:(length(est$x)-1)],est$deri,type="l",ylim=c(0,10),xlab="log(epsilon)",ylab="est. correlation dimension",
     main="Estimated correlation dimension")
abline(trueDim,0,lty=3,col="gray60",lwd=2)

##### decide k_M (plateau log(epsi) from -4.8 to -2.8)
epsi_s = -4.5
epsi_e = -2.5
epsi_seq = seq(epsi_s,epsi_e,length=5)

dis_L2 = as.matrix(dist(V_shape))

for(i in 1:length(epsi_seq)){
  epsi = epsi_seq[i]
  dist_thre = exp(epsi)
  epsi_nbhd = apply(dis_L2,1,function(x) sum(x<dist_thre))
  avg_epsi_nbhd = mean(epsi_nbhd)
  print(epsi)
  print(avg_epsi_nbhd)
}

## parameters for GeoRatio_graph()
K = 2
k_M = 15
k_m = 3
pca_thre = 0.9

set.seed(100)

V_GR = GeoRatio_graph(V_shape,k_max=k_M,k_min=k_m,trueDim,thresh=pca_thre,gmode="max")

dis_L2 = dist(V_shape)
dis_G = as.dist(V_GR$dis_G)
dis_GR = as.dist(V_GR$dis_GR)

hc_L2 = hclust(dis_L2,method="average")
hc_G = hclust(dis_G,method="average")
hc_GR = hclust(dis_GR,method="average")

indi_KML2 = pam(dis_L2,K,diss=T)$clustering
indi_KMG = pam(dis_G,K,diss=T)$clustering
indi_KMGR = pam(dis_GR,K,diss=T)$clustering
indi_HL2 = cutree(hc_L2,k=K)
indi_HG = cutree(hc_G,k=K)
indi_HGR = cutree(hc_GR,k=K)

lpca_V = kplanes(V_shape,K,trueDim,iter.max=100,thresh=1e-5)
indi_lpca = lpca_V$indicator

par(mfrow=c(3,3))
plot(V_shape,col=indi_HGR,main="HC-GR")
plot(V_shape,col=indi_HG,main="HC-G")
plot(V_shape,col=indi_HL2,main="HC-L2")
plot(V_shape,col=indi_KMGR,main="KMC-GR")
plot(V_shape,col=indi_KMG,main="KMC-G")
plot(V_shape,col=indi_KML2,main="KMC-L2")
plot(V_shape,col=indi_lpca,main="LPCA")



####################################################################
##### This section gives the sil info and rep. error for kmeans and
##### hierarchical clustering under different distance measure
##### for K=1 to 10 clusters. Clearly GeoRatio has an advantage 
##### especially under hierarchical clustering
## (reconstruction error using noise data)
#####################################################################
################################# kmeans on L2
sil_KMCL2 = rep(0,10)
e_KMCL2 = rep(0,10)
e_sd_KMCL2 = rep(0,10)
e_noiseless_KMCL2 = rep(0,10)
e_sd_noiseless_KMCL2 = rep(0,10)
for(i in 1:10){
  print(i)
  pamm = pam(V_shape,i)
  if(i==1) {sil_KMCL2[i]=0}
  else
  {sil_KMCL2[i] = pamm$silinfo$avg.width}
  indi = pamm$clustering
  now = lpca(indi,V_shape,d=trueDim)
  
  X.rep = now$X.rep
  error_noiseless.vector = apply((V_shape_noiseless-X.rep)^2,1,mean)/apply(V_shape_noiseless,1,var)
  e_noiseless_KMCL2[i] = mean(error_noiseless.vector)
  e_sd_noiseless_KMCL2[i] = sd(error_noiseless.vector)

  e_KMCL2[i]=now$mean_error
  e_sd_KMCL2[i]=sd(now$all_error)
  print(i)
}
# dev.new()
par(mfrow=c(1,2))
plot(1:10,sil_KMCL2,type="b",main="silhouette (k-means on L2)")
plot(1:10,e_KMCL2,type="b",main="Average representation error")

################################## hierarchical on L2
h_L2 = hclust(as.dist(dis_L2),method="average")

sil_HCL2 = rep(0,10)
e_HCL2 = rep(0,10)
e_sd_HCL2 = rep(0,10)
e_noiseless_HCL2 = rep(0,10)
e_sd_noiseless_HCL2 = rep(0,10)
for(i in 1:10){
  print(i)
  indi = cutree(h_L2,i)
  if(i==1) {sil_HCL2[i]=0}
  else
  {sil_HCL2[i] = summary(silhouette(indi,dist=dis_L2))$avg.width}
  
  now = lpca(indi,V_shape,d=trueDim)
  
  X.rep = now$X.rep
  error_noiseless.vector = apply((V_shape_noiseless-X.rep)^2,1,mean)/apply(V_shape_noiseless,1,var)
  e_noiseless_HCL2[i] = mean(error_noiseless.vector)
  e_sd_noiseless_HCL2[i] = sd(error_noiseless.vector)
  
  e_HCL2[i]=now$mean_error
  e_sd_HCL2[i]=sd(now$all_error)
  print(i)
}
# dev.new()
par(mfrow=c(1,2))
plot(1:10,sil_HCL2,type="b",main="silhouette (hierarchical on L2)")
plot(1:10,e_HCL2,type="b",main="Average representation error")


################################# kmeans on geodesic
sil_KMCG = rep(0,10)
e_KMCG = rep(0,10)
e_sd_KMCG = rep(0,10)
e_noiseless_KMCG = rep(0,10)
e_sd_noiseless_KMCG = rep(0,10)
for(i in 1:10){
  print(i)
  pamm = pam(dis_G,i,diss=T)
  if(i==1) {sil_KMCG[i]=0}
  else
  {sil_KMCG[i] = pamm$silinfo$avg.width}
  indi = pamm$clustering
  now = lpca(indi,V_shape,d=trueDim)
  
  X.rep = now$X.rep
  error_noiseless.vector = apply((V_shape_noiseless-X.rep)^2,1,mean)/apply(V_shape_noiseless,1,var)
  e_noiseless_KMCG[i] = mean(error_noiseless.vector)
  e_sd_noiseless_KMCG[i] = sd(error_noiseless.vector)
  
  e_KMCG[i]=now$mean_error
  e_sd_KMCG[i]=sd(now$all_error)
  print(i)
}
# dev.new()
par(mfrow=c(1,2))
plot(1:10,sil_KMCG,type="b",main="silhouette (k-means on geodesic)")
plot(1:10,e_KMCG,type="b",main="Average representation error")

################################## hierarchical on geodesic
h_G = hclust(as.dist(dis_G),method="average")

sil_HCG = rep(0,10)
e_HCG = rep(0,10)
e_sd_HCG = rep(0,10)
e_noiseless_HCG = rep(0,10)
e_sd_noiseless_HCG = rep(0,10)
for(i in 1:10){
  print(i)
  indi = cutree(h_G,i)
  if(i==1) {sil_HCG[i]=0}
  else
  {sil_HCG[i] = summary(silhouette(indi,dist=dis_G))$avg.width}
  
  now = lpca(indi,V_shape,d=trueDim)
  
  X.rep = now$X.rep
  error_noiseless.vector = apply((V_shape_noiseless-X.rep)^2,1,mean)/apply(V_shape_noiseless,1,var)
  e_noiseless_HCG[i] = mean(error_noiseless.vector)
  e_sd_noiseless_HCG[i] = sd(error_noiseless.vector)
  
  e_HCG[i]=now$mean_error
  e_sd_HCG[i]=sd(now$all_error)
  print(i)
}
# dev.new()
par(mfrow=c(1,2))
plot(1:10,sil_HCG,type="b",main="silhouette (hierarchical on geodesic)")
plot(1:10,e_HCG,type="b",main="Average representation error")

################################# kmeans on GeoRatio
sil_KMCGR = rep(0,10)
e_KMCGR = rep(0,10)
e_sd_KMCGR = rep(0,10)
e_noiseless_KMCGR = rep(0,10)
e_sd_noiseless_KMCGR = rep(0,10)
for(i in 1:10){
  print(i)
  pamm = pam(dis_GR,i,diss=T)
  if(i==1) {sil_KMCGR[i]=0}
  else
  {sil_KMCGR[i] = pamm$silinfo$avg.width}
  indi = pamm$clustering
  now = lpca(indi,V_shape,d=trueDim)
  
  X.rep = now$X.rep
  error_noiseless.vector = apply((V_shape_noiseless-X.rep)^2,1,mean)/apply(V_shape_noiseless,1,var)
  e_noiseless_KMCGR[i] = mean(error_noiseless.vector)
  e_sd_noiseless_KMCGR[i] = sd(error_noiseless.vector)
  
  e_KMCGR[i]=now$mean_error
  e_sd_KMCGR[i]=sd(now$all_error)
  print(i)
}
# dev.new()
par(mfrow=c(1,2))
plot(1:10,sil_KMCGR,type="b",main="silhouette (k-means on GeoRatio)")
plot(1:10,e_KMCGR,type="b",main="Average representation error")

################################## hierarchical on GeoRatio
h_GR = hclust(as.dist(dis_GR),method="average")

sil_HCGR = rep(0,10)
e_HCGR = rep(0,10)
e_sd_HCGR = rep(0,10)
e_noiseless_HCGR = rep(0,10)
e_sd_noiseless_HCGR = rep(0,10)
for(i in 1:10){
  print(i)
  indi = cutree(h_GR,i)
  if(i==1) {sil_HCGR[i]=0}
  else
  {sil_HCGR[i] = summary(silhouette(indi,dist=dis_GR))$avg.width}
  
  now = lpca(indi,V_shape,d=trueDim)
  
  X.rep = now$X.rep
  error_noiseless.vector = apply((V_shape_noiseless-X.rep)^2,1,mean)/apply(V_shape_noiseless,1,var)
  e_noiseless_HCGR[i] = mean(error_noiseless.vector)
  e_sd_noiseless_HCGR[i] = sd(error_noiseless.vector)
  
  e_HCGR[i]=now$mean_error
  e_sd_HCGR[i]=sd(now$all_error)
  
  print(i)
}
# dev.new()
par(mfrow=c(1,2))
plot(1:10,sil_HCGR,type="b",main="silhouette (hierarchical on GeoRatio)")
plot(1:10,e_HCGR,type="b",main="Average representation error")


###################### lpca
set.seed(200)
sil_LPCA = rep(0,10)
e_LPCA = rep(0,10)
e_sd_LPCA = rep(0,10)
e_noiseless_LPCA = rep(0,10)
e_sd_noiseless_LPCA = rep(0,10)
for(i in 1:10){
  print(i)
  KM_test = kplanes(V_shape,i,trueDim,iter.max=200,thresh = 1e-4)
  indi = KM_test$indicator
  if(i==1){
    sil_LPCA[i] = 0
  }else{
    sil_LPCA[i] = summary(silhouette(indi,dist=as.dist(dis_L2)))$avg.width
  }
  
  e_LPCA[i] = KM_test$error
  error.vector = apply((V_shape-KM_test$X.rep)^2,1,mean)/apply(V_shape,1,var)
  error_noiseless.vector = apply((V_shape_noiseless-KM_test$X.rep)^2,1,mean)/apply(V_shape_noiseless,1,var)
  
  e_sd_LPCA[i] = sd(error.vector)
  e_noiseless_LPCA[i] = mean(error_noiseless.vector)
  e_sd_noiseless_LPCA[i] = sd(error_noiseless.vector)
  print(i)
}
# dev.new()
par(mfrow=c(1,2))
plot(1:10,sil_LPCA,type="b",main="silhouette (LPCA)")
plot(1:10,e_LPCA,type="b",main="Average representation error")

