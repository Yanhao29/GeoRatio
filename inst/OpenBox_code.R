############################## example III: Open Box
rm(list = ls())
## package for pam() kmeans clustering
library(cluster)
library(igraph) ## used by GeoRatio_graph
library(GeoRatio)
library(xtable)
library(rgl)

# ## load data
# data(OpenBox)

## SD of noise
nsd = 0.01

## number of data points
n=1200
set.seed(100)

x = rep(0,n)
y = rep(0,n)
z = rep(0,n)
x[1:(n/6)] = runif(n/6)
y[1:(n/6)] = runif(n/6)
z[1:(n/6)] = 0 + rnorm(n/6,0,nsd)
x[(n/6+1):(n/3)] = 1 + rnorm(n/6,0,nsd)
y[(n/6+1):(n/3)] = runif(n/6)
z[(n/6+1):(n/3)] = runif(n/6)
x[(n/3+1):(n/2)] = runif(n/6)
y[(n/3+1):(n/2)] = 1 + rnorm(n/6,0,nsd)
z[(n/3+1):(n/2)] = runif(n/6)
x[(n/2+1):(2*n/3)] = 0 + rnorm(n/6,0,nsd)
y[(n/2+1):(2*n/3)] = runif(n/6)
z[(n/2+1):(2*n/3)] = runif(n/6)
x[(2*n/3+1):(5*n/6)] = runif(n/6)
y[(2*n/3+1):(5*n/6)] = 0 + rnorm(n/6,0,nsd)
z[(2*n/3+1):(5*n/6)] = runif(n/6)
x[(5*n/6+1):n] = runif(n/6,min=1,max=2)
y[(5*n/6+1):n] = runif(n/6)
z[(5*n/6+1):n] = 1 + rnorm(n/6,0,nsd)

xx = rep(0,n)
yy = rep(0,n)
zz = rep(0,n)
xx[1:(n/6)] = runif(n/6)
yy[1:(n/6)] = runif(n/6)
zz[1:(n/6)] = 0 
xx[(n/6+1):(n/3)] = 1 
yy[(n/6+1):(n/3)] = runif(n/6)
zz[(n/6+1):(n/3)] = runif(n/6)
xx[(n/3+1):(n/2)] = runif(n/6)
yy[(n/3+1):(n/2)] = 1 
zz[(n/3+1):(n/2)] = runif(n/6)
xx[(n/2+1):(2*n/3)] = 0 
yy[(n/2+1):(2*n/3)] = runif(n/6)
zz[(n/2+1):(2*n/3)] = runif(n/6)
xx[(2*n/3+1):(5*n/6)] = runif(n/6)
yy[(2*n/3+1):(5*n/6)] = 0 
zz[(2*n/3+1):(5*n/6)] = runif(n/6)
xx[(5*n/6+1):n] = runif(n/6,min=1,max=2)
yy[(5*n/6+1):n] = runif(n/6)
zz[(5*n/6+1):n] = 1 

data = cbind(x,y,z)
data_nl = cbind(xx,yy,zz)

col.f = rep(1:6,each=n/6)
open3d()
plot3d(data,col=col.f,axes=TRUE,box=FALSE,xlab = "", ylab = "", zlab = "",xlim=c(0,2),ylim=c(0,2),zlim=c(0,2),cex.lab=1.5)
open3d()
plot3d(data_nl,col=col.f,axes=TRUE,box=FALSE,xlab = "", ylab = "", zlab = "",xlim=c(0,2),ylim=c(0,2),zlim=c(0,2),cex.lab=1)


## intrinsic dimension of the data
trueDim = 2

## cluster number
K = 6

############## correlation dimension estimation
## number of epsilons
num = 50

## pair-wise Euclidean distance
dis = dist(data,method = "euclidean", diag = FALSE, upper = TRUE)

## correlation dimension estimation
est = correDimEst(num,dis)

## plot correlation dimension estimation
par(mfrow=c(1,1))
plot(est$x[2:(length(est$x)-1)],est$deri,type="l",ylim=c(0,10),xlab="log(epsilon)",ylab="est. correlation dimension",
     main="Estimated correlation dimension")
abline(trueDim,0,lty=3,col="gray60",lwd=2)

##### decide k_M (plateau log(epsi) from -4.8 to -2.8)
epsi_s = -3.5
epsi_e = -1
epsi_seq = seq(epsi_s,epsi_e,length=5)

dis_L2 = as.matrix(dist(data))

for(i in 1:length(epsi_seq)){
  epsi = epsi_seq[i]
  dist_thre = exp(epsi)
  epsi_nbhd = apply(dis_L2,1,function(x) sum(x<dist_thre))
  avg_epsi_nbhd = mean(epsi_nbhd)
  print(epsi)
  print(avg_epsi_nbhd)
}

## parameters for GeoRatio_graph()
k_M = 40
k_m = 3
pca_thre = 0.9

set.seed(100)

data_GR = GeoRatio_graph(data,k_max=k_M,k_min=k_m,trueDim,thresh=pca_thre,gmode="max")

dis_L2 = dist(data)
dis_G = as.dist(data_GR$dis_G)
dis_GR = as.dist(data_GR$dis_GR)

hc_L2 = hclust(dis_L2,method="average")
hc_G = hclust(dis_G,method="average")
hc_GR = hclust(dis_GR,method="average")

indi_KML2 = pam(dis_L2,K,diss=T)$clustering
now_KMCL2 = lpca(indi_KML2, data, d=trueDim)
X_KMCL2 = now_KMCL2$X.rep

indi_KMG = pam(dis_G,K,diss=T)$clustering
now_KMCG = lpca(indi_KMG, data, d=trueDim)
X_KMCG = now_KMCG$X.rep

indi_KMGR = pam(dis_GR,K,diss=T)$clustering
now_KMCGR = lpca(indi_KMGR, data, d=trueDim)
X_KMCGR = now_KMCGR$X.rep

indi_HL2 = cutree(hc_L2,k=K)
now_HCL2 = lpca(indi_HL2, data, d=trueDim)
X_HCL2 = now_HCL2$X.rep

indi_HG = cutree(hc_G,k=K)
now_HCG = lpca(indi_HG, data, d=trueDim)
X_HCG = now_HCG$X.rep

indi_HGR = cutree(hc_GR,k=K)
now_HCGR = lpca(indi_HGR, data, d=trueDim)
X_HCGR = now_HCGR$X.rep

lpca_data = kplanes(data,K,trueDim,iter.max=100,thresh=1e-5)
indi_lpca = lpca_data$indicator
now_LPCA = lpca(indi_lpca, data, d=trueDim)
X_LPCA = now_LPCA$X.rep

open3d()
plot3d(data,col=indi_HGR,main="HC-GR")
open3d()
plot3d(data,col=indi_HG,main="HC-G")
open3d()
plot3d(data,col=indi_HL2,main="HC-L2")
oepn3d()
plot3d(data,col=indi_KMGR,main="KMC-GR")
open3d()
plot3d(data,col=indi_KMG,main="KMC-G")
open3d()
plot3d(data,col=indi_KML2,main="KMC-L2")
open3d()
plot3d(data,col=indi_lpca,main="LPCA")

open3d()
plot3d(X_HCGR,col=indi_HGR,main="HC-GR")
open3d()
plot3d(X_HCG,col=indi_HG,main="HC-G")
open3d()
plot3d(X_HCL2,col=indi_HL2,main="HC-L2")
open3d()
plot3d(X_KMCGR,col=indi_KMGR,main="KMC-GR")
open3d()
plot3d(X_KMCG,col=indi_KMG,main="KMC-G")
open3d()
plot3d(X_KMCL2,col=indi_KML2,main="KMC-L2")
open3d()
plot3d(X_LPCA,col=indi_lpca,main="LPCA")

# data_nl[1,]
# data[1,]
# X_LPCA[1,]
# mean((X_LPCA[1,]-data[1,])^2)/var(data[1,])
# mean((X_LPCA[1,]-data_nl[1,])^2)/var(data_nl[1,])

var_data = apply(data,1,var)
var_data_nl = apply(data_nl,1,var)

re_HCL2 = apply((X_HCL2-data)^2,1,mean)
re_HCL2_nl = apply((X_HCL2-data_nl)^2,1,mean)

re_KMCL2 = apply((X_KMCL2-data)^2,1,mean)
re_KMCL2_nl = apply((X_KMCL2-data_nl)^2,1,mean)

re_HCG = apply((X_HCG-data)^2,1,mean)
re_HCG_nl = apply((X_HCG-data_nl)^2,1,mean)

re_KMCG = apply((X_KMCG-data)^2,1,mean)
re_KMCG_nl = apply((X_KMCG-data)^2,1,mean)

re_HCGR = apply((X_HCGR-data)^2,1,mean)
re_HCGR_nl = apply((X_HCGR-data_nl)^2,1,mean)

re_KMCGR = apply((X_KMCGR-data)^2,1,mean)
re_KMCGR_nl = apply((X_KMCGR-data_nl)^2,1,mean)

re_LPCA = apply((X_LPCA-data)^2,1,mean)
re_LPCA_nl = apply((X_LPCA-data_nl)^2,1,mean)

nre_HCL2 = apply((X_HCL2-data)^2,1,mean)/apply(data,1,var)
nre_HCL2_nl = apply((X_HCL2-data_nl)^2,1,mean)/apply(data_nl,1,var)
nre_KMCL2 = apply((X_KMCL2-data)^2,1,mean)/apply(data,1,var)
nre_KMCL2_nl = apply((X_KMCL2-data_nl)^2,1,mean)/apply(data_nl,1,var)
nre_HCG = apply((X_HCG-data)^2,1,mean)/apply(data,1,var)
nre_HCG_nl = apply((X_HCG-data_nl)^2,1,mean)/apply(data_nl,1,var)
nre_KMCG = apply((X_KMCG-data)^2,1,mean)/apply(data,1,var)
nre_KMCG_nl = apply((X_KMCG-data_nl)^2,1,mean)/apply(data_nl,1,var)
nre_HCGR = apply((X_HCGR-data)^2,1,mean)/apply(data,1,var)
nre_HCGR_nl = apply((X_HCGR-data_nl)^2,1,mean)/apply(data_nl,1,var)
nre_KMCGR = apply((X_KMCGR-data)^2,1,mean)/apply(data,1,var)
nre_KMCGR_nl = apply((X_KMCGR-data_nl)^2,1,mean)/apply(data_nl,1,var)
nre_LPCA = apply((X_LPCA-data)^2,1,mean)/apply(data,1,var)
nre_LPCA_nl = apply((X_LPCA-data_nl)^2,1,mean)/apply(data_nl,1,var)

par(mfrow=c(1,1))
plot(var_data,var_data_nl,xlab='data',ylab='noiseless data', main = 'Var')
abline(a=0,b=1,col=2)

plot(re_HCL2,re_HCL2_nl,xlab='data',ylab='noiseless data', main = 'HC-L2 Recon. Error')
abline(a=0,b=1,col=2)
plot(re_KMCL2,re_KMCL2_nl,xlab='data',ylab='noiseless data', main = 'KMC-L2 Recon. Error')
abline(a=0,b=1,col=2)
plot(re_HCG,re_HCG_nl,xlab='data',ylab='noiseless data', main = 'HC-G Recon. Error')
abline(a=0,b=1,col=2)
plot(re_KMCG,re_KMCG_nl,xlab='data',ylab='noiseless data', main = 'KMC-G Recon. Error')
abline(a=0,b=1,col=2)
plot(re_HCGR,re_HCGR_nl,xlab='data',ylab='noiseless data', main = 'HC-GR Recon. Error')
abline(a=0,b=1,col=2)
plot(re_KMCGR,re_KMCGR_nl,xlab='data',ylab='noiseless data', main = 'KMC-GR Recon. Error')
abline(a=0,b=1,col=2)
plot(re_LPCA,re_LPCA_nl,xlab='data',ylab='noiseless data', main = 'LPCA Recon. Error')
abline(a=0,b=1,col=2)

plot(nre_HCL2,nre_HCL2_nl,xlab='data',ylab='noiseless data', main = 'HC-L2 Norm. Recon. Error')
abline(a=0,b=1,col=2)
plot(nre_KMCL2,nre_KMCL2_nl,xlab='data',ylab='noiseless data', main = 'KMC-L2 Norm. Recon. Error')
abline(a=0,b=1,col=2)
plot(nre_HCG,nre_HCG_nl,xlab='data',ylab='noiseless data', main = 'HC-G Norm. Recon. Error')
abline(a=0,b=1,col=2)
plot(nre_KMCG,nre_KMCG_nl,xlab='data',ylab='noiseless data', main = 'KMC-G Norm. Recon. Error')
abline(a=0,b=1,col=2)
plot(nre_HCGR,nre_HCGR_nl,xlab='data',ylab='noiseless data', main = 'HC-GR Norm. Recon. Error')
abline(a=0,b=1,col=2)
plot(nre_KMCGR,nre_KMCGR_nl,xlab='data',ylab='noiseless data', main = 'KMC-GR Norm. Recon. Error')
abline(a=0,b=1,col=2)
plot(nre_LPCA,nre_LPCA_nl,xlab='data',ylab='noiseless data', main = 'LPCA Norm. Recon. Error')
abline(a=0,b=1,col=2)


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
e_KMCL2_nl = rep(0,10)
e_sd_KMCL2_nl = rep(0,10)
for(i in 1:10){
  print(i)
  pamm = pam(data,i)
  if(i==1) {sil_KMCL2[i]=0}
  else
  {sil_KMCL2[i] = pamm$silinfo$avg.width}
  indi = pamm$clustering
  now = lpca(indi,data,d=trueDim)
  
  X.rep = now$X.rep
  error_noiseless.vector = apply((data_nl-X.rep)^2,1,mean)
  error_noiseless.vector = apply((data_nl-X.rep)^2,1,mean)/apply(data_nl,1,var)
  e_KMCL2_nl[i] = mean(error_noiseless.vector)
  e_sd_KMCL2_nl[i] = sd(error_noiseless.vector)
  
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
e_HCL2_nl = rep(0,10)
e_sd_HCL2_nl = rep(0,10)
for(i in 1:10){
  print(i)
  indi = cutree(h_L2,i)
  if(i==1) {sil_HCL2[i]=0}
  else
  {sil_HCL2[i] = summary(silhouette(indi,dist=dis_L2))$avg.width}
  
  now = lpca(indi,data,d=trueDim)
  
  X.rep = now$X.rep
  error_noiseless.vector = apply((data_nl-X.rep)^2,1,mean)
  error_noiseless.vector = apply((data_nl-X.rep)^2,1,mean)/apply(data_nl,1,var)
  e_HCL2_nl[i] = mean(error_noiseless.vector)
  e_sd_HCL2_nl[i] = sd(error_noiseless.vector)
  
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
e_KMCG_nl = rep(0,10)
e_sd_KMCG_nl = rep(0,10)
for(i in 1:10){
  print(i)
  pamm = pam(dis_G,i,diss=T)
  if(i==1) {sil_KMCG[i]=0}
  else
  {sil_KMCG[i] = pamm$silinfo$avg.width}
  indi = pamm$clustering
  now = lpca(indi,data,d=trueDim)
  
  X.rep = now$X.rep
  error_noiseless.vector = apply((data_nl-X.rep)^2,1,mean)
  error_noiseless.vector = apply((data_nl-X.rep)^2,1,mean)/apply(data_nl,1,var)
  e_KMCG_nl[i] = mean(error_noiseless.vector)
  e_sd_KMCG_nl[i] = sd(error_noiseless.vector)
  
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
e_HCG_nl = rep(0,10)
e_sd_HCG_nl = rep(0,10)
for(i in 1:10){
  print(i)
  indi = cutree(h_G,i)
  if(i==1) {sil_HCG[i]=0}
  else
  {sil_HCG[i] = summary(silhouette(indi,dist=dis_G))$avg.width}
  
  now = lpca(indi,data,d=trueDim)
  
  X.rep = now$X.rep
  error_noiseless.vector = apply((data_nl-X.rep)^2,1,mean)
  error_noiseless.vector = apply((data_nl-X.rep)^2,1,mean)/apply(data_nl,1,var)
  e_HCG_nl[i] = mean(error_noiseless.vector)
  e_sd_HCG_nl[i] = sd(error_noiseless.vector)
  
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
e_KMCGR_nl = rep(0,10)
e_sd_KMCGR_nl = rep(0,10)
for(i in 1:10){
  print(i)
  pamm = pam(dis_GR,i,diss=T)
  if(i==1) {sil_KMCGR[i]=0}
  else
  {sil_KMCGR[i] = pamm$silinfo$avg.width}
  indi = pamm$clustering
  now = lpca(indi,data,d=trueDim)
  
  X.rep = now$X.rep
  error_noiseless.vector = apply((data_nl-X.rep)^2,1,mean)
  error_noiseless.vector = apply((data_nl-X.rep)^2,1,mean)/apply(data_nl,1,var)
  e_KMCGR_nl[i] = mean(error_noiseless.vector)
  e_sd_KMCGR_nl[i] = sd(error_noiseless.vector)
  
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
e_HCGR_nl = rep(0,10)
e_sd_HCGR_nl = rep(0,10)
for(i in 1:10){
  print(i)
  indi = cutree(h_GR,i)
  if(i==1) {sil_HCGR[i]=0}
  else
  {sil_HCGR[i] = summary(silhouette(indi,dist=dis_GR))$avg.width}
  
  now = lpca(indi,data,d=trueDim)
  
  X.rep = now$X.rep
  error_noiseless.vector = apply((data_nl-X.rep)^2,1,mean)
  error_noiseless.vector = apply((data_nl-X.rep)^2,1,mean)/apply(data_nl,1,var)
  e_HCGR_nl[i] = mean(error_noiseless.vector)
  e_sd_HCGR_nl[i] = sd(error_noiseless.vector)
  
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
e_LPCA_nl = rep(0,10)
e_sd_LPCA_nl = rep(0,10)
for(i in 1:10){
  print(i)
  KM_test = kplanes(data,i,trueDim,iter.max=200,thresh = 1e-4)
  indi = KM_test$indicator
  if(i==1){
    sil_LPCA[i] = 0
  }else{
    sil_LPCA[i] = summary(silhouette(indi,dist=as.dist(dis_L2)))$avg.width
  }
  
  e_LPCA[i] = KM_test$error
  error.vector = apply((data-KM_test$X.rep)^2,1,mean)/apply(data,1,var)
  error_noiseless.vector = apply((data_nl-KM_test$X.rep)^2,1,mean)/apply(data_nl,1,var)
  
  e_sd_LPCA[i] = sd(error.vector)
  e_LPCA_nl[i] = mean(error_noiseless.vector)
  e_sd_LPCA_nl[i] = sd(error_noiseless.vector)
  print(i)
}
# dev.new()
par(mfrow=c(1,2))
plot(1:10,sil_LPCA,type="b",main="silhouette (LPCA)")
plot(1:10,e_LPCA,type="b",main="Average representation error")


################################# 
## All tables
################################# 
xtable(rbind(e_HCGR,e_KMCGR,e_HCG,e_KMCG,e_HCL2,e_KMCL2,e_LPCA)*100)
xtable(rbind(e_sd_HCGR,e_sd_KMCGR,e_sd_HCG,e_sd_KMCG,e_sd_HCL2,e_sd_KMCL2,e_sd_LPCA)*100)

xtable(rbind(e_HCGR_nl,e_KMCGR_nl,e_HCG_nl,e_KMCG_nl,e_HCL2_nl,e_KMCL2_nl,e_LPCA_nl)*100)
xtable(rbind(e_sd_HCGR_nl,e_sd_KMCGR_nl,e_sd_HCG_nl,e_sd_KMCG_nl,e_sd_HCL2_nl,e_sd_KMCL2_nl,e_sd_LPCA_nl)*100)

xtable(rbind(sil_HCGR,sil_KMCGR,sil_HCG,sil_KMCG,sil_HCL2,sil_KMCL2))

size_axis = 1.5
size_lab = 1.5
size_main = 2

par(oma=c(0,0.4,0,0),mar=c(5,5,2,2),mfrow=c(1,2))
limi = max(rbind(sil_HCGR,sil_HCG,sil_HCL2,sil_KMCGR,sil_KMCG,sil_KMCL2))
plot(1:10,sil_HCGR,type="b",pch=1,lty=1,col=2,lwd=2,ylim=c(-0.1,limi),xlab="number of clusters",ylab="sil. information"
     ,cex.axis=size_axis,cex.lab=size_lab,cex.main=size_main)
points(1:10,sil_HCG,type="b",pch=1,lty=1,col=3)
points(1:10,sil_HCL2,type="b",pch=1,lty=1,col=1)
points(1:10,sil_KMCGR,type="b",pch=2,lty=1,col=2)
points(1:10,sil_KMCG,type="b",pch=2,lty=1,col=3)
points(1:10,sil_KMCL2,type="b",pch=2,lty=1,col=1)
abline(v=K,col='gray60',lty=3,lwd=2)
# legend(6,1.5,c("HC-GR","HC-G","HC-L2","KMC-GR","KMC-G","KMC-L2"),pch=c(1,1,1,2,2,2),col=c(2,3,1,2,3,1),lwd=c(2,1,1,1,1,1))
# dev.off()

# pdf(paste0(save_path,'M_shape',toString(k_M), '_error.pdf'),width=5.5, height=5)
# par(oma=c(0,0.4,0,0),mfrow=c(1,1))
limi = max(rbind(e_HCGR,e_HCG,e_HCL2,e_KMCGR,e_KMCG,e_KMCL2))
plot(1:10,e_HCGR,type="b",pch=1,lty=1,col=2,lwd=2,ylim=c(0,1),xlab="number of clusters",ylab="error"
     ,cex.axis=size_axis,cex.lab=size_lab,cex.main=size_main)
points(1:10,e_HCG,type="b",pch=1,lty=1,col=3)
points(1:10,e_HCL2,type="b",pch=1,lty=1,col=1)
points(1:10,e_KMCGR,type="b",pch=2,lty=1,col=2)
points(1:10,e_KMCG,type="b",pch=2,lty=1,col=3)
points(1:10,e_KMCL2,type="b",pch=2,lty=1,col=1)
points(1:10,e_LPCA,type="b",pch=3,lty=2,col=4)
abline(v=K,col='gray60',lty=3,lwd=2)
legend(6,1,c("HC-GR","HC-G","HC-L2","KMC-GR","KMC-G","KMC-L2","LPCA"),pch=c(1,1,1,2,2,2,3),col=c(2,3,1,2,3,1,4),lwd=c(2,1,1,1,1,1,1),lty=c(1,1,1,1,1,1,2))
