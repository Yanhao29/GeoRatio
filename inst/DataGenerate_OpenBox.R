##############################
## 3d box
##############################
library(rgl)  ## 3d plot

save_path = paste0(getwd(),"/data/")

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

data = cbind(x,y,z)
col.f = rep(1:6,each=n/6)
open3d()
# plot3d(x,y,z,col=col.f,axes=TRUE,box=FALSE,xlab = "", ylab = "", zlab = "",xlim=c(0,2),ylim=c(0,2),zlim=c(0,2),cex.lab=1.5)
plot3d(data,col=col.f,axes=FALSE,box=FALSE,xlab = "", ylab = "", zlab = "",xlim=c(0,2),ylim=c(0,2),zlim=c(0,2),cex.lab=1)

save(data, file=paste0(save_path,"OpenBox.Rdata"))
