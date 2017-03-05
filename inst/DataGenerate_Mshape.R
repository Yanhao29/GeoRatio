##############################
## M-shape
##############################

save_path = paste0(getwd(),"/data/")

## SD of noise
nsd = 0.003

## number of data points
n=500
set.seed(100)

x = rep(0,n)
y = rep(0,n)
t = runif(n,min=-0.5,max=0)
x[1:(n/4)] = t[1:(n/4)]
x[(n/4+1):(n/2)] = t[(n/4+1):(n/2)]
x[(n/2+1):(3*n/4)] = t[(n/2+1):(3*n/4)]
x[(3*n/4+1):n] = t[(3*n/4+1):n]
y[1:(n/4)] = 0
y[(n/4+1):(n/2)] = 0.1-0.2*x[(n/4+1):(n/2)]
y[(n/2+1):(3*n/4)] = 0.1+0.2*x[(n/2+1):(3*n/4)]
y[(3*n/4+1):n] = 0.2


x = x+rnorm(n,mean=0,sd=nsd)
y = y+rnorm(n,mean=0,sd=nsd)

M_shape = cbind(x,y)
col.f = rep(1:4,each=n/4)
#plot3d(x,y,z,col=col.f,xlab = "x", ylab = "z", zlab="z",xlim=c(-0.5,0),ylim=c(-0.2,0.2),zlim=c(0,0.5))
# dev.new()
par(mfrow=c(1,1))
plot(M_shape,col=col.f,xlab = "x", ylab = "z",xlim=c(0,-0.5),ylim=c(-0.1,0.4))

save(M_shape, file=paste0(save_path,"M_shape.rda"))


