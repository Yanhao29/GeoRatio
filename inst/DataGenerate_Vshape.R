##############################
## V-shape
##############################

save_path = paste0(getwd(),"/data/")

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

x = x+rnorm(n,mean=0,sd=nsd)
y = y+rnorm(n,mean=0,sd=nsd)
x = -x

V_shape = cbind(x,y)

#col.f = rep(1:2,each=n/2)
col.f = rep(1,n)
col.f[i.f]=2

par(mfrow=c(1,1))
## plot true data
plot(x,y,col=col.f,xlab = "", ylab = "",xlim=c(0,0.5),ylim=c(0,0.5))


save(V_shape, file=paste0(save_path,"V_shape.rda"))

