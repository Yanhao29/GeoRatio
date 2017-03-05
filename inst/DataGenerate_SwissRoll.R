##############################
## Swiss roll
##############################
library(rgl)  ## 3d plot

save_path = paste0(getwd(),"/data/")

## SD of noise
nsd = 0.1

## number of data points
n=2000
set.seed(100)

p = runif(n,min=0,max=1)
q = runif(n,min=0,max=1)
t = 3*pi/2*(1+2*p)
x = t*cos(t)+ rnorm(n,0,nsd)
y = t*sin(t)+ rnorm(n,0,nsd)
z = 30*q+ rnorm(n,0,nsd)

SwissRoll = cbind(x,y,z)

rbPal <- colorRampPalette(c('red','yellow'))
col.plot <- rbPal(20)[as.numeric(cut(t,breaks = 20))]
open3d()
plot3d(SwissRoll,col=col.plot,axes=TRUE,box=FALSE,xlab = "", ylab = "", zlab = "")

save(SwissRoll, file=paste0(save_path,"SwissRoll.rda"))





