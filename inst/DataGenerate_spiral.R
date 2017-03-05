##############################
## Spiral
##############################

save_path = paste0(getwd(),"/data/")

## SD of noise
nsd = 0.005

## number of data points
n = 1000
set.seed(100)


t = seq(0,1,length=n)
x = sqrt(t)*cos(10*pi*sqrt(t))+rnorm(n,0,nsd)
y = sqrt(t)*sin(10*pi*sqrt(t))+rnorm(n,0,nsd)

spiral = cbind(x,y)
plot(spiral,xlab='',ylab='')

save(spiral, file=paste0(save_path,"spiral.rda"))


