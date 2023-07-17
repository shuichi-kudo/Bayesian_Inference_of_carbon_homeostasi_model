# Plot the histogram of the posterior distribution of the wild type

# settings
gen = "wt"
batches = c(1,2,3,4,5,6,7)
params = c("a","r","h_LD","h_SD","tL_LD","tL_SD","bmax")
xmins = c(100,0.4,80,40,0.6,0.2,0)
xmaxs = c(270,1.0,180,140,0.9,0.5,0.2)
ymaxs = c(0.12,45,0.2,0.4,80,80,120)

# read library
library(RColorBrewer)
cols_original = brewer.pal(length(batches)+3,"Spectral")
cols = cols_original[c(2,3,4,5,7,8,9,10)]

# plot
for(i in 1:length(params)){
  par = params[i]
  xmin = xmins[i]
  xmax = xmaxs[i]
  ymax = ymaxs[i]
  par(pin = c(2.5,2))
  plot(NA,NA,xlim=c(xmin,xmax),ylim=c(0,ymax),xlab="value",ylab="density",main=par,cex.axis=0.9)
  for(j in 1:length(batches)){
    batch = batches[j]
    # read data
    dataset = read.csv(paste("result/parameter_sample/sample_",batch,"_wt.csv",sep=""))
    sample = as.numeric(dataset[par][,1])
    # density
    den = density(sample)
    den_x = den$x
    den_y = den$y
    # plot
    rgb = col2rgb(cols[j])
    mycol <- rgb(rgb[1], rgb[2], rgb[3], max = 255, alpha = 100)
    lines(den_x,den_y,col=cols[j])
    polygon(den_x,den_y,col=mycol)
  }
}

# tL of both photoperiods ######################################################
xmin = 0
xmax = 24
ymax = 2.8

par(pin = c(4,2))
plot(NA,NA,xlim=c(xmin,xmax),ylim=c(0,ymax),xlab="value",ylab="density",main="tL",cex.axis=0.9,xaxt="n")
axis(side=1,at=seq(xmin,xmax,by=4),label=seq(0,24,4))
for(j in 1:length(batches)){
  batch = batches[j]
  # read data
  dataset = read.csv(paste("result/parameter_sample/sample_",batch,"_wt.csv",sep=""))
  sample1 = as.numeric(dataset["tL_LD"][,1])*24
  sample2 = as.numeric(dataset["tL_SD"][,1])*24
  # density
  den1 = density(sample1)
  den1_x = den1$x
  den1_y = den1$y
  den2 = density(sample2)
  den2_x = den2$x
  den2_y = den2$y
  # plot
  rgb = col2rgb(cols[j])
  mycol <- rgb(rgb[1], rgb[2], rgb[3], max = 255, alpha = 100)
  lines(den1_x,den1_y,col=cols[j])
  polygon(den1_x,den1_y,col=mycol)
  lines(den2_x,den2_y,col=cols[j])
  polygon(den2_x,den2_y,col=mycol)

}
