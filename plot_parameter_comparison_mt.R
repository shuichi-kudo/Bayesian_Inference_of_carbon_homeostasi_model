# Plot a heatmap of the parameter changes from the posterior estimate of the mutants
# Plot a scatter graph of tau L of the mutants

# Settings #####################################################################
ps = read.csv("result/posterior_summary.csv")
mutants = c("gwd","pwd","bam1","bam2","amy3-bam1",
            "sweet1112",
            "prr7","prr7prr9","ztl3",
            "rb5","rb22","kin10","flz14-1","bzip63-1")
Nm = length(mutants)

#log2 fold change ##############################################################
pars = c("a","r","h_LD","h_SD","bmax")
Nm = length(mutants)
Np = length(pars)
log2fc = data.frame(matrix(NA,nrow=Nm,ncol=Np))
colnames(log2fc) = pars
rownames(log2fc) = mutants
ps_wt = subset(ps,ps$Genotype=="wt")
for(i in 1:Nm){
  gen = mutants[i]
  ps_mt = subset(ps,ps$Genotype==gen)
  bat = ps_mt$Batch
  ps_wt_b = subset(ps_wt,ps_wt$Batch==bat)
  for(j in 1:Np){
    par_wt = ps_wt_b[2+j]
    par_mt = ps_mt[2+j]
    log2fc[i,j] = log2(par_mt/par_wt)
  }
}

library(gplots)
colors = seq(-4,4,length=100)
my_palette <- colorRampPalette(c("navy","blue", "deepskyblue","white", "orange","red","deeppink"))
heatmap.2(as.matrix(log2fc),col=my_palette,
          Rowv=NA,Colv=NA,
          breaks=colors, density.info="none",
          dendrogram="none",symm=F,symkey=F,symbreaks=T, scale="none",
          margins=c(7,7))
          # 4 x 6 inches
# Scatter plot of tL ###########################################################

#difference
dps = data.frame(matrix(NA,nrow=Nm,ncol=2))
colnames(dps) = c("tL_LD","tL_SD")
rownames(dps) = mutants
for(i in 1:length(mutants)){
  gen = mutants[i]
  ps_mt = subset(ps,ps$Genotype==gen)
  bat = ps_mt$Batch
  ps_wt = subset(ps,(ps$Batch==bat)&(ps$Genotype=="wt"))
  dps[i,1] = (ps_mt$tL_LD - ps_wt$tL_LD)*24
  dps[i,2] = (ps_mt$tL_SD - ps_wt$tL_SD)*24
}

# scatter plot
dps1 = dps[1:5,] # starch degrading enzymes
dps2 = dps[6,] # sucrose transporter
dps3 = dps[7:9,] # clock components
dps4 = dps[10:14,] # sugar sensors
x1 = dps1$tL_LD
y1 = dps1$tL_SD
x2 = dps2$tL_LD
y2 = dps2$tL_SD
x3 = dps3$tL_LD
y3 = dps3$tL_SD
x4 = dps4$tL_LD
y4 = dps4$tL_SD

par(pin = c(2,2))
plot(NA,NA,xlim=c(-3,7),ylim=c(-2,7),xlab="difference of tL (LD)",ylab="difference of tL (SD)",asp=1)
abline(h=0,col="black")
abline(v=0,col="black")
abline(h=c(2,-2,4,-4,6,-6),col="gray",lty=2)
abline(v=c(2,-2,4,-4,6,-6),col="gray",lty=2)
points(x1,y1,col="orange",pch=19,cex=1.5)
points(x3,y3,col="deepskyblue",pch=15,cex=1.5)
points(x4,y4,col="green3",pch=17,cex=1.5)
points(x2,y2,col="pink",pch=18,cex=1.5)



par(pin = c(4,4))
plot(NA,NA,xlim=c(-2.3,1.3),ylim=c(-1,1),xlab="difference of tL (LD)",ylab="difference of tL (SD)",asp=1)
abline(h=0,col="black")
abline(v=0,col="black")
abline(h=c(2,-2,1,-1),col="gray",lty=2)
abline(v=c(2,-2,1,-1),col="gray",lty=2)
points(x1,y1,col="orange",pch=19,cex=1.5)
points(x3,y3,col="deepskyblue",pch=15,cex=1.5)
points(x4,y4,col="green3",pch=17,cex=1.5)
points(x2,y2,col="pink",pch=18,cex=1.5)

# Scatter plot of tD ###########################################################

#difference
dps = data.frame(matrix(NA,nrow=Nm,ncol=2))
colnames(dps) = c("tD_LD","tD_SD")
rownames(dps) = mutants
for(i in 1:length(mutants)){
  gen = mutants[i]
  ps_mt = subset(ps,ps$Genotype==gen)
  bat = ps_mt$Batch
  ps_wt = subset(ps,(ps$Batch==bat)&(ps$Genotype=="wt"))
  dps[i,1] = ((1-ps_mt$tL_LD) - (1-ps_wt$tL_LD))*24
  dps[i,2] = ((1-ps_mt$tL_SD) - (1-ps_wt$tL_SD))*24
}

# scatter plot
dps1 = dps[1:5,] # starch degrading enzymes
dps2 = dps[6,] # sucrose transporter
dps3 = dps[7:9,] # clock components
dps4 = dps[10:14,] # sugar sensors
x1 = dps1$tD_LD
y1 = dps1$tD_SD
x2 = dps2$tD_LD
y2 = dps2$tD_SD
x3 = dps3$tD_LD
y3 = dps3$tD_SD
x4 = dps4$tD_LD
y4 = dps4$tD_SD

par(pin = c(2,2))
plot(NA,NA,xlim=c(-7,3),ylim=c(-7,3),xlab="difference of tD (LD)",ylab="difference of tD (SD)",asp=1)
abline(h=0,col="black")
abline(v=0,col="black")
abline(h=c(2,-2,4,-4,6,-6),col="gray",lty=2)
abline(v=c(2,-2,4,-4,6,-6),col="gray",lty=2)
points(x1,y1,col="orange",pch=19,cex=1.5)
points(x3,y3,col="deepskyblue",pch=15,cex=1.5)
points(x4,y4,col="green3",pch=17,cex=1.5)
points(x2,y2,col="pink",pch=18,cex=1.5)



par(pin = c(4,4))
plot(NA,NA,xlim=c(-1.3,2.3),ylim=c(-1,1),xlab="difference of tD (LD)",ylab="difference of tD (SD)",asp=1)
abline(h=0,col="black")
abline(v=0,col="black")
abline(h=c(2,-2,1,-1),col="gray",lty=2)
abline(v=c(2,-2,1,-1),col="gray",lty=2)
points(x1,y1,col="orange",pch=19,cex=1.5)
points(x3,y3,col="deepskyblue",pch=15,cex=1.5)
points(x4,y4,col="green3",pch=17,cex=1.5)
points(x2,y2,col="pink",pch=18,cex=1.5)


# dtL with error bar ###########################################################
library(HDInterval)

# difference
dps = data.frame(matrix(NA,nrow=Nm,ncol=6))
colnames(dps) = c("tL_LD_mu","tL_LD_low","tL_LD_upp","tL_SD_mu","tL_SD_low","tL_SD_upp")
rownames(dps) = mutants

for(i in 1:length(mutants)){
  #target
  gen = mutants[i]
  ps_mt = subset(ps,ps$Genotype==gen)
  bat = ps_mt$Batch
  
  sample_wt = read.csv(paste("result/parameter_sample/sample_",bat,"_wt.csv",sep=""))
  sample_mt = read.csv(paste("result/parameter_sample/sample_",bat,"_",gen,".csv",sep=""))
  sample_tL_LD_wt0 = sample_wt$tL_LD*24
  sample_tL_SD_wt0 = sample_wt$tL_SD*24
  sample_tL_LD_mt0 = sample_mt$tL_LD*24
  sample_tL_SD_mt0 = sample_mt$tL_SD*24
  
  #sampling
  sample_tL_LD_wt = sample(sample_tL_LD_wt0,1000,replace=TRUE)
  sample_tL_SD_wt = sample(sample_tL_SD_wt0,1000,replace=TRUE)
  sample_tL_LD_mt = sample(sample_tL_LD_mt0,1000,replace=TRUE)
  sample_tL_SD_mt = sample(sample_tL_SD_mt0,1000,replace=TRUE)
  sample_dtL_LD = sample_tL_LD_mt - sample_tL_LD_wt
  sample_dtL_SD = sample_tL_SD_mt - sample_tL_SD_wt
  
  write.csv(sample_tL_LD_wt,paste("result/sample_dtL/sample_dtL_LD_",gen,"_","bat",".csv",sep=""))
  
  #write
  dps[i,1] = mean(sample_dtL_LD)
  dps[i,2] = (hdi(sample_dtL_LD,0.95))[1]
  dps[i,3] = (hdi(sample_dtL_LD,0.95))[2]
  dps[i,4] = mean(sample_dtL_SD)
  dps[i,5] = (hdi(sample_dtL_SD,0.95))[1]
  dps[i,6] = (hdi(sample_dtL_SD,0.95))[2]
}

write.csv(dps,"result/summary_dtL.csv")

# scatter plot
dps1 = dps[1:5,] # starch degrading enzymes
dps2 = dps[6,] # sucrose transporter
dps3 = dps[7:9,] # clock components
dps4 = dps[10:14,] # sugar sensors

px1 = dps1$tL_LD_mu
lx1 = dps1$tL_LD_low
ux1 = dps1$tL_LD_up
py1 = dps1$tL_SD_mu
ly1 = dps1$tL_SD_low
uy1 = dps1$tL_SD_up

px2 = dps2$tL_LD_mu
lx2 = dps2$tL_LD_low
ux2 = dps2$tL_LD_up
py2 = dps2$tL_SD_mu
ly2 = dps2$tL_SD_low
uy2 = dps2$tL_SD_up

px3 = dps3$tL_LD_mu
lx3 = dps3$tL_LD_low
ux3 = dps3$tL_LD_up
py3 = dps3$tL_SD_mu
ly3 = dps3$tL_SD_low
uy3 = dps3$tL_SD_up

px4 = dps4$tL_LD_mu
lx4 = dps4$tL_LD_low
ux4 = dps4$tL_LD_up
py4 = dps4$tL_SD_mu
ly4 = dps4$tL_SD_low
uy4 = dps4$tL_SD_up


par(pin = c(2,2))
plot(NA,NA,xlim=c(-3,7.2),ylim=c(-2,7.2),xlab="difference of tL (LD)",ylab="difference of tL (SD)",asp=1)
abline(h=0,col="black")
abline(v=0,col="black")
abline(h=c(2,-2,4,-4,6,-6),col="gray",lty=2)
abline(v=c(2,-2,4,-4,6,-6),col="gray",lty=2)
points(px1,py1,col="orange",pch=19,cex=1.5)
points(px3,py3,col="deepskyblue",pch=15,cex=1.5)
points(px4,py4,col="green3",pch=17,cex=1.5)
points(px2,py2,col="hotpink",pch=18,cex=1.5)

arrows(px1,ly1,px1,uy1, length = 0.05, angle = 90, code = 3,col="orange")
arrows(lx1,py1,ux1,py1, length = 0.05, angle = 90, code = 3,col="orange")
arrows(px3,ly3,px3,uy3, length = 0.05, angle = 90, code = 3,col="deepskyblue")
arrows(lx3,py3,ux3,py3, length = 0.05, angle = 90, code = 3,col="deepskyblue")
arrows(px4,ly4,px4,uy4, length = 0.05, angle = 90, code = 3,col="green3")
arrows(lx4,py4,ux4,py4, length = 0.05, angle = 90, code = 3,col="green3")
arrows(px2,ly2,px2,uy2, length = 0.05, angle = 90, code = 3,col="hotpink")
arrows(lx2,py2,ux2,py2, length = 0.05, angle = 90, code = 3,col="hotpink")



par(pin = c(4,4))
plot(NA,NA,xlim=c(-2.3,1.3),ylim=c(-1,1),xlab="difference of tL (LD)",ylab="difference of tL (SD)",asp=1)
abline(h=0,col="black")
abline(v=0,col="black")
abline(h=c(2,-2,1,-1),col="gray",lty=2)
abline(v=c(2,-2,1,-1),col="gray",lty=2)
points(px1,py1,col="orange",pch=19,cex=1.5)
points(px3,py3,col="deepskyblue",pch=15,cex=1.5)
points(px4,py4,col="green3",pch=17,cex=1.5)
points(px2,py2,col="hotpink",pch=18,cex=1.5)

arrows(px1,ly1,px1,uy1, length = 0.05, angle = 90, code = 3,col="orange")
arrows(lx1,py1,ux1,py1, length = 0.05, angle = 90, code = 3,col="orange")
arrows(px3,ly3,px3,uy3, length = 0.05, angle = 90, code = 3,col="deepskyblue")
arrows(lx3,py3,ux3,py3, length = 0.05, angle = 90, code = 3,col="deepskyblue")
arrows(px4,ly4,px4,uy4, length = 0.05, angle = 90, code = 3,col="green3")
arrows(lx4,py4,ux4,py4, length = 0.05, angle = 90, code = 3,col="green3")
arrows(px2,ly2,px2,uy2, length = 0.05, angle = 90, code = 3,col="hotpink")
arrows(lx2,py2,ux2,py2, length = 0.05, angle = 90, code = 3,col="hotpink")
