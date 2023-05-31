# This code visualizes the posterior distribution of the model prediction

library(HDInterval)

# function #####################################################################
light_function = function(ts,dusk,Tday){
  lts = numeric(length=length(ts))
  for(i in 1:length(ts)){
    t = ts[i]%%Tday
    if (t <= dusk){
      lts[i] = 1
    } else {
      lts[i] = 0
    }
  }
  return(lts)
}

light_background = function(dusk,Tday,xlim,ylim){
  xmin = xlim[1]
  xmax = xlim[2]
  dt = (xmax-xmin)/1000
  ts = seq(xmin,xmax,dt)
  lts = light_function(ts,dusk,Tday)
  interval = length(ts)
  xleft = ts[1]-dt
  ymin=ylim[1]
  ymax=ylim[2]
  
  for (i in 1:(interval-1)){
    lt = lts[i]
    if(lt==1){
      col = "snow"
    } else {
      col="ivory2"
    }
    
    xleft = ts[i]
    xright = ts[i+1]
    rect(xleft,ymin-10,xright,ymax+10,col=col,border=col)
  }
}

# plot by one genotype #########################################################

# target genotype
gen = "wt"
bat = "1"

# target metabolite and photoperiod 
met = 3 # Starch:1, Sucrose:2, Beta:3
pho = 2 # LD:1, SD:2

# dataset
dataset = read.csv("dataset20220721.csv")
dataset["Beta"] = (dataset$Maltose)/((dataset$Starch)^(2/3))

met_list = c("Starch","Sucrose","Beta")
pho_list = c("LD","SD")
label_list = c("Concentration [umol/g]","Concentration [umol/g]","Relative Degradation rate")
xmin = 0
xmax = 1
ymins = c(0,0,0)
ymaxs = c(70,3,2)

# preparing a plot window
if(pho == 1){
  dusk = 16/24
} else {
  dusk = 8/24
}
par(mai=c(1,1,0.8,0.3))
par(xpd=FALSE)
plot(NA,NA,xlim=c(xmin,xmax),ylim=c(ymins[met],ymaxs[met]),xlab="Time",ylab=label_list[met],main=paste(met_list[met],pho_list[pho],sep=" "),xaxt="n")
axis(side=1,at=seq(xmin,xmax,by=1/6),label=seq(0,24,4))
light_background(dusk,1,c(xmin-1,xmax+1),c(ymins[met],ymaxs[met]))
col = "deepskyblue"
pch = 19

# Data 

subdata = subset(dataset,(dataset$Batch==bat)&(dataset$Genotype==gen))
data = subset(subdata,subdata$Photoperiod==pho_list[pho])
Ts = data$Time
X = data[met_list[met]][,1]

if(met==3){
  timeset = unique(T_wt)
  b0_mean = mean((subset(X_wt,T_wt==timeset[0])))
  X_wt = X_wt/b0_mean
}

# model prediction 

sim_result = read.csv(paste("result/posterior_prediction/posterior_prediction_",bat,"_",gen,".csv",sep=""))
Nt = length(sim_result[1,]) - 3
t = (subset(sim_result,sim_result$Label=="Time"))[1,4:(Nt+3)]
x = subset(sim_result,(sim_result$Label==met_list[met])&(sim_result$Photoperiod==pho_list[pho]))[4:(Nt+3)]

# mean and 95% interval (HDI) 

stats = data.frame(matrix(NA,ncol=4,nrow=Nt))
colnames(stats) = c("Time","Mean","2.5%","97.5%")
for(j in 1:Nt){
  stats[j,1] = t[j]
  stats[j,2] = mean(x[,j])
  stats[j,3] = (hdi(x[,j]))[1]
  stats[j,4] = (hdi(x[,j]))[2]
}

if(met == 3){
  coef = max(stats[,2])
  stats[,c(2,3,4)] = stats[,c(2,3,4)]/coef
}

# plot 
x = stats[,1]-1
y = stats[,2]
y_low = stats[,3]
y_upp = stats[,4]
rgb = col2rgb(col)
mycol <- rgb(rgb[1], rgb[2], rgb[3], max = 255, alpha = 125)

polygon(c(x,rev(x)),c(y_upp,rev(y_low)),col=mycol,border=NA)
lines(x,y,"l",col=col,lwd=2)
points(Ts,X,pch=pch,col=mycol,cex=1)


# plot all mutants #############################################################

gens = c("gwd","pwd","bam1","bam2","amy3-bam1","sweet1112","prr7","prr7prr9","ztl3","rb5","rb22","kin10","flz14-1","bzip63-1")

# plot options
col_wt = "deepskyblue"
col_mt = "deeppink"
pch_wt = 19
pch_mt = 19

# dataset
dataset = read.csv("dataset20220721.csv")
dataset["Beta"] = (dataset$Maltose)/((dataset$Starch)^(2/3))
met_list = c("Starch","Sucrose","Beta")
pho_list = c("LD","SD")
label_list = c("Concentration [umol/g]","Concentration [umol/g]","Relative Degradation rate")

xmin = 0
xmax = 1
ymins = c(0,0,0)
ymaxs = c(70,5,2)

for(gen in gens){
  print(gen)
  for(pho in 1:2){
    for(met in 1:3){
      
      # preparing a plot window
      if(pho == 1){
        dusk = 16/24
      } else {
        dusk = 8/24
      }
      
      par(pin = c(2.6,1.6))
      par(xpd=FALSE)
      plot(NA,NA,xlim=c(xmin,xmax),ylim=c(ymins[met],ymaxs[met]),xlab="Time",ylab=label_list[met],main=paste(met_list[met],pho_list[pho],sep=" "),xaxt="n",cex.axis=1.2)
      axis(side=1,at=seq(xmin,xmax,by=1/6),label=seq(0,24,4),cex.axis=1.2)
      light_background(dusk,1,c(xmin-1,xmax+1),c(ymins[met],ymaxs[met]))
      box(lwd=1.3) 
      
      # reading and processing the data 
      gen_bat_list = read.csv("mutant_batch_matching.csv")
      bat = (subset(gen_bat_list,gen_bat_list$Genotype==gen))$Batch
      
      # data of wild-type
      subdata_wt = subset(dataset,(dataset$Batch==bat)&(dataset$Genotype=="wt"))
      data_wt = subset(subdata_wt,subdata_wt$Photoperiod==pho_list[pho])
      T_wt = data_wt$Time
      X_wt = data_wt[met_list[met]][,1]
      
      if(met==3){
        timeset = unique(T_wt)
        b0_mean =mean((subset(X_wt,T_wt==timeset[0]))) 
        coef_b = b0_mean
        X_wt = X_wt/coef_b
      }
      
      # data of mutant
      subdata_mt = subset(dataset,(dataset$Batch==bat)&(dataset$Genotype==gen))
      data_mt = subset(subdata_mt,subdata_mt$Photoperiod==pho_list[pho])
      T_mt = data_mt$Time
      X_mt = data_mt[met_list[met]][,1]
      
      if(met==3){
        X_mt = X_mt/coef_b
      }
      
      # model prediction of wild-type
      sim_result_wt = read.csv(paste("result/simulation2/simulation_",bat,"_","wt",".csv",sep=""))
      Nt_wt = length(sim_result_wt[1,]) - 3
      t_wt = (subset(sim_result_wt,sim_result_wt$Label=="Time"))[1,4:(Nt_wt+3)]
      x_wt = subset(sim_result_wt,(sim_result_wt$Label==met_list[met])&(sim_result_wt$Photoperiod==pho_list[pho]))[4:(Nt_wt+3)]
      
      # model prediction of mutants
      sim_result_mt = read.csv(paste("result/simulation2/simulation_",bat,"_",gen,".csv",sep=""))
      Nt_mt = length(sim_result_mt[1,]) - 3
      t_mt = (subset(sim_result_mt,sim_result_mt$Label=="Time"))[1,4:(Nt_mt+3)]
      x_mt = subset(sim_result_mt,(sim_result_mt$Label==met_list[met])&(sim_result_mt$Photoperiod==pho_list[pho]))[4:(Nt_wt+3)]
      
      # mean and HDI of wild-type
      stats_wt = data.frame(matrix(NA,ncol=4,nrow=Nt_wt))
      colnames(stats_wt) = c("Time","Mean","2.5%","97.5%")
      for(j in 1:Nt_wt){
        stats_wt[j,1] = t_wt[j]
        stats_wt[j,2] = mean(x_wt[,j])
        stats_wt[j,3] = (hdi(x_wt[,j]))[1]
        stats_wt[j,4] = (hdi(x_wt[,j]))[2]
      }
      
      if(met == 3){
        coef = stats_wt[0,2]
        stats_wt[,c(2,3,4)] = stats_wt[,c(2,3,4)]/coef
      }
      
      # mean and HDI of mutant
      stats_mt = data.frame(matrix(NA,ncol=4,nrow=Nt_wt))
      colnames(stats_mt) = c("Time","Mean","2.5%","97.5%")
      for(j in 1:Nt_mt){
        stats_mt[j,1] = t_mt[j]
        stats_mt[j,2] = mean(x_mt[,j])
        stats_mt[j,3] = (hdi(x_mt[,j]))[1]
        stats_mt[j,4] = (hdi(x_mt[,j]))[2]
      }
      
      if(met == 3){
        stats_mt[,c(2,3,4)] = stats_mt[,c(2,3,4)]/coef
      }
      
      # plot wild-type
      x = stats_wt[,1]-1
      y = stats_wt[,2]
      y_low = stats_wt[,3]
      y_upp = stats_wt[,4]
      col = col_wt
      rgb = col2rgb(col)
      mycol <- rgb(rgb[1], rgb[2], rgb[3], max = 255, alpha = 125)
      
      polygon(c(x,rev(x)),c(y_upp,rev(y_low)),col=mycol,border=NA)
      lines(x,y,"l",col=col,lwd=2)
      points(T_wt,X_wt,pch=pch_wt,col=mycol,cex=1)
      
      # plot mutant
      x = stats_mt[,1]-1
      y = stats_mt[,2]
      y_low = stats_mt[,3]
      y_upp = stats_mt[,4]
      col = col_mt
      rgb = col2rgb(col)
      mycol <- rgb(rgb[1], rgb[2], rgb[3], max = 255, alpha = 125)
      
      polygon(c(x,rev(x)),c(y_upp,rev(y_low)),col=mycol,border=NA)
      lines(x,y,"l",col=col,lwd=2)
      points(T_mt,X_mt,pch=pch_mt,col=mycol,cex=1)
      dev.off()
    }
  }
  
}


