#plot the posterior prediction of the model with estimated parameter
#based on posterior_prediction _b_wt.csv

#read library
library(HDInterval)
library(RColorBrewer)

#functions
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

#Settings
gen = "wt"
batches = c(1,2,3,4,5,6,7,9,12,13,14)
mets = c("Starch","Sucrose","Beta")
phos = c("LD","SD")
xmin = 0
xmax = 1
ymins = c(0,0,0)
ymaxs = c(70,3,2)
label_list = c("Concentration [umol/g]","Concentration [umol/g]","Relative Degradation rate")
#Data
dataset = read.csv("dataset.csv")
#Beta
dataset["Beta"] = dataset["Maltose"]/((dataset["Starch"])^(2/3))


#Main
cols = brewer.pal(length(batches),"Spectral")
for(i in 1:length(mets)){
  for(j in 1:length(phos)){
    met = mets[i]
    pho = phos[j]
    
    #prepare a plot window
    if(pho == "LD"){dusk = 16/24} else {dusk = 8/24}
    if(met == "Beta"){
      par(pin = c(3,2.5),ps=15) #Beta
    } else {
      par(pin = c(3,1.8)) # Starch and sucrose
    }
    plot(NA,NA,xlim=c(xmin,xmax),ylim=c(ymins[i],ymaxs[i]),xlab="Time",ylab=label_list[i],main=paste(met,pho,sep=" "),xaxt="n",cex.axis=1.1)
    axis(side=1,at=seq(xmin,xmax,by=1/6),label=seq(0,24,4),cex.axis=1.1)
    light_background(dusk,1,c(xmin-1,xmax+1),c(ymins[i],ymaxs[i]))
    box(lwd=1.3) 
    for(k in 1:length(batches)){
      batch = batches[k] 
      #read data
      subdata = subset(dataset,(dataset$Batch==batch)&(dataset$Genotype==gen))
      data = subset(subdata,subdata$Photoperiod==pho)
      Ts = data$Time
      X = data[met][,1]
      if(met=="Beta"){
        #mean at each time point
        timeset = unique(Ts)
        nt = length(timeset)
        b_mean = numeric(nt)
        for(t in 1:nt){
          b_mean[t] = mean((subset(X,Ts==timeset[t])))
        }
        coef_b = max(b_mean)
        X = X/coef_b
      }
      #simulation result
      sim_result = read.csv(paste("result/simulation/simulation_",batch,"_",gen,".csv",sep=""))
      Nt = length(sim_result[1,]) - 3
      t = (subset(sim_result,sim_result$Label=="Time"))[1,4:(Nt+3)]
      x = subset(sim_result,(sim_result$Label==met)&(sim_result$Photoperiod==pho))[4:(Nt+3)]
      #statistics (HDI,mean)
      stats = data.frame(matrix(NA,ncol=4,nrow=Nt))
      colnames(stats) = c("Time","Mean","2.5%","97.5%")
      for(j in 1:Nt){
        stats[j,1] = t[j]
        stats[j,2] = mean(x[,j])
        stats[j,3] = (hdi(x[,j]))[1]
        stats[j,4] = (hdi(x[,j]))[2]
      }
      if(met == "Beta"){
        coef = max(stats[,2])
        stats[,c(2,3,4)] = stats[,c(2,3,4)]/coef
      }
      x = stats[,1]-1
      y = stats[,2]
      y_low = stats[,3]
      y_upp = stats[,4]
      col = cols[k]
      rgb = col2rgb(col)
      mycol <- rgb(rgb[1], rgb[2], rgb[3], max = 255, alpha = 100)
      polygon(c(x,rev(x)),c(y_upp,rev(y_low)),col=mycol,border=NA)
      lines(x,y,"l",col=col,lwd=2)
      points(Ts,X,pch=19,col=mycol,cex=0.8)
    }
    #legend("topleft",legend=batches,col=cols,lwd=2)
  }
}



