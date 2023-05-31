#This file execute a fitting by stan and return sampling file of parameters

# Settings ######################################################################

#target
batches = c(1,2,3,4,5,6,7,9,12,13)

#detail settings
warmup_day = 10
fit_day = 3
dt = 0.01
iteration = 2000
chains = 4
fitseed = 42

dataset_file = "dataset.csv"
fit_file = "fitting_wt.stan"
parameters = c("a","cdwn","h_LD","h_SD","tL_LD","tL_SD","r")

#read rstan package
library(rstan)
rstan_options(auto_write=TRUE)
options(mc.cores = parallel::detectCores())

#make a directory to save the result
cur_dir = getwd()
new_dir1 = paste(cur_dir,"/result",sep="")
new_dir2 = paste(cur_dir,"/result/stan_wt",sep="")
dir.create(new_dir1)
dir.create(new_dir2)

# functions ####################################################################
data_replication=function(dataset,fit_day){
  rpt = fit_day-1
  dataset0 = subset(dataset,dataset$Time < 1)
  repdata = data.frame(dataset0)
  for(i in 1:rpt){
    adddata = data.frame(dataset0)
    adddata$Time = dataset0$Time + i
    addeddata = rbind(repdata,adddata)
    repdata = addeddata
  }
  rtrn = data.frame(repdata)
  return(rtrn)
}

time_index = function(t_data,fit_day,dt){
  Times = t_data
  Nd = length(Times)
  Ni = as.integer(fit_day/dt)
  I = numeric(length=Nd)
  for(i in 1:Nd){
    tdata = Times[i]
    diffs = numeric(Ni)
    for(j in 1:Ni){
      tsim = 0 + (j-1)*dt
      diff = abs(tdata - tsim)
      diffs[j] = diff
    }
    I[i] = which(diffs==min(diffs))
  }
  return (I)
}

data_processing = function(data,fit_day,dt){
  Ns = dim(data)[1]
  cols = c("Time","Photoperiod","Starch","Sucrose")
  Data_ret = data.frame(matrix(NA,nrow=Ns,ncol=length(cols)))
  colnames(Data_ret) = cols
  Timeindex = time_index(data$Time,fit_day,dt)
  Pho = data$Photoperiod
  Ct =data$Starch
  St = data$Sucrose
  for(i in 1:Ns){
    if(Pho[i]=="LD"){
      pho = 1
    } else {
      pho = 2
    }
    Data_ret[i,1] = Timeindex[i]
    Data_ret[i,2] = pho
    Data_ret[i,3] = Ct[i]
    Data_ret[i,4] = St[i]
  }
  return(Data_ret)
}

# sampling ####################################################################
gen = "wt"

#read dataset
Data0 = read.csv(dataset_file)
Nb = length(batches)

for(i in 1:Nb){
  bat = batches[i]
  
  #data processing
  Data1 = subset(Data0,(Data0$Batch==bat)&(Data0$Genotype==gen))
  Data_rep = data_replication(Data1,fit_day)
  Data = data_processing(Data_rep,fit_day,dt)
  Data_LD = subset(Data,Data$Photoperiod==1)
  Data_SD = subset(Data,Data$Photoperiod==2)
  I_LD = Data_LD[,1:2]
  I_SD = Data_SD[,1:2]
  X_LD = Data_LD[,3:4]
  X_SD = Data_SD[,3:4]
  Ns_LD = dim(I_LD)[1]
  Ns_SD = dim(I_SD)[1]
  N_warmup = as.integer(warmup_day/dt)
  N_fit = as.integer(fit_day/dt)
  
  #initialization
  set.seed(42)
  initial_pars = list(
    a      = runif(1, 100, 200),
    cdwn   = runif(1,   0,  20),
    h_LD   = runif(1,  50, 100),
    h_SD   = runif(1,  50, 100),
    tL_LD  = runif(1, 0.5, 0.9),
    tL_SD  = runif(1, 0.3, 0.6),
    r      = runif(1, 0.6, 0.9)
  )
  
  #stan input
  fitdata = list(
    N_warmup = N_warmup,
    N_fit = N_fit,
    Ns_LD = Ns_LD,
    Ns_SD = Ns_SD,
    dt = dt,
    I_LD = I_LD,
    I_SD = I_SD,
    X_LD = X_LD,
    X_SD = X_SD
  )
  
  #stan compile
  stanmodel = stan_model(file=fit_file)
  
  #stan sampling
  st_fit = sampling(stanmodel,
                    data   = fitdata,
                    init   = function()initial_pars,
                    chains = chains,
                    iter   = iteration,
                    seed   = fitseed,
                    sample_file=paste("result/stan_wt/sampling_",bat,"_","wt.csv",sep=""))
  
  trace = (stan_trace(st_fit,pars=parameters))$data
  Np = length(parameters)
  Ns = as.integer(length(trace$parameter)/Np)
  par_set = data.frame(matrix(NA,nrow=Ns,ncol=Np))
  colnames(par_set)=parameters
  for(i in 1:Np){
    par = parameters[i]
    par_data = subset(trace,trace$parameter==par)
    for(j in 1:Ns){
      par_set[j,i] = (par_data$value)[j]
    }
  }
  
  #Summary
  posterior_summary = data.frame(matrix(NA,nrow=3,ncol=(Np+1)))
  colnames(posterior_summary) = c("method",parameters)
  for(i in 1:Np){
    par_i = par_set[,i]
    mu = mean(par_i)
    md = median(par_i)
    dens = density(par_i)
    densx = dens$x
    densy = dens$y
    mo = densx[which.max(densy)]
    posterior_summary[1,(i+1)]=mu
    posterior_summary[2,(i+1)]=md
    posterior_summary[3,(i+1)]=mo
  }
  posterior_summary[1,1] = "mean"
  posterior_summary[2,1] = "median"
  posterior_summary[3,1] = "mode"
  
  print(posterior_summary)
  #output
  write.csv(posterior_summary,paste("result/posterior_summary_",bat,"_",gen,".csv",sep=""),row.names=F)
}



#convergence check #############################################################
for(i in 1:Nb){
  batch = batches[i]
  sampling_file_name = paste("result/stan_wt/sampling_",batch,"_",gen,sep="")
  chains = c(1,2,3,4)
  csvfiles = sprintf(paste(sampling_file_name,"%s.csv",sep="_"),chains)
  st_fit = read_stan_csv(csvfiles=csvfiles)
  sm = (summary(st_fit))$summary
  sm_par = sm[1:8,10]
  max_Rhat = max(sm_par)
  print(paste("Batch ",batch," : maximam Rhat is ",max_Rhat))
  if(max_Rhat > 1.1){
    print("maximam Rhat is too large!")
  }
}






