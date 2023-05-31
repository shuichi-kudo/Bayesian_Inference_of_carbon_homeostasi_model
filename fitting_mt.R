# This program execute stan fitting of all mutants using estimated parameter distributions of WT 
# Return result of model selection and MCMC samples of each model and mutants
# Developed :2022/03/10 
# Final update : 2023/01/31 
# Before executing this mutant analysis, you need to complete the wild type analysis using "fitting_wt.R"
# You need enough memory size 20~30GB
# Settings ######################################################################

model_list = read.csv("model_list.csv")
gens = c("pwd","gwd","bam1","bam2","amy3-bam1",
         "rb5","rb22","kin10","bzip63-1","sweet1112",
         "prr7","prr7prr9","ztl3","flz14-1")

# Sampling setting

warmup_day = 3
fit_day = 3
dt = 0.01
iteration = 2000
chains = 4
fitseed = 10000
estimation_method = "mode"
result_file = "result/model_comp_result.csv"
sample_file_name = "result/stan_mt/sampling_bat_gen_model.csv"
dataset_file = "dataset.csv"

t0 = 0
te = warmup_day + fit_day

cur_dir = getwd()
new_dir = paste(cur_dir,"/stan_mt/stan_result",sep="")
dir.create(new_dir)


# Functions ####################################################################
data_replication=function(dataset,fit0,fite){
  rpt = fite - fit0
  dataset0 = subset(dataset,dataset$Time < 1)
  repdata = data.frame(dataset0)
  for(i in 1:(rpt)){
    adddata = data.frame(dataset0)
    adddata$Time = dataset0$Time + i
    addeddata = rbind(repdata,adddata)
    repdata = addeddata
  }
  rtrn = data.frame(repdata)
  rtrn$Time = repdata$Time + fit0 - 1
  return(rtrn)
}

time_index = function(dataset,t0,te,dt){
  Times = dataset$Time
  Nd = length(Times)
  Ni = as.integer((te-t0)/dt)
  I = numeric(length=Nd)
  for(i in 1:Nd){
    tdata = Times[i]
    for(j in 1:Ni){
      tsim = t0 + (j-1)*dt
      diff = tdata - tsim
      if(abs(tdata-tsim) <= (dt/2)){
        I[i] = j
        break
      }
    }
  }
  return (I)
}

data_processing = function(data,t0,te,dt,fit0,fite){
  Ns = dim(data)[1]
  cols = c("Time","Photoperiod","Starch","Sucrose")
  Data_ret = data.frame(matrix(NA,nrow=Ns,ncol=length(cols)))
  colnames(Data_ret) = cols
  Timeindex = time_index(data,t0,te,dt)
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

my_waic = function(st_fit){
  log_lik = (extract(st_fit))$log_likelihood
  lppd = sum(log10(colMeans(exp(log_lik))))
  p_waic = sum(apply(log_lik,2,var))
  waic = -2*lppd + 2*p_waic
  return(waic)
}

# main #######################################################################

# Batch of the specific Genotype
gen_bat_list = read.csv("mutant_batch_matching.csv")

# Model list
model_names = model_list$Model
Nps = model_list$Np

#output files
result = data.frame(matrix(data=NA,nrow=(length(gens)*length(model_names)),ncol=17))
colnames(result) = c("Genotype","Model","a","cdwn","h_LD","h_SD","tL_LD","tL_SD","r",
                     "Np","ESS","Rhat","Lpd","AIC","WAIC","rank_AIC","rank_WAIC")
#read package
library(rstan)
rstan_options(auto_write=TRUE)
options(mc.cores = parallel::detectCores())

#read dataset
Data0 = read.csv(dataset_file)

fit0 = warmup_day
fite = warmup_day + fit_day


for(i in 1:length(gens)){
  # target mutant
  gen = gens[i]
  bat = batchs[i]
  print(gen)
  
  # fix parameters
  post_sum_wt = read.csv(paste("result/posterior_summary_",bat,"_wt.csv",sep=""))
  par_wt = subset(post_sum_wt,post_sum_wt$method==estimation_method)
  a = par_wt$a
  cdwn = par_wt$cdwn
  h_LD = par_wt$h_LD
  h_SD = par_wt$h_SD
  tL_LD = par_wt$tL_LD
  tL_SD = par_wt$tL_SD
  r = par_wt$r
  
  # data processing
  Data1 = subset(Data0,Data0$Genotype==gen)
  Data_rep = data_replication(Data1,fit0,fite)
  Data = data_processing(Data_rep,t0,te,dt,fit0,fite)
  Data_LD = subset(Data,Data$Photoperiod==1)
  Data_SD = subset(Data,Data$Photoperiod==2)
  I_LD = Data_LD[,1:2]
  I_SD = Data_SD[,1:2]
  X_LD = Data_LD[,3:4]
  X_SD = Data_SD[,3:4]
  Ns_LD = dim(I_LD)[1]
  Ns_SD = dim(I_SD)[1]
  Ni = as.integer((te-t0)/dt)
  
  # initialization
  set.seed(100)
  initial_pars = list(
    a      = runif(1, 100, 200),
    cdwn   = runif(1,   0,  20),
    h_LD   = runif(1,  50, 100),
    h_SD   = runif(1,  50, 100),
    tL_LD  = runif(1, 0.5, 0.9),
    tL_SD  = runif(1, 0.3, 0.6),
    r      = runif(1, 0.6, 0.9)
  )
  
  # stan input
  fitdata = list(
    Ni = Ni,
    Ns_LD = Ns_LD,
    Ns_SD = Ns_SD,
    t0 = t0,
    te = te,
    dt = dt,
    a = a,
    cdwn = cdwn,
    h_LD = h_LD,
    h_SD = h_SD,
    tL_LD = tL_LD,
    tL_SD = tL_SD,
    r = r,
    I_LD = I_LD,
    I_SD = I_SD,
    X_LD = X_LD,
    X_SD = X_SD,
    C0_LD = C0_LD,
    C0_SD = C0_SD
  )
  
  # model selection
  Nm = length(model_names)
  for(j in 1:Nm){
    # stan compile
    model_name = model_names[j]
    Np = Nps[j]
    model_j = subset(model_list,model_list$Model==model_name)
    fit_file = paste("mutant_models_stan/","fitting_",model_name,".stan",sep="")
    sample_file = paste("result/stan_mt/","sampling_",bat,"_",gen,"_",model_name,".csv",sep="")
    stanmodel = stan_model(file=fit_file)
    st_fit = sampling(stanmodel,
                      data   = fitdata,
                      init   = function()initial_pars,
                      chains = chains,
                      iter   = iteration,
                      seed   = fitseed,
                      sample_file=sample_file)
    waic = my_waic(st_fit)
    st_sm = (summary(st_fit))$summary
    lp = st_sm[dim(st_sm)[1],1]
    aic = -2*lp + 2*Np
    rhats = st_sm[1:Np,10]
    esss = st_sm[1:Np,9]
    k = j + Nm*(i-1)
    result[k,1] = gen
    result[k,2] = model_name
    result[k,3] = model_j$a
    result[k,4] = model_j$cdwn
    result[k,5] = model_j$h_LD
    result[k,6] = model_j$h_SD
    result[k,7] = model_j$tL_LD
    result[k,8] = model_j$tL_SD
    result[k,9] = model_j$r
    result[k,10] = model_j$Np
    result[k,11] = min(esss)
    result[k,12] = max(rhats)
    result[k,13] = lp
    result[k,14] = aic
    result[k,15] = waic
  }
  #rank
  k = (i-1)*Nm
  result[(k+1):(k+Nm),16] = rank(result[(k+1):(k+Nm),13])
  result[(k+1):(k+Nm),17] = rank(result[(k+1):(k+Nm),14])
  write.csv(result,paste(result_file,sep=""),row.names=F)
}
# output
write.csv(result,paste(result_file,sep=""),row.names=F)
