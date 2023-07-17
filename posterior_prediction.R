#This code calculates the posterior distribution of the prediction from the parameter samples by MCMC

gen_bat_list = read.csv("genotype_batch_matching.csv")
gens = gen_bat_list$Genotype
bats = gen_bat_list$Batch

cur_dir = getwd()
new_dir = paste(cur_dir,"/result/posterior_prediction",sep="")
dir.create(new_dir)

for(j in 1:length(gens)){
  
  gen = gens[j]
  bat = bats[j]
  print(gen)
  
  #reading the parameter samples
  par_sample = read.csv(paste("result/parameter_sample/sample_",bat,"_",gen,".csv",sep=""))
  
  #stan input
  t0 = 0
  te = 5
  dt = 0.01
  Ni = as.integer((te-t0)/dt)
  fit0 = 3
  fite = 5
  thin = 5
  sample_N = 1000
  warmup = as.integer(fit0/dt)
  Nrec = as.integer((Ni-warmup)/thin)
  sample_index = sample(1:Ns,sample_N)
  
  data_size = 2*Nrec*4*sample_N
  print(data_size)
  
  
  simdata = list(
    iter = 1000,
    Nrec = Nrec,
    warmup = warmup,
    thin = thin,
    Ni = Ni,
    t0 = t0,
    dt = dt,
    a      = (par_sample$a)[sample_index],
    cdwn   = 1/(par_sample$bmax)[sample_index],
    h_LD   = (par_sample$h_LD)[sample_index],
    h_SD   = (par_sample$h_SD)[sample_index],
    tL_LD  = (par_sample$tL_LD)[sample_index],
    tL_SD  = (par_sample$tL_SD)[sample_index],
    r      = (par_sample$r)[sample_index]
  )
  
  #stan 
  
  simmodel = stan_model(file="posterior_prediction.stan")
  
  st_sim   = sampling(simmodel,
                      data  = simdata,
                      chains= 1,
                      iter  = 1,
                      algorithm   = "Fixed_param"
  )
  
  sim = extract(st_sim)
  
  Y_LD = (sim$Y_LD)[1,,]
  Y_SD = (sim$Y_SD)[1,,]
  
  t_LD = Y_LD[seq(1,dim(Y_LD)[1],4),]
  t_SD = Y_SD[seq(1,dim(Y_LD)[1],4),]
  c_LD = Y_LD[seq(2,dim(Y_LD)[1],4),]
  c_SD = Y_SD[seq(2,dim(Y_LD)[1],4),]
  s_LD = Y_LD[seq(3,dim(Y_LD)[1],4),]
  s_SD = Y_SD[seq(3,dim(Y_LD)[1],4),]
  b_LD = Y_LD[seq(4,dim(Y_LD)[1],4),]
  b_SD = Y_SD[seq(4,dim(Y_LD)[1],4),]
  
  #summarizing the stan output
  
  sim_result = data.frame(matrix(NA,nrow=4*sample_N,ncol=(Nrec+3)))
  tNrecs = paste(rep("t",Nrec),seq(1:Nrec),sep="")
  colnames(sim_result) = c("ID","Photoperiod","Label",tNrecs)
  
  for(i in 1:sample_N){
    sim_result[(4*(i-1)+1):(4*i),1] = rep(i,4)
    sim_result[(4*(i-1)+1):(4*i),2] = rep("LD",4)
    sim_result[(4*(i-1)+1),3:(Nrec+3)] = c("Time",t_LD[i,])
    sim_result[(4*(i-1)+2),3:(Nrec+3)] = c("Starch",c_LD[i,])
    sim_result[(4*(i-1)+3),3:(Nrec+3)] = c("Sucrose",s_LD[i,])
    sim_result[(4*(i-1)+4),3:(Nrec+3)] = c("Beta",b_LD[i,])
  }
  for(i in 1:sample_N){
    sim_result[(4*(i-1)+1):(4*i),1] = rep(i,4)
    sim_result[(4*(i-1)+1):(4*i),2] = rep("SD",4)
    sim_result[(4*(i-1)+1),3:(Nrec+3)] = c("Time",t_SD[i,])
    sim_result[(4*(i-1)+2),3:(Nrec+3)] = c("Starch",c_SD[i,])
    sim_result[(4*(i-1)+3),3:(Nrec+3)] = c("Sucrose",s_SD[i,])
    sim_result[(4*(i-1)+4),3:(Nrec+3)] = c("Beta",b_SD[i,])
  }
  write.csv(sim_result,paste("result/posterior_prediction/posterior_prediction_",bat,"_",gen,".csv",sep=""),row.names=F)
  
}



