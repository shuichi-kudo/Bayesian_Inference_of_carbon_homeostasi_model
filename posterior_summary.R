#This code calculates a summary table of the posterior estimates from the parameter samples

#Target
gen_bat_list = read.csv("genotype_batch_matching.csv")
gens = gen_bat_list$Genotype
bats = gen_bat_list$Batch
pars = c("a","r","h_LD","h_SD","bmax","tL_LD","tL_SD")

#Function: mode for MAP estimates
my_mode = function(samples){
  if(samples[1]==samples[2]){
    mode = samples[1]
  } else {
    dens = density(samples)
    x = dens$x
    y = dens$y
    mode = (x[which(y==max(y))])[1]
  }
  return(mode)
}

cols = c("Genotype","Batch",pars)
ps = data.frame(matrix(NA,ncol=length(cols),nrow=length(gens)))

for(i in 1:length(gens)){
  gen = gens[i]
  bat = bats[i]
  samples = read.csv(paste("result/parameter_sample/sample_",bat,"_",gen,".csv",sep=""))
  ps[i,1] = gen
  ps[i,2] = bat
  ps[i,3] = my_mode(samples$a)
  ps[i,4] = my_mode(samples$r)
  ps[i,5] = my_mode(samples$h_LD)
  ps[i,6] = my_mode(samples$h_SD)
  ps[i,6] = my_mode(samples$bmax)
  ps[i,7] = my_mode(samples$tL_LD)
  ps[i,8] = my_mode(samples$tL_SD)
}

colnames(ps) = cols

write.csv(ps,"posterior_summary.csv")
