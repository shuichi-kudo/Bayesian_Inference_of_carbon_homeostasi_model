#Return MCMC samples from stan sampling file
#Output:　result/parameter_sample/sample_bat_gen.csv

#Target
gen_bat_list = read.csv("genotype_batch_matching.csv")
gens = gen_bat_list$Genotype
bats = gen_bat_list$Batch

#Rename parameters
pars = c("a","cdwn","h_LD","h_SD","tL_LD","tL_SD","r")
new_pars = c("a","bmax","h_LD","h_SD","tL_LD","tL_SD","r")
rank = 1

#Read the result of the model selection
model_comp_result = read.csv("result/model_comp_result.csv")

#Read library
library("rstan")

for(i in 1:length(gens)){
  gen = gens[i]
  bat = bats[i]
  if(gen == "wt"){
    # wt
    sample_name = paste("result/wt/sampling_",bat,"_",gen,sep="")
    sample_files = sprintf(paste(sample_name,"%s.csv",sep="_"),1:4)
    sample = read_stan_csv(csvfiles=sample_files)

    trace = (stan_trace(sample,pars=pars))$data
    Np = length(pars)
    Ns = as.integer(length(trace$parameter)/Np)
    par_sample = data.frame(matrix(0,nrow=Ns,ncol=Np))
    colnames(par_sample)=new_pars
    for(j in 1:length(pars)){
      par = pars[j]
      par_data = subset(trace,trace$parameter==par)
      if(par == "cdwn"){
        for(k in 1:Ns){
          par_sample[k,j] = 1/(par_data$value)[k]
        }
      } else {
        for(k in 1:Ns){
          par_sample[k,j] = (par_data$value)[k]
        }
      }
    }
    
  } else {
    # mutant
    # Selecting the best model based on WAIC and R_hat (test of convergence)
    flag = F
    rank_candidate = rank
    while(flag=F){
      best_model_candidate = subset(model_comp_result,(model_comp_result$Genotype==gen)&(model_comp_result$rank_WAIC==rank_candidate))
      R_hat = best_model_candidate$Rhat
      if (R_hat < 1.1){
        flag = T
      } else {
        flag = F
        rank_candidate = rank_candidate + 1
      }
    }
    best_model = best_model_candidate
    model_name = best_model$Model
    inc_par_vec = best_model[1,3:8]
    sample_name = paste("result/stan_result/sampling_",bat,"_",gen,"_",model_name,sep="")
    sample_files = sprintf(paste(sample_name,"%s.csv",sep="_"),1:4)
    sample = read_stan_csv(csvfiles=sample_files)
    post_sum_wt = read.csv(paste("result/posterior_summary_",bat,"_wt.csv",sep=""))
    
    inc_pars = c()
    for(j in 1:length(inc_par_vec)){
      if(inc_par_vec[j]==1){
        inc_pars = c(inc_pars,pars[j])
      }
    }
    trace = (stan_trace(sample,pars=inc_pars))$data
    Np = length(pars)
    Ns = as.integer(length(trace$parameter)/length(inc_pars))
    par_sample = data.frame(matrix(0,nrow=Ns,ncol=Np))
    colnames(par_sample)=new_pars
    
    for(j in 1:length(pars)){
      par = pars[j]
      if(inc_par_vec[j]==1){ # = = = = = = = = = = = = = = = = = = Delta != 0
        par_data = subset(trace,trace$parameter==par)
        if(par=="cdwn"){
          for(k in 1:Ns){
            par_sample[k,j] = 1/(par_data$value)[k]
          }
        } else {
          for(k in 1:Ns){
            par_sample[k,j] = (par_data$value)[k]
          }
        }
      } else { # = = = = = = = = = = = = = = = = = = = = = = = = = Delta == 0
        if(par=="cdwn"){
          for(k in 1:Ns){
            par_sample[k,j] = 1/post_sum_wt[3,1+j] #median of the wild-type sample
          }
        }
        for(k in 1:Ns){
          par_sample[k,j] = post_sum_wt[3,1+j]
        }
      }
    }
  }
  
  write.csv(par_sample,paste("result/parameter_sample/sample_",bat,"_",gen,".csv",sep=""),row.names=F)
  
}
