functions{
  real beta_L(real t,real a, real cdwn, real tL, real r, real kappa){
    real tD;
    real beta;
    tD = 1 - tL;
    if(tD > r){
      beta = 0;
    } else {
      beta = a*(r-tD)/((a*tD*t+cdwn)^kappa);
    }
    return beta;
  }
  real beta_D(real t,real a, real cdwn, real tL, real r, real kappa){
    real tD;
    real beta;
    tD = 1 - tL;
    if(tD > r){
      beta = a*r*(tD^(kappa-1))*(1-tD)/((a*r*(1-tD)*(1-t)+cdwn*tD)^kappa);
    } else {
      beta = a*tL/((a*tL*(1-t)+cdwn)^kappa);
    }
    return beta;
  }
  real dual_beta_L(real t,real a, real cdwn, real tL, real r, real kappa){
    real tD;
    real beta;
    real tday;
    real trev;
    tD = 1 - tL;
    tday = fmod(t,1);
    if(tday < tL){
      beta = beta_L(tday,a,cdwn,tL,r,kappa);
    } else {
      trev = (1-tday)*tL/tD;
      beta = beta_L(trev,a,cdwn,tL,r,kappa);
    }
    return beta;
  }
  real dual_beta_D(real t,real a, real cdwn, real tL, real r, real kappa){
    real tD;
    real beta;
    real tday;
    real trev;
    tD = 1- tL;
    tday = fmod(t,1);
    if(tday < tL){
      trev = 1 - (tD/tL)*tday;
      beta = beta_D(trev,a,cdwn,tL,r,kappa);
    } else {
      beta = beta_D(tday,a,cdwn,tL,r,kappa);
    }
    return beta;
  }
  real[] ODEs(real t,real c,real s,real a,real cdwn,real dusk,real h,real lambda,real tL,real r){
    real dydt[2];
    real Ct;
    real St;
    real dcdt;
    real dsdt;
    real beta;
    int lt;
    real tday;
    real tD;
    real kappa;
    //Time and Phase in a day
    tday = fmod(t,1);
    //External condition
    if(tday <= dusk){
      lt = 1;
    } else {
      lt = 0;
    }
    kappa = 2.0/3.0;
    tD = 1-tL;
    //Starch degradation rate
    beta = lt*dual_beta_L(t,a,cdwn,tL,r,kappa)+(1-lt)*dual_beta_D(t,a,cdwn,tL,r,kappa);
    //Starch-Sucrose dynamics
    if(c<0){
      Ct = 0;
    } else {
      Ct = c;
    }
    if(s<0){
      St = 0;
    } else {
      St = s;
    }
    dcdt = a*r*lt - lambda*beta*Ct^kappa;
    dsdt = a*(1-r)*lt + lambda*beta*Ct^kappa - h*St;
    //return
    dydt[1] = dcdt;
    dydt[2] = dsdt;
    return dydt;
  }
  real[,] rungekutta(int N_warmup,int N_fit, real dt,real a,real cdwn,real dusk,real h,real lambda,real tL,real r){
    real y[3,N_fit];
    real k1[2];
    real k2[2];
    real k3[2];
    real k4[2];
    real k5[2];
    real t;
    real c;
    real s;
    //initial set
    t = 0;
    c = cdwn;
    s = a*tL/h;
    //warm up
    for(i in 1:N_warmup){
      //gradients
      k1 = ODEs(t       , c             , s             , a, cdwn, dusk, h, lambda, tL, r);
      k2 = ODEs(t+0.5*dt, c+0.5*dt*k1[1], s+0.5*dt*k1[2], a, cdwn, dusk, h, lambda, tL, r);
      k3 = ODEs(t+0.5*dt, c+0.5*dt*k2[1], s+0.5*dt*k2[2], a, cdwn, dusk, h, lambda, tL, r);
      k4 = ODEs(t+    dt, c+    dt*k3[1], s+    dt*k3[2], a, cdwn, dusk, h, lambda, tL, r);
      //upgrade
      t = t + dt;
      c = c + dt*(k1[1] + 2*k2[1] + 2*k3[1] + k4[1])/6;
      s = s + dt*(k1[2] + 2*k2[2] + 2*k3[2] + k4[2])/6;
    }
    //calculation
    for(i in 1:N_fit){
      y[1,i] = t - dt*N_warmup;
      y[2,i] = c;
      y[3,i] = s;
      //gradients
      k1 = ODEs(t       , c             , s             , a, cdwn, dusk, h, lambda, tL, r);
      k2 = ODEs(t+0.5*dt, c+0.5*dt*k1[1], s+0.5*dt*k1[2], a, cdwn, dusk, h, lambda, tL, r);
      k3 = ODEs(t+0.5*dt, c+0.5*dt*k2[1], s+0.5*dt*k2[2], a, cdwn, dusk, h, lambda, tL, r);
      k4 = ODEs(t+    dt, c+    dt*k3[1], s+    dt*k3[2], a, cdwn, dusk, h, lambda, tL, r);
      //upgrade
      t = t + dt;
      c = c + dt*(k1[1] + 2*k2[1] + 2*k3[1] + k4[1])/6;
      s = s + dt*(k1[2] + 2*k2[2] + 2*k3[2] + k4[2])/6;
    }
    return y;
  }
}
data {
  int N_warmup;
  int N_fit;
  int Ns_LD;  //number of data
  int Ns_SD;
  real dt;
  //real a;
  //real cdwn;
  real h_LD;
  //real h_SD;
  //real tL_LD;
  real tL_SD;
  //real r;
  int I_LD[Ns_LD,2];
  int I_SD[Ns_SD,2];
  real X_LD[Ns_LD,2];
  real X_SD[Ns_SD,2];
}
transformed data {
  //fixed paraameters
  real dusk_LD;
  real dusk_SD;
  real lambda;
  dusk_LD = 16.0/24.0;
  dusk_SD =  8.0/24.0;
  lambda = 1.0;
}
parameters {
  real<lower=0> a;
  real<lower=0> cdwn;
  //real<lower=0> h_LD;
  real<lower=0> h_SD;
  real<lower=0,upper=1> tL_LD;
  //real<lower=0,upper=1> tL_SD;
  real<lower=0,upper=1> r;
  real<lower=0> sC;
  real<lower=0> sS;
}
model {
  real y_LD[3,N_fit];
  real y_SD[3,N_fit];
  int i_LD;
  int i_SD;
  y_LD = rungekutta(N_warmup, N_fit, dt, a, cdwn, dusk_LD, h_LD, lambda, tL_LD, r);
  y_SD = rungekutta(N_warmup, N_fit, dt, a, cdwn, dusk_SD, h_SD, lambda, tL_SD, r);
  for (i in 1:Ns_LD){
    i_LD = I_LD[i,1];
    X_LD[i,1] ~ normal(y_LD[2,i_LD],sC);
    X_LD[i,2] ~ normal(y_LD[3,i_LD],sS);
  }
  for (i in 1:Ns_SD){
    i_SD = I_SD[i,1];
    X_SD[i,1] ~ normal(y_SD[2,i_SD],sC);
    X_SD[i,2] ~ normal(y_SD[3,i_SD],sS);
  }
  //prior probability for batch1
  //a ~ normal(,);
  //cdwn ~ normal(,);
  //h ~ normal(,);
  //lambda ~ normal(1,3);
  //tL_LD ~ normal(,);
  //tL_SD ~ normal(,);
  //r ~ normal(,);
  //target += lognormal_lpdf(r-|,);
}
generated quantities{
  vector[(Ns_LD+Ns_SD)] log_likelihood;
  {
    //I don't want to save following intermediate variables
    real y_LD[3,N_fit];
    real y_SD[3,N_fit];
    int i_LD;
    int i_SD;
    y_LD = rungekutta(N_warmup, N_fit, dt, a, cdwn, dusk_LD, h_LD, lambda, tL_LD, r);
    y_SD = rungekutta(N_warmup, N_fit, dt, a, cdwn, dusk_SD, h_SD, lambda, tL_SD, r);
    for (i in 1:Ns_LD){
      i_LD = I_LD[i,1];
      log_likelihood[i] = normal_lpdf(X_LD[i,1]|y_LD[2,i_LD],sC)+normal_lpdf(X_LD[i,2]|y_LD[3,i_LD],sS);
    }
    for (i in 1:Ns_SD){
      i_SD = I_SD[i,1];
      log_likelihood[i+Ns_LD] = normal_lpdf(X_SD[i,1]|y_SD[2,i_SD],sC)+normal_lpdf(X_SD[i,2]|y_SD[3,i_SD],sS);
    }
  }
}
