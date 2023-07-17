functions{
  
  real[] ODEs(real t,real c,real s,real a,real cdwn,real dusk,real h,real lambda,real tL,real r){
    real dydt[3];
    real dcdt;
    real dsdt;
    real betaL;
    real betaD;
    real beta;
    real lt;
    real lf;
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
    
    //Inner condition
    if(tday <= tL){
      lf = 1;
    } else {
      lf = 0;
    }
    
    kappa = 2.0/3.0;
    tD = 1-tL;
    //Starch degradation rate
    if(tD > r){
      beta = (1-lt)*(a*r*(tD^(kappa-1))*(1-tD))/(cdwn*tD+a*r*(1-tD)*(1-tday))^kappa;
    } else {
      beta = a*(lt*(r-tD)+(1-lt)*tL)*((lf/(a*tD*tday+cdwn)^kappa)+((1-lf)/(a*tL*(1-tday)+cdwn)^kappa));
    }
    
    //Starch-Sucrose dynamics
    dcdt = a*r*lt - lambda*beta*(fabs(c)^kappa);
    dsdt = a*(1-r)*lt + lambda*beta*(fabs(c)^kappa) - h*s;
    
    //return
    dydt[1] = dcdt;
    dydt[2] = dsdt;
    dydt[3] = beta;
    
    return dydt;
  }
  
  real[,] rungekutta(real t0, real dt,int Ni,real c0,real a,real cdwn,real dusk,real h,real lambda,real tL,real r){
    real y[4,Ni];
    real k1[3];
    real k2[3];
    real k3[3];
    real k4[3];
    real k5[3];
    real t;
    real c;
    real s;
    //initial set
    t = t0;
    c = c0;
    s = a*tL/h;
    
    //calculation
    for(i in 1:Ni){
      //gradients
      k1 = ODEs(t, c, s, a, cdwn, dusk, h, lambda, tL, r);
      k2 = ODEs(t+0.5*dt, c+0.5*dt*k1[1], s+0.5*dt*k1[2], a, cdwn, dusk, h, lambda, tL, r);
      k3 = ODEs(t+0.5*dt, c+0.5*dt*k2[1], s+0.5*dt*k2[2], a, cdwn, dusk, h, lambda, tL, r);
      k4 = ODEs(t+dt,     c+dt*k3[1],     s+dt*k3[2],     a, cdwn, dusk, h, lambda, tL, r);
      
      //record
      y[1,i] = t;
      y[2,i] = c;
      y[3,i] = s;
      y[4,i] = (k1[3] + 2*k2[3] + 2*k3[3] + k4[3])/6; //beta
      
      //upgrade
      t = t + dt;
      c = c + dt*(k1[1] + 2*k2[1] + 2*k3[1] + k4[1])/6;
      s = s + dt*(k1[2] + 2*k2[2] + 2*k3[2] + k4[2])/6;
      
    }
    
    return y;
  }
}

data {
  int iter;
  
  int Nrec;
  int warmup;
  int thin;
  
  int Ni;
  real t0;
  real dt;
  real C0_LD;
  real C0_SD;
  //parameters
  real a[iter];
  real cdwn[iter];
  real h[iter];
  real<lower=0,upper=1> tL_LD[iter];
  real<lower=0,upper=1> tL_SD[iter];
  real r[iter];
}

transformed data {
  //fixed paraameters
  real dusk_LD;
  real dusk_SD;
  real lambda;
  dusk_LD = 16.0/24.0;
  dusk_SD = 8.0/24.0;
  lambda = 1.0;
}

model {
}

generated quantities {
  real Y_LD[4*iter,Nrec];
  real Y_SD[4*iter,Nrec];
  for(i in 1:iter){
    real y_LD[4,Ni];
    real y_SD[4,Ni];
    y_LD = rungekutta(t0, dt, Ni, C0_LD, a[i], cdwn[i], dusk_LD, h[i], lambda, tL_LD[i], r[i]);
    y_SD = rungekutta(t0, dt, Ni, C0_SD, a[i], cdwn[i], dusk_SD, h[i], lambda, tL_SD[i], r[i]);
    for(j in 1:Nrec){
      int J;
      J = warmup + thin*j;
      Y_LD[(4*i-3):(4*i),j] = y_LD[1:4,J];
      Y_SD[(4*i-3):(4*i),j] = y_SD[1:4,J];
    }
  }
}
