/**
 * Model for the UHECR arrival directions.
 *
 * @author Francesca Capel
 * @date August 2018
 */

functions {

#include vMF.stan
#include observatory_exposure.stan
#include utils.stan
  
}

data {

  /* sources */
  int<lower=0> Ns;
  unit_vector[3] varpi[Ns]; 
  vector[Ns] D;
  
  /* uhecr */
  int<lower=0> N; 
  unit_vector[3] arrival_direction[N]; 
  vector[N] zenith_angle;
  vector[N] A;
  
  /* observatory */
  real<lower=100, upper=10000> kappa_d;  
  real<lower=0> alpha_T;
  int Ngrid;
  vector[Ngrid] eps[Ns];
  vector[Ngrid] kappa_grid;
  
}

parameters { 

  /* source luminosity */
  real<lower=0, upper=(1e5 / Ns)> L;
  
  /* background flux */
  real<lower=0, upper=1e3> F0;

  /* deflection */
  real<lower=1, upper=1000> kappa;  
  
}

transformed parameters {

  /* total source flux */
  real<lower=0> Fs;   
  
  /* source flux */
  vector[Ns + 1] F;

  /* associated fraction */
  real<lower=0, upper=1> f; 

  real<lower=0> FT;

  Fs = 0;
  for (k in 1:Ns) {
    F[k] = L / (4 * pi() * pow(D[k], 2));
    Fs += F[k];
  }
  F[Ns + 1] = F0;

  FT = (F0 + Fs);

  f = Fs / FT;
 
}

model {

  vector[Ns + 1] log_F;
  real Nex;
  
  log_F = log(F);

  /* Nex */
  Nex = get_Nex_arr(F, eps, kappa_grid, kappa, alpha_T);

  /* rate factor */
  for (i in 1:N) {

    vector[Ns + 1] lps = log_F;

    for (k in 1:Ns + 1) {

      if (k < Ns + 1) {
	lps[k] += fik_lpdf(arrival_direction[i] | varpi[k], kappa, kappa_d);
      }
      else {
	lps[k] += log(1 / ( 4 * pi() ));
      }
      
      lps[k] += log(A[i] * cos(zenith_angle[i]));
    }
    
    target += log_sum_exp(lps);

  }
  
  /* normalise */
  target += -Nex; 

  /* priors */
  kappa ~ lognormal(log(100.), 1.);

  L ~ normal(0, 1.0e3 / Ns);
  F0 ~ normal(0, 1.0e2);

}
