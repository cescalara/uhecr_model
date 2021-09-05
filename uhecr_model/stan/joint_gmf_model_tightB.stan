/**
 * Joint model for UHECR energies and arrival directions.
 * Has one background component and option to use full 
 * energy calculation or interpolation.
 *
 * @author Francesca Capel
 * @date October 2018
 */

functions {

#include energy_spectrum.stan
#include uhecr_propagation.stan
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
  real<lower=0> Edet[N];
  vector[N] zenith_angle;
  vector[N] A;
  vector[N] kappa_gmf;  
  int Z;
  
  /* observatory */ 
  real<lower=0> alpha_T;
  int Ngrid;
  vector[Ngrid] eps[Ns];
  vector[Ngrid] kappa_grid;

  /* energy */
  real<lower=0> Eth;
  real<lower=0> Eerr;
  vector[Ngrid] Earr_grid[Ns+1];
  vector[Ngrid] E_grid;
  
}

transformed data {

  real x_r[0];
  int x_i[0];
  vector[Ns] Eth_src;
  real D_in[Ns, 1];
  vector[Ns] D_kappa;
  
  /* D in Mpc for ODE solver */
  for (k in 1:Ns) {
    D_in[k, 1] = (D[k] / 3.086) * 100; // Mpc
  }

  /* D in Mpc / 10 for kappa calculation */
  for (k in 1:Ns) {
    D_kappa[k] = (D[k] / 3.086) * 10; // Mpc / 10
  }
  
  /* Equivalent Eth at sources */
  Eth_src = get_Eth_src(Eth, D_in, x_r, x_i);
  
}


parameters { 

  /* source luminosity */
  real<lower=0, upper=(1e5 / Ns)> L;
  
  /* background flux */
  real<lower=0, upper=1e3> F0;
  
  /* magnetic field strength */
  real<lower=1, upper=1e2> B;
  
  /* energy spectrum */
  real<lower=1, upper=10> alpha;  
  vector<lower=Eth, upper=1e4>[N] E;

}


transformed parameters {

  /* flux */
  real<lower=0> Fs;
  vector[Ns+1] F;
  real<lower=0> FT;
  
  /* associated fraction */
  real<lower=0, upper=1> f; 
    
  /* association probability */
  vector[Ns+1] lp[N];
  real Earr[N];
  vector[N] kappa;
  vector[Ns+1] log_F;
  
  /* Nex */
  vector[Ns] Eex;
  vector[Ns] kappa_ex;
  real<lower=0> Nex;

  /* define transformed paramaters */
  for (k in 1:Ns) {
    F[k] = L / (4 * pi() * pow(D[k], 2));
  }
  Fs = sum(F[1:Ns]); 
  F[Ns+1] = F0;
  log_F = log(F);

  FT = F0 + Fs;
  f = Fs / FT;

  /* likelihood calculation */
  /* rate factor */
  for (i in 1:N) {

    lp[i] = log_F;
 
    for (k in 1:Ns+1) {

      /* sources */
      if (k < Ns+1) {

	kappa[i] = get_kappa(E[i], B, D_kappa[k], Z);
	lp[i, k] += fik_lpdf(arrival_direction[i] | varpi[k], kappa[i], kappa_gmf[i]);

	/* choose full energy calculation or interpolation for speed */
	//Earr[i] = get_arrival_energy(E[i], D_in[k], x_r, x_i); // full calc
	Earr[i] = interpolate(E_grid, Earr_grid[k], E[i]); // interpolation
	
	lp[i, k] += pareto_lpdf(E[i] | Eth, alpha - 1);
      
      }
      
      /* background */
      if (k == Ns + 1) {

	lp[i, k] += log(1 / ( 4 * pi() ));
	Earr[i] = E[i];
	lp[i, k] += pareto_lpdf(E[i] | Eth, alpha - 1);
      
      }

      /* truncated gaussian */
      lp[i, k] += normal_lpdf(Edet[i] | Earr[i], Eerr * Earr[i]);
      if (Edet[i] < Eth) {
	lp[i, k] += negative_infinity();
      }
      else {
	lp[i, k] += -normal_lccdf(Eth | Earr[i], Eerr * Earr[i]);
      }

      /* exposure factor */
      lp[i, k] += log(A[i] * cos(zenith_angle[i]));
            
    }
  }

  /* Nex */
  Eex = get_Eex(alpha, Eth_src);
  kappa_ex = get_kappa_ex(Eex, B, D_kappa, Z);
  Nex = get_Nex(F, eps, kappa_grid, kappa_ex, alpha_T, Eth_src, Eth, alpha); 
  
}


model {

  /* rate factor */
  for (i in 1:N) {
    target += log_sum_exp(lp[i]);
  }
  
  /* normalise */
  target += -Nex; 

  /* priors */
  alpha ~ normal(3, 2);
  B ~ normal(1, 3);
  L ~ normal(0, 1.0e3 / Ns);
  F0 ~ normal(0, 1.0e2);

}

generated quantities {

  int<lower=1, upper=Ns+1> lambda[N];

  /* used in calculating the source-UHECR association probabilities */
  for (i in 1:N) {

    lambda[i] = categorical_logit_rng(lp[i]);

  }

}

