/**
 * Joint model for UHECR energies and arrival directions.
 * Generative model (i.e. forward simulation)
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

  /* source positions */
  int<lower=0> Ns;
  unit_vector[3] varpi[Ns];
  real D[Ns];

  /* source spectrum */
  real alpha;
  real<lower=0> Eth;
  real<lower=0> Eerr;
  
  /* flux */
  real<lower=0> L[Ns];
  real<lower=0> F0;
  
  /* deflection */
  real<lower=0> B;
  real<lower=0> kappa_d;

  /* observatory parameters */
  real A;
  real a0;
  real<lower=0> theta_m;
  real<lower=0> alpha_T;
  vector[Ns] eps;
  
}

transformed data {
  
  /* definitions */
  real<lower=0> Fs = get_Fs(L, D);
  real<lower=0> F_T = F0 + Fs;
  real<lower=0, upper=1> f = Fs / F_T;
  simplex[Ns] w = get_source_weights(L, D);  
  vector[Ns+1] F;

  real x_r[1];
  int x_i[0];
  vector[Ns] Eth_src;
  real D_in[Ns, 1];
  vector[Ns] D_kappa;
  
  simplex[Ns+1] w_exposure;
  real<lower=0> Nex;
  int<lower=0> N;

  /* flux and distance */
  for (k in 1:Ns) {
    F[k] = w[k] * Fs;
    D_in[k, 1] = (D[k] / 3.086) * 100; // Mpc
  }
  F[Ns+1] = F0;

  /* D in Mpc / 10 for kappa calculation */
  for (k in 1:Ns) {
    D_kappa[k] = (D[k] / 3.086) * 10; // Mpc / 10
  }
 
  /* Eth_src */
  x_r[1] = 1.0e4; // approx inf
  Eth_src = get_Eth_src_sim(Eth, D_in, x_r, x_i);
  
  /* N */
  w_exposure = get_exposure_weights(F, eps, alpha_T, Eth_src, Eth, alpha);

  Nex = get_Nex_sim(F, eps, alpha_T, Eth_src, Eth, alpha);
  
  N = poisson_rng(Nex);

}

generated quantities {

  int lambda[N];
  unit_vector[3] omega;
  unit_vector[3] arrival_direction[N];
  real Nex_sim = Nex;
  
  real E[N];
  real kappa[N];
  real Earr[N];
  real Edet[N];
  
  for (i in 1:N) {

    /* label */
    lambda[i] = categorical_rng(w_exposure);

    /* source */
    if (lambda[i] < Ns+1) {

      E[i] = spectrum_rng(alpha, Eth_src[lambda[i]]);
      kappa[i] = get_kappa(E[i], B, D_kappa[lambda[i]]);
      omega = exposure_limited_vMF_rng(varpi[lambda[i]], kappa[i], a0, theta_m);      
      Earr[i] = get_arrival_energy_sim(E[i], D_in[lambda[i]], x_r, x_i);

    }
  
    /* background */
    if (lambda[i] == Ns+1) {

      E[i] = spectrum_rng(alpha, Eth);
      omega = exposure_limited_sphere_rng(a0, theta_m);
      Earr[i] = E[i];

    }
    
    /* detection */
    arrival_direction[i] = vMF_rng(omega, kappa_d);  	  
    Edet[i] = normal_rng(Earr[i], Eerr * Earr[i]);

  }
  
}

