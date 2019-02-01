/**
 * Simple simulation to show model effects.  
 *
 * @author Francesca Capel
 * @date November 2018
 */

functions {
  
#include energy_spectrum.stan
#include uhecr_propagation.stan
#include vMF.stan
#include observatory_exposure.stan
#include utils.stan
  
}


data {

  /* source position */
  unit_vector[3] varpi;
  
  /* source spectrum */
  real alpha;
  real<lower=0> Eth; // EeV
  real<lower=0> Eerr; // EeV
  
  /* flux */
  int<lower=0> N;
 
  /* deflection */
  real<lower=0> B; // nG
  real<lower=0> kappa_d;

}

transformed data {
  
  real x_r[1];
  int x_i[0];

  /* approximate infinity */
  x_r[1] = 1.0e4;
 
}

generated quantities {

  unit_vector[3] omega[N];
  unit_vector[3] omega_det[N];
  
  real E[N];
  real kappa[N];
  real Earr[N];
  real Edet[N];
  real Eth_src[N];
  real D[N];
  real D_in[N, 1];
  
  for (i in 1:N) {

    /* Uniformly distributed sources */
    D[i] = uniform_rng(0, 400);
    D_in[i, 1] = D[i];
    Eth_src[i] = get_source_threshold_energy_sim(Eth, D_in[i], x_r, x_i);

    /* Energy sampling and propagation */
    E[i] = spectrum_rng(alpha, Eth_src[i]);
    kappa[i] = get_kappa(E[i], B, D[i] / 10);
    omega[i] = vMF_rng(varpi, kappa[i]);      
    Earr[i] = get_arrival_energy_sim(E[i], D_in[i], x_r, x_i);
      
    /* Detection effects */
    omega_det[i] = vMF_rng(omega[i], kappa_d);  	  
    Edet[i] = normal_rng(Earr[i], Eerr * Earr[i]);

  }
  
}

