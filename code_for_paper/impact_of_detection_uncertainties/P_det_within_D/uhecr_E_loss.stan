/**
 * Simulation of UHECR energy losses using the 
 * continuous loss approximation (see uhecr_propagation.stan)
 * 
 * Used to investigation impact of detection uncertainties
 */


functions {

#include uhecr_propagation.stan
  
}

data {

  int<lower=0> N;
  real alpha;
  real<lower=0> Eth_sim;
  real<lower=0> f_E;
  real D;
  
}

transformed data {

  /* for integral in beta_bh calc */
  real x_r[1];
  int x_i[0];
  real params[0];
  real D_in[1];
  real Eth_src;
  
  /* approximate infinity */
  x_r[1] = 1.0e4;

  /* find corresponding energy threshold at the source */
  D_in[1] = D;
  if (D == 0) {
    Eth_src = Eth_sim;
  }
  else {
    Eth_src = get_source_threshold_energy(Eth_sim, D_in, x_r, x_i);
  }
  
}

generated quantities {

  real E[N];
  real Earr[N];
  real Edet[N];

  /* loop over UHECR */
  for (i in 1:N) {

    /* sample from power law spectrum */
    E[i] = spectrum_rng(alpha, Eth_src);

    /* propagate */
    if (D == 0) {
      
      Earr[i] = E[i];

    }
    else {
      
      Earr[i] = get_arrival_energy(E[i], D_in, x_r, x_i);

    }

    /* detection effects */
    Edet[i] = normal_rng(Earr[i], f_E * Earr[i]);
      
  }
  
}
