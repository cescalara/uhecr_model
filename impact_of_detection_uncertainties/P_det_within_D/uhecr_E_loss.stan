/**
 * Simple simulation to test out energy modelling.
 * Propagation losses based on Khanin et al. 2016.
 * Investigating sample selection effects.
 */

functions {


  /**
   * Shape of the energy spectrum: a power law.
   */
  real dNdE_pl(real E, real alpha) {

    real spec = pow(E, -alpha);
    return spec; 

  }
    
  /**
   * Sample an energy from a spectrum defined by alpha.
   * Uses rejection sampling.
   * Sampled energy is in units of EeV.
   */
  real spectrum_rng(real alpha, real Eth) {

    real E;
    real d;
    real d_upp_lim = dNdE_pl(Eth, alpha);

    int accept = 0;
    while(accept != 1) {
      E = uniform_rng(Eth, 1e4);
      d = uniform_rng(0, d_upp_lim);
          
      if (d < dNdE_pl(E, alpha)) {
	accept = 1;
      }
    }

    return E;
  }


  /**
   * Analytic approx of GZK losses due to photomeson production.
   * Originally parametrised by Anchordorqui et al. (1997)
   * Based on Berezinsky and Grigor'eva (1988) suggestion of Aexp(-B/E) form.
   * @param z the redshift
   * @apram E the energy in eV
   */
  real beta_pi(real z, real E) {

    real output;
    real check = 6.86 * exp(-0.807 * z) * 1e20;
    real p[3];
    p[1] = 3.66e-8;
    p[2] = 2.87e20;
    p[3] = 2.42e-8;

    if (E <= check) {
      output = p[1] * pow(1 + z, 3) * exp(-p[2] / ((1 + z) * E));
    }
    
   if (E > check) {
      output = p[3] * pow(1 + z, 3);
    }
   
   return output;
  }

  /**
   * Losses due to adiabatic expansion. 
   * Makes use of the hubble constant converted into units of [yr^-1].
   * @param z the redshift
   */
  real beta_adi(real z) {

    real lCDM[2];
    real a;
    real b;
    real H0 = 70.4; // s^-1 km/Mpc
    H0 = ((H0 / 3.086e22) * 1e3) * 3.154e7; // yr^-1
    lCDM[1] = 0.272;
    lCDM[2] = 0.728;
    
    a = lCDM[1] * pow(1 + z, 3);
    b = 1 - sum(lCDM) * pow(1 + z, 2);
    
    return H0 * pow(a + lCDM[2] + b, 0.5);
  }

  /**
   * Approximation of the phi(xi) function as xi->inf.
   * Used in calculation of beta_bh.
   * Described in Chodorowski et al. (1992).
   */
  real phi_inf(real xi) {

    real d[4]; 
    real sum_term = 0;
    d[1] = -86.07;
    d[2] = 50.96;
    d[3] = -14.45;
    d[4] = 8.0 / 3.0;

    for (i in 1:num_elements(d)) {
      sum_term += d[i] * pow(log(xi), i - 1);
    }

    return xi * sum_term;
  }

  /**
   * Approximation of the phi(xi) function for different regimes
   * of xi.
   * Used in calculation of beta_bh.
   * Described in Chodorowski et al. (1992). 
   */ 
  real phi(real xi) {

    real output;
    real sum_term;
    real phi_inf_term;
    real c[4];
    real f[3];
    
    if (xi == 2) { 
      output = (pi() / 12); // * pow(xi - 2, 4);
    }
    
    else if (xi < 25) {

      c[1] = 0.8048;
      c[2] = 0.1459;
      c[3] = 1.137e-3;
      c[4] = -3.879e-6;
      sum_term = 0;

      for (i in 1:num_elements(c)) {
	  sum_term += c[i] * pow(xi - 2, i); 
	}
      
      output = (pi() / 12) * pow(xi - 2, 4) / (1 + sum_term);
    }
    
    else if (xi >= 25) {

      f[1] = 2.910;
      f[2] = 78.35;
      f[3] = 1837;
	
      phi_inf_term = phi_inf(xi);
      sum_term = 0;
      for (i in 1:num_elements(f)) {
	sum_term += f[i] * pow(xi, -i);
      }

      output = phi_inf_term / (1 - sum_term);
    }
    
    return output;
  }

  real phi_integrand(real xi, real E, real z) {

    real num = phi(xi);
    real B = 1.02;
    real denom = exp( (B * xi) / ((1 + z) * E) ) - 1;

    return num / denom;
  }
  
  /**
   * Integrand as a functional of phi(xi) used in calcultion of beta_bh.
   * Described in de Domenico and Insolia (2013).
   */
  real[] integrand(real xi, real[] state, real[] params, real[] x_r, int[] x_i) { 

    real dstatedxi[1];

    real E = params[1];
    real z = params[2];

    dstatedxi[1] = phi_integrand(xi, E, z);
    return dstatedxi;
  }

  /**
   * Losses due to Bethe-Heitler pair production.
   * Described in de Domenico amnd Insolia (2013).    
   */
  real beta_bh(real z, real E, real[] xiout, real[] x_r, int[] x_i) {
    
    real params[2];
    real integration_result[1,1];
    real state0[1];
    real integ;
    real A = 3.44e-18;  
    
    state0[1] = phi_integrand(2.0, E, z);
    if (is_nan(state0[1])) {
      print("phi_integrand(2.0, E, z) is nan");
      print("E: ", E);
      print("z: ", z);
    }
    
    params[1] = E;
    params[2] = z;

    integration_result = integrate_ode_rk45(integrand, state0, 2.0, xiout, params, x_r, x_i);
    integ = integration_result[1,1];

    return (A  / pow(E, 3)) * integ;
  }

  /**
   * Total energy losses as a function of distance.
   * @param E the energy in eV (!)
   */
  real dEdr(real r, real E, real[] xiout, real[] x_r, int[] x_i) {

    real c = 3.066e-7; // Mpc/yr
    real DH = 4285.714; // Mpc
    real z = r / DH;
    real beta_pp;
    real Ltot;

    beta_pp = 3.154e7 * beta_bh(z, E / 1.0e18, xiout, x_r, x_i);
    Ltot = c / (beta_adi(z) + beta_pi(z, E) + beta_pp);
  
    return - E / Ltot;
    
  }

  /**
   * ODE system to be solved for arrival energy.
   * NB: E in eV (!)
   */
  real[] E_ode(real r, real[] state, real[] params, real[] x_r, int[] x_i) {

    real dstatedr[1];

    real E = state[1];
    
    dstatedr[1] = dEdr(r, E, x_r, x_r, x_i);
    return dstatedr;
    
  }

  /**
   * ODE system for be solved for equivalent source energy.
   * Nb: E in eV (!)
   */
  real[] E_ode_rev(real r, real[] state, real[] params, real[] x_r, int[] x_i) {

    real dstatedr[1];
    real D = params[1];
    real r_rev = D - r;
    real E = state[1];
    
    dstatedr[1] = - dEdr(r_rev, E, x_r, x_r, x_i);
    return dstatedr;
    
  }

  /**
   * Calculate equivalent energy at the source for an arrival energy of Eth.
   * Solves the ODE system  dE/dr = E / Lloss.
   */
  real get_source_threshold_energy(real Eth, data real[] D, data real[] x_r, int[] x_i) {

    real Eth_src;
    real params[1];
    real integration_result[1, 1];
    real Eth_in[1];

    Eth_in[1] = Eth * 1e18; // eV
    params[1] = D[1];
    
    integration_result = integrate_ode_rk45(E_ode_rev, Eth_in, 0.0, D, params, x_r, x_i);
    Eth_src = integration_result[1, 1] / 1e18; // EeV   

    return Eth_src;
  }

  /**
   * Get the vector of source threshold energies for all sources, 
   * including the background component.
   */
  vector get_Eth_src(real Eth, data real[,] D, data real[] x_r, int[] x_i) {

    int N = num_elements(D);
    vector[N] Eth_src;
    
    for (k in 1:N) {
      Eth_src[k] = get_source_threshold_energy(Eth, D[k], x_r, x_i);
    }    

    return Eth_src;

  }

  /**
   * Calculate the arrival energy taking into account all propagation losses
   * for a given intial energy E.
   * Solves the ODE system dE/dr = - E / Lloss
   */
  real get_arrival_energy(real E, data real[] D, data real[] x_r, int[] x_i) {
    
    real Earr;
    real params[0];
    real integration_result[1, 1];
    real E_in[1];

    E_in[1] = E * 1e18; // eV

    integration_result = integrate_ode_rk45(E_ode, E_in, 0.0, D, params, x_r, x_i);
    Earr = integration_result[1, 1] / 1e18; // EeV   

    return Earr;   
  }
  
}

data {

  int<lower=0> N;
  real alpha;
  real<lower=0> Eth_sim;
  real D;

}

transformed data {

  int Nobs;

  /* for integral in beta_bh calc */
  real x_r[1];
  int x_i[0];
  real params[0];
  real D_in[1];
  real Eth_src;
  
  /* approximate infinity */
  x_r[1] = 1.0e4;

  D_in[1] = D;
  if (D == 0) {
    Eth_src = Eth_sim;
  }
  else {
    Eth_src = get_source_threshold_energy(Eth_sim, D_in, x_r, x_i);
  }
  Nobs = N;
  
}

generated quantities {

  real E[Nobs];
  real Earr[Nobs];
  real Edet[Nobs];
  
  for (i in 1:Nobs) {

    E[i] = spectrum_rng(alpha, Eth_src);
    if (D == 0) {
      Earr[i] = E[i];
    }
    else {
      Earr[i] = get_arrival_energy(E[i], D_in, x_r, x_i);
    }
    Edet[i] = normal_rng(Earr[i], 0.12 * Earr[i]);
      
  }
  
}
