/**
 * Energy spectrum sampling.
 * 
 * @author Francesca Capel
 * @date June 2018
 */

/**
 * Shape of the energy spectrum: a power law.
 */
real dNdE_pl(real E, real alpha) {
  
  real spec = pow(E, -alpha);
  return spec; 
  
}

/**
 * Sample an energy from a power law spectrum defined by alpha.
 * Uses rejection sampling.
 * Sampled energy is in units of EeV.
 * @param alpha the spectral index
 * @param Emin the minimum energy
 */  
real spectrum_rng(real alpha, real Emin) {
  
  real E;
  real d;
  real d_upp_lim = dNdE_pl(Emin, alpha);

  int accept = 0;
  
  while(accept != 1) {
    
    E = uniform_rng(Emin, 1e4);
    d = uniform_rng(0, d_upp_lim);
    
    if (d < dNdE_pl(E, alpha)) {
      accept = 1;
    }
  }

  return E;
}

/**
 * Log pdf for the power law distribution, bounded
 * below by Emin.
 * @param E energy in EeV
 * @param alpha spectral index
 * @param Emin minimum energy in EeV
 */
real spectrum_lpdf(real E, real alpha, real Emin) {
  
  real norm = log(alpha - 1) - log(Emin);
  real lprob = -alpha * log(E / Emin);
  
  return lprob + norm;
}
