/**
 * Functions for modelling the exposure of a UHECR observatory.
 *
 * @author Francesca Capel
 * @date May 2018
 */


/**
 * Calculate xi part of exposure.
 * @param theta from 0 to pi.
 * @param p observatory dependent parameters.
 */
real xi_exp(real theta, real[] p) { 
  return (p[3] - (p[2] * cos(theta))) / (p[1] * sin(theta));
}

/**
 * Calculate alpha_m part of exposure.
 * @param theta from 0 to pi.
 * @param p observatory dependent parameters.
 */
real alpha_m(real theta, real[] p) {
  
  real am;
  
  real xi_val = xi_exp(theta, p);
  if (xi_val > 1) {
    am = 0;
  }
  else if (xi_val < -1) {
    am = pi();
  }
  else {
    am = acos(xi_val);
  }
  
  return am;
}

/**
 * Calculate the exposure factor for a given position on the sky. 
 * @param theta from 0 to pi.
 * @param p observatory dependent parameters.
 */
real m(real theta, real[] p) {
  return (p[1] * sin(theta) * sin(alpha_m(theta, p)) 
	  + alpha_m(theta, p) * p[2] * cos(theta));
}
