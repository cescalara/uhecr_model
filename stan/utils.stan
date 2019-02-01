/**
 * Utils for use with Stan simulations and models.
 *
 * @author Francesca Capel
 * @date October 2018
 */
  
/**
 * Calculate weights from source distances.
 */
vector get_source_weights(real[] L, real[] D) {
  
  int N = num_elements(D);
  vector[N] weights;
  
  real normalisation = 0;
  
  for (k in 1:N) {
    normalisation += (L[k] / pow(D[k], 2));
  }
  for (k in 1:N) {
    weights[k] = (L[k] / pow(D[k], 2)) / normalisation;
  }
  
  return weights;
}

/**
 * Calculate weights for each source accounting for exposure 
 * and propagation effects.
 */
vector get_exposure_weights(vector F, vector eps, real alpha_T, vector Eth_src, real Eth, real alpha) {
  
  int N = num_elements(F);
  vector[N] weights;
  
  real normalisation = 0;
  
  for (k in 1:N-1) {
    normalisation += F[k] * eps[k] * pow(Eth_src[k] / Eth, 1 - alpha);
  }
  normalisation += F[N] * (alpha_T / (4 * pi()));
  
  for (k in 1:N-1) {
    weights[k] = (F[k] * eps[k] * pow(Eth_src[k] / Eth, 1 - alpha)) / normalisation;
  }
  weights[N] = (F[N] * (alpha_T / (4 * pi()))) / normalisation;
  
  return weights;
}

/**
 * Convert from unit vector omega to theta of spherical coordinate system.
 * @param omega a 3D unit vector.
 */
real omega_to_theta(vector omega) {
  
  real theta;
  
  int N = num_elements(omega);
  
  if (N != 3) {
    print("Error: input vector omega must be of 3 dimensions");
  }
  
  theta = acos(omega[3]);
  
  return theta;
}

/**
 * Calculate the expected value of N for the generative model.
 */
real get_Nex_sim(vector F, vector eps, real alpha_T, vector Eth_src, real Eth, real alpha) {
  
  int N = num_elements(F);
  real Nex = 0;
  
  for (k in 1:N-1) {
    Nex += F[k] * eps[k] * pow(Eth_src[k] / Eth, 1 - alpha);
  }
  Nex += F[N] * (alpha_T / (4 * pi()));

  return Nex;
}
  
/**
 * Calculate the total source flux.
 * @param L the luminosity in s^-1
 * @param D the distance in Mpc
 */
real get_Fs(real[] L, real[] D) {
  
  int N = num_elements(D);
  real Fs = 0;

  for (k in 1:N) {
    Fs += L[k] / (4 * pi() * pow(D[k], 2));
  }
  
  return Fs;
}

/**
 * Sample from the vMF centred on varpi with spread kappa, 
 * accounting for detector exposure effects.
 * Uses rejection sampling.
 */
vector exposure_limited_vMF_rng(vector varpi, real kappa, real a0, real theta_m) {
  
  real params[3];
  real m_max;
  real accept;
  int count;
  vector[3] omega;
  real theta;
  real pdet;
  vector[2] p;
  
  /* exposure */
  params[1] = cos(a0);
  params[2] = sin(a0);
  params[3] = cos(theta_m);
  m_max = m(pi(), params);
  
  accept = 0;
  count = 0;

  while (accept != 1) {	  

    omega = vMF_rng(varpi, kappa);
    theta = omega_to_theta(omega);
    pdet = m(theta, params) / m_max;
    p[1] = pdet;
    p[2] = 1 - pdet;
    accept = categorical_rng(p);
    count += 1;
    if (count > 1.0e7) {
      
      print("Was stuck in exposure_limited_rng");
      accept = 1;

    }
    
  } 
  
    return omega;
  }

/**
 * Sample uniformly from the unit sphere, 
 * accounting for detector exposure effects.
 * Uses rejection sampling.
 */
vector exposure_limited_sphere_rng(real a0, real theta_m) {
  
  real params[3];
  real m_max;
  real accept;
  vector[3] omega;
  real theta;
  real pdet;
  vector[2] p;
  
  params[1] = cos(a0);
  params[2] = sin(a0);
  params[3] = cos(theta_m);
  m_max = m(pi(), params);
  
  accept = 0;
  while (accept != 1) {	  

    omega = sphere_rng(1);
    theta = omega_to_theta(omega);
    pdet = m(theta, params) / m_max;
    p[1] = pdet;
    p[2] = 1 - pdet;
    accept = categorical_rng(p);

  } 

    return omega;
  }


/**
 * Interpolate x from a given set of x and y values.
 */
real interpolate(vector x_values, vector y_values, real x) {
  real x_left;
  real y_left;
  real x_right;
  real y_right;
  real dydx;
  
  int Nx = num_elements(x_values);
  real xmin = x_values[1];
  real xmax = x_values[Nx];
  int i = 1;
  
  if (x > xmax || x < xmin) {
    
    if(x > xmax) {
      return y_values[Nx];
    }
    else if (x < xmin) {
      return y_values[1];
    }
  }
  
  if( x >= x_values[Nx - 1] ) {
    i = Nx - 1;
  }
  else {
    while( x > x_values[i + 1] ) { i += 1; }
  }
  
  x_left = x_values[i];
  y_left = y_values[i];
  x_right = x_values[i + 1];
  y_right = y_values[i + 1];
  
  dydx = (y_right - y_left) / (x_right - x_left);
  
  return y_left + dydx * (x - x_left);
}

/**
 * Calculate the Nex for a given kappa by
 * interpolating over a vector of eps values
 * for each source.
 */
real get_Nex(vector F, vector[] eps, vector kappa_grid, vector kappa, real alpha_T, vector Eth_src, real Eth, real alpha) {
  
  int Ns = num_elements(F);
  vector[Ns] N;
  real eps_from_kappa;
  
  for (k in 1:Ns-1) {
    eps_from_kappa = interpolate(kappa_grid, eps[k], kappa[k]);
    N[k] = F[k] * eps_from_kappa * pow(Eth_src[k] / Eth, 1 - alpha); 
  }
  N[Ns] = F[Ns] * (alpha_T / (4 * pi()));
  
  return sum(N);
}

/**
 * Calculate the Nex for a given kappa by
 * interpolating over a vector of eps values
 * for each source.
 * This version ignores energy information and is 
 * to be used in the arrival direction only models.
 */
real get_Nex_arr(vector F, vector[] eps, vector kappa_grid, real kappa, real alpha_T) {
  
  int Ns = num_elements(F);
  vector[Ns] N;
  real eps_from_kappa;
  
  for (k in 1:Ns-1) {
    eps_from_kappa = interpolate(kappa_grid, eps[k], kappa);
    N[k] = F[k] * eps_from_kappa; 
  }
  N[Ns] = F[Ns] * (alpha_T / (4 * pi()));
  
  return sum(N);
}


/**
 * Define the fik PDF.
 * NB: Cannot be vectorised.
 * Uses sinh(kappa) ~ exp(kappa)/2 
 * approximation for kappa > 100.
 */
real fik_lpdf(vector v, vector mu, real kappa, real kappa_d) {
  
  real lprob;
  real inner = abs_val((kappa_d * v) + (kappa * mu));
  
  if (kappa > 100 || kappa_d > 100) {
    lprob = log(kappa * kappa_d) - log(4 * pi() * inner) + inner - (kappa + kappa_d) + log(2);
  }
  else {   
    lprob = log(kappa * kappa_d) - log(4 * pi() * sinh(kappa) * sinh(kappa_d)) + log(sinh(inner)) - log(inner);
  }
  
  return lprob;   
}

/**
 * Calculate the deflection parameter of the vMF distribution.
 * @param E energy in EeV
 * @param B rms magnetic field strength in nG
 * @param D distance in Mpc / 10
 */
real get_kappa(real E, real B, real D) {
  
  return 2.3 * inv_square( 0.0401 * inv(E / 50) * B * sqrt(D) );
}

/**
 * Calculate the vector of kappa values for different sources.
 * @param E energy in EeV
 * @param B rms magnetic field strength in nG
 * @param D distance in Mpc / 10
 */
vector get_kappa_ex(vector E, real B, vector D) {
  
  int Ns = num_elements(E);
  vector[Ns] kappa_ex = 2.3 * inv_square( 0.0401 * inv(E / 50) * B .* sqrt(D) );
  return kappa_ex;
}

/**
 * Calculate the vector of expected E values for different sources.
 */
vector get_Eex(real alpha, vector Eth_src) {
  
  int N = num_elements(Eth_src);
  vector[N] Eex = pow(2, 1 / (alpha - 1)) * Eth_src[1:N];
  
  return Eex;
}

vector get_eps_from_kappa(vector kappa_grid, vector[] eps, vector kappa_ex) {
  
  int N = num_elements(kappa_ex);
  vector[N] eps_from_kappa;
  
  for (k in 1:N) {
    eps_from_kappa[k] = interpolate(kappa_grid, eps[k], kappa_ex[k]);
  }
  
  return eps_from_kappa;
}  
