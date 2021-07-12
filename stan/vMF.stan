/**
 * Functions for sampling and fitting the vMF distribution in Stan.
 *
 * @author Francesca Capel
 * @date May 2018
 */

/**
 * Compute the absolute value of a vector. 
 */
real abs_val(vector input_vector) {
  
  real av;
  int n = num_elements(input_vector);
  
  real sum_squares = 0;
  for (i in 1:n) {
    sum_squares += (input_vector[i] * input_vector[i]);
  }
  av = sqrt(sum_squares);
  return av;
  
}

/**
 * Sample point on sphere orthogonal to mu.
 */
vector sample_orthonormal_to_rng(vector mu) {
  
  int dim = num_elements(mu);
  vector[dim] v;
  vector[dim] proj_mu_v;
  vector[dim] orthto;
  
  for (i in 1:dim) {
    v[i] = normal_rng(0, 1);
  }
  
  proj_mu_v = mu * dot_product(mu, v) / abs_val(mu);
  orthto = v - proj_mu_v;
  
  return (orthto / abs_val(orthto));
  
}

/**
 * Rejection sampling scheme for sampling distance from center on
 * surface of the sphere.
 */
real sample_weight_rng(real kappa, int dim) {
  
  real sdim = dim - 1; /* as S^{n-1} */
  real b = sdim / (sqrt(4. * pow(kappa, 2) + pow(sdim, 2)) + 2 * kappa);
  real x = (1 - b) / (1 + b);
  real c = kappa * x + sdim * log(1 - pow(x, 2));
  
  int i = 0;
  real z;
  real w;
  real u;

  while (i == 0) {

    z = beta_rng(sdim / 2, sdim / 2);
    w = (1 - (1 + b) * z) / (1 - (1 - b) * z);
    u = uniform_rng(0, 1);

    if (kappa * w + sdim * log(1 - x * w) - c >= log(u)) {
      i = 1;
    }

  }
  
  return w;
}

/**
 * Generate an N-dimensional sample from the von Mises - Fisher
 * distribution around center mu in R^N with concentration kappa.
 * 
 * KW: Based on the following code:
 * https://github.com/jasonlaska/spherecluster/blob/develop/spherecluster/util.py
 */
vector vMF_rng(vector mu, real kappa) {
  
  int dim = num_elements(mu);
  vector[dim] result;
  
  real w = sample_weight_rng(kappa, dim);
  vector[dim] v = sample_orthonormal_to_rng(mu);
  
  result = ( v * sqrt(1 - pow(w, 2)) ) + (w * mu);
  return result;
  
}

/**
 * Sample a point uniformly from the surface of a sphere of 
 * a certain radius.
 */
vector sphere_rng(real radius) {
  
  vector[3] result;
  real u = uniform_rng(0, 1);
  real v = uniform_rng(0, 1);
  real theta = 2 * pi() * u;
  real phi = acos( (2 * v) - 1 );
  
  result[1] = radius * cos(theta) * sin(phi); 
  result[2] = radius * sin(theta) * sin(phi); 
  result[3] = radius * cos(phi);
  
  return result;
  
}
