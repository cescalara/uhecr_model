/**
 * Model for UHECR arrival directions.
 * Reparametrised in terms of the source luminosity and
 * background flux.
 * Some extensions in order to faciliate easy comparison with the
 * current joint_energy_arrival models.
 * @author Francesca Capel
 * @date August 2018
 */


functions {

  /**
   * compute the absolute value of a vector
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
      print("Warning, x is outside of interpolation range!");
      print("Returning edge values.");
      print("x:", x);
      print("xmax", xmax);
      
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
      while( x > x_values[i + 1] ) { i = i+1; }
    }

    x_left = x_values[i];
    y_left = y_values[i];
    x_right = x_values[i + 1];
    y_right = y_values[i + 1];

    dydx = (y_right - y_left) / (x_right - x_left);
    
    return y_left + dydx * (x - x_left);
  }

  /**
   * Calculate the N_ex for a given kappa by
   * interpolating over a vector of eps values
   * for each source.
   */
  real get_Nex(vector F, vector[] eps, vector kappa_grid, real kappa, real alpha_T) {

    real eps_from_kappa;
    int N = num_elements(F);
    real Nex = 0;

    for (i in 1:N-1) {
      eps_from_kappa = interpolate(kappa_grid, eps[i], kappa);
      Nex += F[i] * eps_from_kappa;
    }
    Nex += F[N] * (alpha_T / (4 * pi()));

    return Nex;
  }
  
  /**
   * Define the fik PDF.
   * NB: Cannot be vectorised.
   * Uses sinh(kappa) ~ exp(kappa)/2 
   * approximation for kappa > 100.
   */
  real fik_lpdf(vector v, vector mu, real kappa, real kappa_c) {

    real lprob;
    real inner = abs_val((kappa_c * v) + (kappa * mu));
    
    if (kappa > 100 || kappa_c > 100) {
      lprob = log(kappa * kappa_c) - log(4 * pi() * inner) + inner - (kappa + kappa_c) + log(2);
    }
    else {   
      lprob = log(kappa * kappa_c) - log(4 * pi() * sinh(kappa) * sinh(kappa_c)) + log(sinh(inner)) - log(inner);
    }

    return lprob;   
  }
  
}

data {

  /* sources */
  int<lower=0> Ns;
  unit_vector[3] varpi[Ns]; 
  vector[Ns] D;
  real Dbg;
  
  /* uhecr */
  int<lower=0> N; 
  unit_vector[3] arrival_direction[N]; 
  vector[N] zenith_angle;
  vector[N] A;
  
  /* observatory */
  real<lower=100, upper=10000> kappa_c;  
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
  //F0 = L0 / (4 * pi() * pow(Dbg, 2));
  F[Ns + 1] = F0;

  FT = (F0 + Fs);

  f = Fs / FT;

  if (is_nan(f)) {
    print("Fs: ", Fs);
    print("F0: ", F0);
  }
  
}

model {

  vector[Ns + 1] log_F;
  real Nex;
  
  log_F = log(F);

  /* Nex */
  Nex = get_Nex(F, eps, kappa_grid, kappa, alpha_T);

  /* rate factor */
  for (i in 1:N) {

    vector[Ns + 1] lps = log_F;

    for (k in 1:Ns + 1) {
      if (k < Ns + 1) {
	lps[k] += fik_lpdf(arrival_direction[i] | varpi[k], kappa, kappa_c);
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
