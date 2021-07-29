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
