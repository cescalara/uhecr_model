import numpy as np
from scipy import integrate, interpolate
from scipy.optimize import bisect

from fancy.propagation.energy_loss import get_arrival_energy


def P_x_z(z, x, Esrc, D_grid, E_grid, Earr_grid, B, alpha):
    """
    Function for calculation of P(x, z | E, alpha, B).
    E in EeV, B in nG, D in Mpc.

    For details, see Equation 13 in Capel & Mortlock (2019).

    @author Francesca Capel
    @date October 2018
    """

    # x, z ->  D, theta
    D = np.sqrt(x**2 + z**2)
    costheta = z / D
    sintheta = x / D

    # Constant part of \bar{\theta} dependence
    theta0 = np.deg2rad(2.3) * 50 / np.sqrt(10)
    kappa0 = 2.3 / theta0**2
    
    # Set up interpolation
    get_Esrc = interpolate.interp1d(D_grid, Esrc)
    get_Earr = interpolate.RectBivariateSpline(E_grid, D_grid, Earr_grid)

    # Initialise output
    P_xz = np.zeros(np.shape(D))

    # Handle single values (useful for testing)
    if isinstance(D, float) or isinstance(D, int):
        _d = D
        e_src = get_Esrc(_d)
        ct = costheta
        st = sintheta
        inner1 = kappa0 *  e_src**2 * (1 / _d) * ct * B**-2
        inner2 = kappa0 * e_src**2 * (1 / _d) * B**-2

        term1 = _d**-3 * e_src**(2 - alpha) * B**-2

        # Handle overflow in sinh
        if inner2 > 100:
            term2 = np.exp(inner1 - inner2 + np.log(2))
        else:
            term2 = np.exp(inner1) / np.sinh(inner2)

        term3 = 1

        P_xz = term1 * term2 * term3

    # Handle grid of values
    else:

        # Loop over elements
        for i, row in enumerate(D):
            
            for j, _d in enumerate(row):

                e_src = get_Esrc(_d)
                ct = costheta[i][j]
                st = sintheta[i][j]
                inner1 = kappa0 *  e_src**2 * (1 / _d) * ct * B**-2
                inner2 = kappa0 * e_src**2 * (1 / _d) * B**-2

                # Find dEarr/dE
                if (e_src < 5e3):
                    dEarrdE = get_Earr(e_src, _d, dx = 1, grid = False)
                else:
                    dEarrdE = 1
                    
                term1 = _d**-3 * e_src**(2 - alpha) * B**-2
                
                # Handle overflow in sinh
                if inner2 > 100:
                    term2 = np.exp(inner1 - inner2 + np.log(2))
                else:
                    term2 = np.exp(inner1) / np.sinh(inner2)

                term3 = 1 / abs(dEarrdE)

                P_xz[i][j] = term1 * term2 * term3
            
    return P_xz
