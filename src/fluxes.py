import numpy as np

from scipy.interpolate import interp1d
from scipy.optimize import minimize_scalar

def flux(a, b, u):
    """
    Godunov flux function
    Equation 2.5
    """
    if np.abs(a-b) < 1e-10:
        return f(a)
    elif a < b:
        u_cont = interp1d(u, f(u), fill_value="extrapolate")
        return minimize_scalar(u_cont, bounds=(a, b), method='bounded').x
    elif a > b:
        u_cont = interp1d(u, -f(u), fill_value="extrapolate")
        return minimize_scalar(u_cont, bounds=(b, a), method='bounded').x