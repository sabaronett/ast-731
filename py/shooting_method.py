#!/usr/bin/env python3
#==============================================================================
# shooting_method.py
#
# Use the shooting method with the fourth-order Runge--Kutta integration scheme
# to solve the Lane--Emden equation (Hansen, Kawaler & Trimble, 2004, eq. 7.44)
# for polytropic parameters x = ξ, y = θₙ, and z = dy/dz = dθₙ/dξ.
#
# Author: Stanley A. Baronett
# Created: 2022-11-19
# Updated: 2022-12-09
#==============================================================================

# From Hansen, Kawaler & Trimble (2004, eqs. 7.45 and 7.48)
yfunc = lambda x, n : 1 - (1/6)*x**2 + (n/120)*x**4 - (n*(8*n - 5)/15120)*x**6
zfunc = lambda x, n : -(1/3)*x + (n/30)*x**3 - (n*(8*n - 5)/2520)*x**5
yprime = lambda x, y, z, n : z
zprime = lambda x, y, z, n : -y**n - (2/x)*z

def rk4(n, xi, yi, zi, h):
    """
    ... from Hansen, Kawaler & Trimble (2004, eq. 7.46)
    """
    k1 = h*yprime(xi, yi, zi, n)
    l1 = h*zprime(xi, yi, zi, n)
    k2 = h*yprime(xi + h/2, yi + k1/2, zi + l1/2, n)
    l2 = h*zprime(xi + h/2, yi + k1/2, zi + l1/2, n)
    k3 = h*yprime(xi + h/2, yi + k2/2, zi + l2/2, n)
    l3 = h*zprime(xi + h/2, yi + k2/2, zi + l2/2, n)
    k4 = h*yprime(xi + h, yi + k3, zi + l3, n)
    l4 = h*zprime(xi + h, yi + k3, zi + l3, n)
    xi1 = xi + h
    yi1 = yi + k1/6 + k2/3 + k3/3 + k4/6
    zi1 = zi + l1/6 + l2/3 + l3/3 + l4/6

    return xi1, yi1, zi1
