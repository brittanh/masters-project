#!/opt/local/bin/python
# -*- encoding: ascii -*-
"""
    @purpose: Setting up collocation
    @author: Brittany Hall (based on Joel Anderson's Matlab script)
    @date: 18.10.2017
    @version: 0.1
    @updates:
"""
from casadi import *
from numpy import zeros, convolve, polyval, polyder, polyint, array, append
def collocationSetup():
    #Degree of interpolating polynomial
    d = 3
    #Get collocation points
    tau_root = collocation_points(d, 'legendre')
    tau_root = append(0,tau_root)
    #Coefficients of the collocation equation
    C = zeros((d+1, d+1))
    #Coefficients of the continuity equation
    D = zeros((d+1, 1))
    #Coefficients of the quadrature function
    B = zeros((d+1, 1))

    #Construct polynomial basis
    for j in range(0,d+1):
        #Construct Lagrange polynomials to get the polynomial basis at the collocation point
        coeff = 1
        for r in range(0,d+1):
            if r != j:
                coeff = convolve(coeff, [1, -tau_root[r]])
                coeff = coeff/(tau_root[j]-tau_root[r])
        #Evaluate the polynomial at the final time to get the coefficients of the continuity equation
        D[j] = polyval(coeff,1.0)
        #Evaluate the time derivative of the polynomial at all collocation points to get the coefficients of the continuity equation
        pder = polyder(coeff)
        for r in range(0,d+1):
            C[j,r] = polyval(pder,tau_root[r])
        #Evaluate the integral of the polynomial to get the coefficients of the quadrature function
        pint = polyint(coeff)
        B[j] = polyval(pint, 1.0)
    return B,C,D,d
