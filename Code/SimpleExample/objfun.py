#!/opt/local/bin/python
# -*- encoding: ascii -*-
"""
@purpose: Simple problem to test optimization using path-following algorithm
@author: Brittany Hall
@date: 12.09.2017
@version: 1
@updates:
"""

from casadi import *
from numpy import array, empty
#problem to solve
# min p1x1^3+x2^2
# st x2-exp(-x1)=>0, x1=> p2
#want to find approximate solution at p=(8,1)

#Defining Constants
t = 0
delta_t = 1e-4
N = 1/delta_t
alpha1 = 0.66 #0<alpha1<1

#Approximate solution
x0 = array([0.5,0.6]) #numerical matrix
p = empty([1, N])
p[0] = array([-1, 4]) #initial parameter value
p[N] = array([8,1]) #final parameter value

#initial variables from NLP
chi = MX.sym('X',2,1) #primal variable
lamb = MX.sym('lam',1) #dual variable
mu = MX.sym(',mu',2,1) #dual variable

#Calculate Lagrangian

#Calculate Jacobian of Lagrangian

#Calculate gradient of F

#Path Following Algorithm
##Algorithm 2 from 2016_Suwartadi_etal

for k in range(1,N):
    delta_p = p[k]-p[k-1] #calculating change in parameter p
    #solve QP

    #if QP is feasible:
    if :
        chi = chi + delta_chi
        lamb = delta_lamb
        mu = delta_mu
        t = t + delta_t
        k = k +1
    else:
        delta_t = alpha1*delta_t
        t = t- delta_t
