#!/opt/local/bin/python
# -*- encoding: ascii -*-
"""
    @purpose: Solving an NLP problem using a path following algorithm
    @author: Brittany Hall
    @date: 18.09.2017
    @version: 0.1
    @updates:
"""
from casadi import *
from numpy import array, zeros
import matplotlib.pyplot as plt
from problem import prob, obj
from nlp_solve import *
from pathfollowing import *

#Parameters
delta_t = 0.1                                                      #t step size
N = 1/delta_t                                             #number of iterations
alpha1 = 0.5

#Initial Values
t = 0
p_init = array([1,-4])                                 #initial parameter value
p_final = array([8,1])                                   #final parameter value
x_init = array([0.5,0.6])                              #initial primal variable
y_init = array([1.2])                                    #initial dual variable


#Solving NLP at p0 to get initial values

x_opt, lam_opt, mu_opt, con = nlp_solve(prob, obj, p_init, x_init, y_init)

chi = x_opt#initial chi
#print(chi)
#print(lam_opt)
#print(mu_opt)
#define method to use (predictor or predictor-corrector)
case = 'corrector'

#Solving the NLP to get optimal parameters using path-following algorithm
x_init, y_init, t_list, x_list_0, x_list_1 = pathfollowing(t, delta_t, alpha1, p_init, p_final, x_init, y_init, chi, lam_opt, mu_opt, case)

print(x_list_0)
print(x_list_1)
print(t_list)
