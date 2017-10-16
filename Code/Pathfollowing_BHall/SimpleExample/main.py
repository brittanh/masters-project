#!/opt/local/bin/python
# -*- encoding: ascii -*-
"""
    @purpose: Solving an NLP problem using a path following algorithm
    @author: Brittany Hall
    @date: 18.09.2017
    @version: 0.1
    @updates:
"""
from numpy import array, zeros, linspace, meshgrid, arange
import matplotlib.pyplot as plt
from problem import prob, obj
from nlp_solve import *
from pathfollowing import *

#Initial Values
p_init = array([1,-4])                                 #initial parameter value
p_final = array([8,1])                                   #final parameter value
x_init = array([[0.5],[0.6]])                              #initial primal variable
y_init = array([1.2])                                    #initial dual variable

"""
Plotting objective function and constraints (at initial and final state)
"""
#initial state
x1_list = linspace(0,1.6)
x2_list = linspace(0,1)
X1, X2 = meshgrid(x1_list, x2_list)
OBJ = p_init[0]*X1**3 + X2**2
plt.figure()
plt.contour(X1, X2, OBJ)
plt.plot(x1_list, exp(-x1_list),'b')
plt.plot(0.5,0.6,'r*')
plt.title('Contour Plot ($p_0$)')
plt.xlabel('$x_1$')
plt.ylabel('$x_2$')
plt.xlim(0,1.6)
plt.ylim(0,1)
plt.show()
#final state
OBJ = p_final[0]*X1**3 + X2**2
plt.figure()
plt.contour(X1, X2, OBJ)
plt.plot(x1_list, exp(-x1_list),'b')
plt.axvline(x=1, color='b')
plt.plot(1,exp(-1),'r*')
plt.title('Contour Plot ($p_f$)')
plt.xlabel('$x_1$')
plt.ylabel('$x_2$')
plt.xlim(0,1.6)
plt.ylim(0,1)
plt.show()

"""
Solving the problem
"""
#Solving NLP at p0 to get initial values
x_opt, lam_opt, mu_opt, con = nlp_solve(prob, obj, p_init, x_init, y_init)

#define method to use (predictor or predictor-corrector)
case = 'predictor-corrector'

#Solving the NLP to get optimal parameters using path-following algorithm
x_init, y_init, t_list, x_list_0, x_list_1, lam_list, mu_list, p = pathfollowing(p_init, p_final, x_init, x_opt, y_init, lam_opt, mu_opt, case)

#print(x_list_0)
#print(x_list_1)
#print(t_list)
#print(p)

plt.figure()
plt.plot(t_list, x_list_0,'.')
plt.title('$x_1$ as a function of $t$')
plt.xlabel('$t$')
plt.ylabel('$x_1$')
plt.xlim(0,1)
plt.ylim(0.2,1)
plt.xticks(arange(0,1.1,0.1))
plt.yticks(arange(0.2,1.1,0.1))
plt.show()
