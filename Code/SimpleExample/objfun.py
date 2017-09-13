#!/opt/local/bin/python
# -*- encoding: ascii -*-
"""
@purpose: Simple problem to test path-following algorithm
          Taken from Kungurtsev et al 
          ("SQP methods for parametric nonlinear optimization")
@author: Brittany Hall
@date: 12.09.2017
@version: 1
@updates:
"""

from casadi import *

#problem to solve
# min p1x1^3+x2^2
# st x2-exp(-x1)=>0, x1=> p2
#want to find approximate solution at p=(8,1)

#Approximate solution
x0 = SX([0.5,0.6]) #numerical matrix
y0 = SX([1.2])

p = SX([1,-4])

#Defining Constants and Variables
t = 0
x = MX.sym('x',2,1)
y = MX.sym('y',1)
p = MX.sym('p',2,1)

#Defining objective function
J = objfun('J',[x,p], [mtimes(p[0],x[0]**3)+x[1]**2])

#Defining constraints
cont = constraint('cont',[x,p],[x[1]-exp(-x[0]),x[0]-p[1]])

#Calculating Derivative/gradient
def df(x,y,p):
    return gradient(objfun(x,y,p))

#Calculating Jacobian
def Jacobian(x,y,p):
    return jacobian()
#Calculating optimality residual calculation


#Parametric SQP Method (path-following)
def algo_pathfollowing:
    while t<1:
    #Calculate Active set

