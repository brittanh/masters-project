#!/opt/local/bin/python
# -*- encoding: ascii -*-
"""
    @purpose: Defining the problem to be solved
    @author: Brittany Hall
    @date: 18.09.2017
    @version: 0.1
    @updates:
"""

from casadi import SX, Function, vertcat
from numpy import array, ones, zeros, exp


#Defining the problem
def prob():
    """
    Information on the problem to be solved
    """
    n = 2                                         #number of variables [x1, x2]
    np = 2                                       #number of parameters [p1, p2]
    neq = 0                                     #number of equality constraints
    niq = 2                                   #number of inequality constraints
    name = "Problem 1"
    return n, np, neq, niq, name

def obj(x, y, p, neq, niq, n, np):
    """
    Problem to be solved
    """
    p = SX.sym('p',np)                                              #Parameters
    x = SX.sym('x',n)                                                 #Variable
    f = p[0]*x[0]**3 + x[1]**2                                 #Objective array
    f_fun = Function('f_fun',[x,p],[p[0]*x[0]**3+x[1]**2])  #Objective function
    
    con = vertcat(exp(-x[0])-x[1],p[1]-x[0])                  #Constraint array
    conf = Function('conf',[x,p],[exp(-x[0])-x[1],p[1]-x[0]])#Constraint function
    
    
    #Specifying Bounds
    ubx = 1e16*ones([1,n])                                #Variable upper bound
    lbx = -1e16*ones([1,n])                               #Variable lower bound
    ubg = zeros([1,niq+neq])                            #Constraint upper bound
    lbg= -1e16*ones([1,niq+neq])                        #Constraint lower bound
    return x, p, f, f_fun, con, conf, ubx, lbx, ubg, lbg


