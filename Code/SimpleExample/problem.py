#!/opt/local/bin/python
# -*- encoding: ascii -*-
"""
    @purpose: Defining the problem to be solved
    @author: Brittany Hall
    @date: 18.09.2017
    @version: 1
    @updates:
"""

from casadi import *
from numpy import array


#Defining the problem
def prob():
    """
        Information on the problem to be solved
        """
    n = 2                                                  #number of variables
    neq = 0                                     #number of equality constraints
    niq = 2                                   #number of inequality constraints

def obj(x,y,p):
    """
        Problem to be solved and calculation of derivates, Lagrangians,etc
        """
    x = SX.sym('x',2)                             #x variable (primal solution)
    p = SX.sym('p',2)                               #p variable (dual solution)
    f = p[0]*x[0]**3+x[1]                                   #objective function
    con = vertcat([exp(-x[0])-x[1],p[1]-x[0]])
    
    #Constructing the Lagrangian
    lag = f+ mul((y.T),con)
    
    #Calculating gradient of F
    opf = {'input':['x','p'],\
        'output': ['f']}
    f = SXFunction('f',[x,p],[f],opf) #making the objective function a function
    g = SXFunction(f.gradient())            #gradient of the objective function

#Calculating derivative of Lagrangian and constraints
op = {'input':['x','p'],\
    'output':['lag']}
    lagr = SXFunction('lagr',[x,p],[lag],op)
    H = SXFunction(lagr.hessian('x','lag'))                            #Hessian
    J_lagr_x = SXFunction(lagr.jacobian('x','lag'))     #Jacobian of Lagrangian
    Lxp = J_lagr_x.jacobian(1,0)
    
    opc = {'input':['x','p'],\
        'output':['con']}
    cst = SXFunction('cst',[x,p],[con],opc)
    J_cst = cst.jacobian('x','con')
    cp = cst.jacobian('p','con')

#Defining Input values and evaluating
f.SetInput(x,'x')
    f.setInput(p,'p')
    g.setInput(x,'x')
    g.setInput(p,'p')
    H.setInput(x,'x')
    H.setInput(p,'p')
    Lxp.setInput(x,'x')
    Lxp.setInput(p,'p')
    J_cst.setInput(x,'x')
    J_cst.setInput(p,'p')
    cp.setInput(x,'x')
    cp.setInput(p,'p')
    cst.setInput(x,'x')
    cst.setInput(p,'x')
    
    f.evaluate()
    g.evaluate()
    H.evaluate()
    Lxp.evaluate()
    J_cst.evaluate()
    cp.evaluate()
    Lxp.evaluate()
    cst.evaluate()
    
    f = array(f.getOutput())
    g   = array(g.getOutput())
    H   = array(H.getOutput())
    Lxp = array(Lxp.getOutput())
    J_cst   = array(J_cst.getOutput())
    cp  = array(cp.getOutput())
    cst = array(cst.getOutput())
    
    #Equality Constraints
    Jeq = array([])
    dpe = array([])
    
    return f,g,H,Lxp,cst,J_cst,cp,Jeq,dpe


