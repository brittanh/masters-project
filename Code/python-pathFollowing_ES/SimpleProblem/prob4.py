'''
Created on 18. nov. 2015

@author: suwartad
'''
from numpy import array
from casadi import *
#from scipy.sparse import csc_matrix

def prob4():
    '''
    information on the problem
    '''
    n   = 2
    neq = 0
    niq = 2
    name = "Problem 4"
    return n, neq, niq, name

def objective(x,y,p):
    '''
    definition of a problem and its derivatives
    '''
    # construct Lagrangian
    v = SX.sym('v',2)  # x variable resemble
    t = SX.sym('t',2) # p variable resemble
    f = t[0]*v[0]**3 + v[1]**2
    c_expr   = vertcat(exp(-v[0])-v[1],t[1]-v[0])
    lag_expr = f + mul(y.T,c_expr)
    
    # objective function related derivatives
    opf = {'input_scheme': ['v','t'], \
           'output_scheme': ['f']}
    f = SXFunction('f',[v,t],[f],opf)
    g = SXFunction(f.gradient());
    
    # Lagrangian and constraint related derivatives
    op = {'input_scheme': ['v','t'], \
          'output_scheme': ['lag_expr']}
    lagr       = SXFunction('lagr',[v,t],[lag_expr],op)
    H          = SXFunction(lagr.hessian('v','lag_expr'))
    jac_lagr_x = SXFunction(lagr.jacobian('v','lag_expr'))
    Lxp        = jac_lagr_x.jacobian(1,0)
    
    # inequality constraint
    opc = {'input_scheme': ['v','t'], \
           'output_scheme': ['c_expr']}
    cst = SXFunction('cst',[v,t],[c_expr],opc)
    J   = cst.jacobian('v','c_expr')
    cp  = cst.jacobian('t','c_expr')
    
    # set input values, evaluate, and obtain their values
    f.setInput(x,'v')
    f.setInput(p,'t')
    g.setInput(x,'v')
    g.setInput(p,'t')
    H.setInput(x,'v')
    H.setInput(p,'t')
    Lxp.setInput(x,'v')
    Lxp.setInput(p,'t')
    J.setInput(x,'v')
    J.setInput(p,'t')
    cp.setInput(x,'v')
    cp.setInput(p,'t')
    cst.setInput(x,'v')
    cst.setInput(p,'t') 

    f.evaluate()
    g.evaluate()
    H.evaluate()
    J.evaluate()
    cp.evaluate()
    Lxp.evaluate()
    cst.evaluate()
    
    f   = array(f.getOutput())
    '''
    g   = csc_matrix(g.getOutput())
    H   = csc_matrix(H.getOutput())
    Lxp = csc_matrix(Lxp.getOutput())
    J   = csc_matrix(J.getOutput())
    cp  = csc_matrix(cp.getOutput())
    '''
    g   = array(g.getOutput())
    H   = array(H.getOutput())
    Lxp = array(Lxp.getOutput())
    J   = array(J.getOutput())
    cp  = array(cp.getOutput())
    cst = array(cst.getOutput())
    
    # equality constraint
    Jeq = array([])
    dpe = array([])
    
    return f,g,H,Lxp,cst,J,cp,Jeq,dpe
