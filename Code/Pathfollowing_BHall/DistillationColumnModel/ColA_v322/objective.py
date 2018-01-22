#!/opt/local/bin/python
# -*- encoding: ascii -*-
"""
    @purpose: Objective function to be solved
    @author: Brittany Hall
    @date: 08.10.2017
    @version: 0.1
    @updates:
"""
from casadi import MX, SX, vertcat, Function, mtimes, DM, jacobian, hessian, transpose
from numpy import size, append,multiply, shape
import scipy.io as spio
from itPredHorizon_pf import itPredHorizon_pf
from ColCSTR_model import ColCSTR_model
from collocationSetup import collocationSetup

def objective(x,y,p,N,params):
    
    nPrimal = x.numel()               #number of primal sln
    nDual = y['lam_g'].numel()               #number of dual
    nParam = p.numel()                  #number of parameters
    
    #Loading initial states and controls
    data = spio.loadmat('CstrDistXinit.mat', squeeze_me = True)
    Xinit = data['Xinit']
    xf = Xinit[0:84]
    u_opt = Xinit[84:]

    #Model parameters
    NT = params['dist']['NT']             #Stages in column
    Uf = 0.3                          #Feed rate to CSTR F_0
    _,state,xdot,inputs = ColCSTR_model(Uf,params)
    sf = Function('sf',[state,inputs],[xdot])

    params['model']['sf'] = sf
    params['model']['xdot_val_rf_ss'] = xf
    params['model']['x'] = x
    params['model']['u_opt'] = u_opt
    
    nx = params['prob']['nx']
    nu = params['prob']['nu']
    nk = params['prob']['nk']
    tf = params['prob']['tf']
    ns = params['prob']['ns']
    h = params['prob']['h']

    #Preparing collocation matrices
    _, C, D, d = collocationSetup()
    params['prob']['d'] = d
    colloc = {'C':C,'D':D, 'h':h}
    params['colloc'] = colloc
    
    #NLP variable vector
    V = MX()                              #Decision variables (control + state)
    obj = 0                                                 #Objective function
    cons = MX()                                          #Nonlinear Constraints
    
    delta_time = 1
    alpha = 1
    beta = 1
    gamma = 1

    params['weight']['delta_time'] = delta_time
    params['weight']['alpha'] = alpha
    params['weight']['beta'] = beta
    params['weight']['gamma'] = gamma
    
    #Initial states and Controls
    data_init = spio.loadmat('CstrDistXinit.mat', squeeze_me = True)
    xf = data_init['Xinit'][0:84]
    u_opt = data_init['Xinit'][84:89]
    
    #"Lift" Initial conditions
    X0 = MX.sym('X0', nx)
    V = vertcat(V,X0)                     #Decision variables
    cons = vertcat(cons, X0 - x[0:nx,0]) #Nonlinear constraints
    cons_x0 = X0 - x[0:nx,0]

    #Formulating the NLP
    Xk = X0
    
    data = spio.loadmat('Qmax.mat', squeeze_me = True)
    params['Qmax'] = data['Qmax']
    ssoftc = 0
    
    for i in range(0,N):
       obj, cons, V, Xk, params, ssoftc = itPredHorizon_pf(Xk, V, cons, obj, params, i, ssoftc)

    V = vertcat(V[:])
    cons = vertcat(cons[:])

    #Objective function and constraint functions
    f = Function('f', [V], [obj], ['V'], ['objective'])
    c = Function('c', [V], [cons], ['V'], ['constraint'])
    cx0 = Function('cx0',[X0], [cons_x0], ['X0'], ['constraint'])

    #Constructing Lagrangian
    lag_expr = obj + mtimes(transpose(y['lam_g']),cons)
    g = Function('g',[V],[jacobian(obj,V),obj])
    lagr = Function('lagr', [V], [lag_expr], ['V'], ['lag_expr'])
    [H,gg] = hessian(lag_expr,V)
    H = Function('H',[V],[H,gg])
    [Hobj,gobj] = hessian(obj,V)
    Hobj = Function('Hobj', [V], [Hobj,gobj])
    J = Function('J',[V],[jacobian(cons,V),cons])
    Jp = Function('Jp',[X0],[jacobian(cons_x0,X0),cons_x0])

    #Evaluating functions at current point
    f = f(x)
    g = g(x)
    g = g[0]
    g = transpose(g)
    H = H(x)
    H = H[0]
    Lxp = H[0:nPrimal,0:nParam]
    J = J(x)
    J = J[0]
    Jtemp = DM.zeros((nDual,nParam))
    cp = Jp(x[0:nParam])
    cp = cp[0]
    Jtemp[0:nParam, 0:nParam] = cp.full()
    cp = Jtemp.sparse()
    cst = c(x)

    #Evaluation of objective function used for Greshgorin bound
    Hobj = Hobj(x)
    Hobj = Hobj[0].sparse()
    f = f.full()
    g = g.sparse()
    Lxp = Lxp.sparse()
    cst = cst.full()
   
    #Equality constraint
    Jeq = J
    dpe = cp

    return f, g, H, Lxp, cst, J, cp, Jeq, dpe, Hobj
                               
                               
                   

