#!/opt/local/bin/python
# -*- encoding: ascii -*-
"""
    @purpose: Solving the optimal control problem
    @author: Brittany Hall
    @date: 07.10.2017
    @version: 0.1
    @updates:
"""
from casadi import Function, MX, SX, vertcat
from collocationSetup import collocationSetup
from ColCSTR_model import ColCSTR_model
from numpy import zeros, ones, array, transpose, matlib, tile, reshape, shape, savetxt
import scipy.io as spio
from itPredHorizon import itPredHorizon

def optProblem(x, u, x0_measure, N, params):
    
    NT = params['dist']['NT']
    Uf = params['dist']['F_0']
    
    #Modeling the system
    _, state, xdot, inputs = ColCSTR_model(Uf,params)
    f = Function('f', [state, inputs], [xdot])

    #Unpacking parameters
    x_min = params['bounds']['x_min']
    x_max = params['bounds']['x_max']

    #Loading steady state data
    data = spio.loadmat('CstrDistXinit.mat', squeeze_me = True)
    Xinit = data['Xinit']
    xf = Xinit[0:84]
    u_opt = Xinit[84:89]

    #Problem dimensions
    nx = params['prob']['nx']          #Number of states
    nu = params['prob']['nu']          #Number of inputs
    nk = params['prob']['nk']
    tf = params['prob']['tf']
    h =  params['prob']['h']
    ns = params['prob']['ns']

    #Collecting model variables
    u = tile(u,nk)
    model = {'NT': NT, 'f': f, 'xdot_val_rf_ss': xf,
            'x': x, 'u_opt': u_opt, 'u':u}
    params['model'] = model

    #Preparing collocation matrices
    _, C, D, d = collocationSetup() 
    params['prob']['d'] = d

    #Collecting collocation variables
    colloc = {'C': C, 'D': D, 'h': h}
    params['colloc'] = colloc
    
    #Empty NLP
    w = MX()               #Decision variables (control + state)
    w0 = []                                       #Initial guess
    lbw = []                  #Lower bound for decision variable
    ubw = []                  #Upper bound for decision variable
    g = MX()                               #Nonlinear constraint
    lbg = []               #Lower bound for nonlinear constraint
    ubg = []               #Upper bound for nonlinear constraint
    J = 0                         #Initialize objective function
    
    #Weight variables
    delta_t = 1
    alpha = 1
    beta = 1
    gamma = 1
    weight = {'delta_t': delta_t, 'alpha': alpha,
        'beta': beta, 'gamma': gamma}
    params['weight'] = weight
   
    #Initial conditions
    X0 = MX.sym('X0', nx)
    w =  vertcat(w,X0)
    w0 = [i for i in x[0,0:nx]]
    lbw = [i for i in x_min]
    ubw = [i for i in x_max]
    g = vertcat(g, X0-x0_measure)
    lbg = params['bounds']['lbg']
    ubg = params['bounds']['ubg']

    Xk = X0
    data = spio.loadmat('Qmax.mat', squeeze_me = True)
    Qmax = data['Qmax']
    params['Qmax'] = Qmax
  
    count = 2                         #Counter for state variable
    ssoftc = 0
    for iter in range(0,N):
        J, g, w0, w, lbg, ubg, lbw, ubw, Xk, params, count, ssoftc = itPredHorizon(Xk, w, w0, lbw, ubw, lbg, ubg, g, J, params, iter, count, ssoftc, d)

    return J, g, w0, w, lbg, ubg, lbw, ubw, params
