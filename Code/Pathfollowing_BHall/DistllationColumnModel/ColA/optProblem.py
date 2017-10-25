#!/opt/local/bin/python
# -*- encoding: ascii -*-
"""
    @purpose: Solving the optimal control problem
    @author: Brittany Hall
    @date: 07.10.2017
    @version: 0.1
    @updates:
"""
from casadi import *
from collocationSetup import *
from ColCSTR_model import *
from numpy import zeros, ones, array, transpose, matlib, tile
import scipy.io as spio
from itPredHorizon import *

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
    nx = params['prob']['nx']    #Number of states (CSTR + Distillation Column)
    nu = params['prob']['nu']               #Number of inputs (LT, VB, F, D, B)
    nk = params['prob']['nk']
    tf = params['prob']['tf']                                            #[min]
    h =  params['prob']['h']
    ns = params['prob']['ns']

    #Collecting model variables
    u = tile(u,nk)
    model = {'NT': NT, 'f': f, 'xdot_val_rf_ss': xf, 'x': x, 'u_opt': u_opt, 'u':u}
    params['model'] = model

    #Preparing collocation matrices
    _, C, D, d = collocationSetup() #function from casadi

    #Collecting collocation variables
    colloc = {'C': C, 'D': D, 'h': h}
    params['colloc'] = colloc
    
    #Empty NLP
    w = MX()                              #Decision variables (control + state)
    w0 = array([])                                               #Initial guess
    lbw = array([])                          #Lower bound for decision variable
    ubw = array([])                          #Upper bound for decision variable
    g = MX()                                              #Nonlinear constraint
    lbg = array([])                       #Lower bound for nonlinear constraint
    ubg = array([])                       #Upper bound for nonlinear constraint
    J = 0                                                   #Objective function
    
    delta_t = 1
    alpha = 1
    beta = 1
    gamma = 1

    #Weight variables
    weight = {'delta_t': delta_t, 'alpha': alpha, 'beta': beta, 'gamma': gamma}
    params['weight'] = weight
   
    #Initial conditions
    X0 = MX.sym('X0', nx)
    w =  vertcat(w,X0)
    lbw = append(lbw,x_min)
    ubw = append(ubw,x_max)
    w0 = append(w0, x[0,0:nx])
    g = vertcat(g, X0-x0_measure)
    lbg = append(lbg, params['bounds']['lbg'])
    ubg = append(ubg, params['bounds']['ubg'])

    #Formulating the NLP
    Xk = X0
    data = spio.loadmat('Qmax.mat', squeeze_me = True)
    Qmax = data['Qmax']
    params['Qmax'] = Qmax
  
    count = 2                                       #Counter for state variable
    ssoftc = 0
    for iter in range(0,N):
      J, g, w0, w, lbg, ubg, lbw, ubw, Xk, params, count, ssoftc = itPredHorizon(Xk, w, w0, lbw, ubw, lbg, ubg, g, J, params, iter, count, ssoftc, d)

    return J, g, w0, w, lbg, ubg, lbw,ubw, params
