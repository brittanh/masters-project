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
from numpy import zeros, ones, array, transpose, append, matlib, tile
import scipy.io as spio
from itPredHorizon import *

def optProblem(x, u, x0_measure, N, params):
    
    global nx, nu, nk, d, tf, ns
    NT = params['dist']['NT']
    Uf = params['dist']['F_0']
    
    #Modeling the system
    _, state, xdot, inputs = ColCSTR_model(Uf,params)
    f = Function('f', [state, inputs], [xdot])

    #Unpacking parameters
    x_min = params['bounds']['x_min']
    x_max = params['bounds']['x_max']
    u_min = params['bounds']['lbu']
    u_max = params['bounds']['ubu']

    #Loading steady state data
    data = spio.loadmat('CstrDistXinit.mat', squeeze_me = True)
    Xinit = data['Xinit']
    xf = Xinit[0:84]
    u_opt = Xinit[84:89]

    #Problem dimensions
    nx = 2*NT+2                  #Number of states (CSTR + Distillation Column)
    nu = 5                                  #Number of inputs (LT, VB, F, D, B)
    nk = 1
    tf = 1                                                               #[min]
    h = tf/nk
    ns = 0

    #Collecting model variables
    u = tile(u,nk)
    model = {'NT': NT, 'f': f, 'xdot_val_rf_ss': xf, 'x': x, 'u_opt': u_opt, 'u':u}
    params['model'] = model

    #Preparing collocation matrices
    _, C, D, d = collocationSetup() #function from casadi

    #Collecting collocation variables
    colloc = {'C': C, 'D': D, 'h': h}
    params['colloc'] = colloc

    delta_t = 1
    alpha = 1
    beta = 1
    gamma = 1

    #Weight variables
    weight = {'delta_t': delta_t, 'alpha': alpha, 'beta': beta, 'gamma': gamma}
    params['weight'] = weight

    #Initial conditions
    X0 = MX.sym('X0', nx, 1)
    w = {}                                 #Decision variables (control + state)
    w['X0'] =  X0
    lbw = x_min                             #Lower bound for decision variables
    ubw = x_max                             #Upper bound for decision variables
    w0 = x[0,0:nx]                                               #Initial guess
    J = 0                                                   #Objective function
    g = X0-x0_measure                                     #Nonlinear constraint
    lbg = params['bounds']['lbg']         #Lower bound for nonlinear constraint
    ubg = params['bounds']['ubg']         #Upper bound for nonlinear constraint

    Xk = X0
    
    data = spio.loadmat('Qmax.mat', squeeze_me = True)
    Qmax = data['Qmax']
    params['Qmax'] = Qmax
  
    count = 2                                       #Counter for state variable
    ssoftc = 0
    for iter in range(0,N):
      J, g, w0, w, lbg, ubg, lbw, ubw, Xk, params, count, ssoftc = itPredHorizon(Xk, w, w0, lbw, ubw, lbg, ubg, g, J, params, iter, count, ssoftc, nk, d, nu, nx)

    return J, g, w0, w, lbg, ubg, lbw,ubw, params
