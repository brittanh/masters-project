#!/opt/local/bin/python
# -*- encoding: ascii -*-
"""
    @purpose: solving optimal control problem
    @author: Brittany Hall
    @date: 07.10.2017
    @version: 0.1
    @updates:
"""
from casadi import *
from ColCSTR_model import *
from numpy import zeros, ones, array, transpose, append
import scipy.io as spio
from itPredHorizon import *

def optProblem(x, u, N, x0_measure, params):
    #global nx, nu, nk, d, tf, ns
    
    NT = params['dist']['NT']
    Uf = 0.3                                           #Feed rate to CSTR (F_0)

    #Calling the model
    _, state, xdot, inputs = ColCSTR_model(Uf,params)
    raw_input()
    f = Function('f', [state, inputs], [xdot])

    #loading initial data
    data = spio.loadmat('CstrDistXinit.mat', squeeze_me = True)
    Xinit = data['Xinit']
    xf = Xinit[0:84]
    u_opt = Xinit[85:89]

    #problem dimensions
    nx = 2*NT+2 #number of states (CSTR + Distillation Column)
    nu = 5 #number of inputs (LT, VB, F, D, B)
    nk = 1
    tf = 1 #[min]
    h = tf/nk
    ns = 0

    #Collecting model variables
    model = {'NT': NT, 'f': f, 'xdot_val_rf_ss': xf, 'x': x, 'u_opt': u_opt, 'u': tile(u,(1,nk))}
    params['model'] = model

    #Preparing collocation matrices
    _, C, D, d = collocationSetup() #function from casadi

    #Collecting collocation variables
    colloc = {'C': C, 'D': D, 'h': h}
    params['colloc'] = colloc

    #Creating empty NLP
    w = []                                #Decision variables (control + state)
    w0 = []                                                      #Initial guess
    lbw = []                                #Lower bound for decision variables
    ubw = []                                #Upper bound for decision variables
    J = 0                                                   #Objective function
    g = []                                                #Nonlinear constraint
    lbg = []                              #Lower bound for nonlinear constraint
    ubg = []                              #Upper bound for nonlinear constraint

    delta_t = 1
    alpha = 1
    beta = 1
    gamma = 1

    #Weight variables
    weight = {'delta_t': delta_t, 'alpha': alpha, 'beta': beta, 'gamma': gamma}
    params['weight'] = weight

    #Initial conditions
    X0 = MX.sym('X0', nx)
    w = append(w[:], X0)
    lbw = append(lbw, x_min)
    ubw = append(ubw, x_max)
    w0 = append(w0,transpose(x[0,0:nx]))
    g = append(g,X0-x0_measure)
    lbg = append(lbg, zeros((nx,1)))
    ubg = apppend(ubg, zeros((nx,1)))

    Xk = X0
    data = spio.loadmat('Qmax.mat', squeeze_me = True)
    Qmax = data['Qmax']
    params['Qmax'] = Qmax
  
    count = 2 #counter for state variable
    ssoftc = 0
    for i in range(0,N):
      J, g, w0, w, lbg, ubg, lbw, ubw, Xk, params, count, ssoftc = itPredHorizon(Xk, w, w0, lbw, ubw, lbg, ubg, g, J, params, i, count, ssoftc)
              
    return J, g, w0, w, lbg, ubg, lbw,ubw, params
