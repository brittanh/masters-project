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

def optProblem(x, u, N, x0_measure):
    #NT = 41                                                   #Number of stages
    Uf = 0.3                                           #Feed rate to CSTR (F_0)

    #Calling the model
    _, state, xdot, inputs = ColCSTR_model(Uf)
    f = Function('f', [state, inputs], [xdot])

    #constraints
    VB_max = 4.008
    #State bounds and initial guess
    x_min = zeros((84,1))
    x_max = ones((84,1))
    #Control bounds
    u_min = array([[0.1],[0.1],[0.1],[0.1],[0.1]])
    u_max = array([[10],[VB_max],[10],[1.0],[1.0]])

    #Collecting all bounds
    bounds = {'x_min': x_min, 'x_max': x_max, 'u_min': u_min, 'u_max': u_max}
    params = {'bounds': bounds}

    #loading initial data
    data = spio.loadmat('CstrDistXinit.mat', squeeze_me = True)
    Xinit = data['Xinit']
    xf = Xinit[0:84]
    u_opt = Xinit[85:89]

    #Prices to be used
    pf = 1
    pV = 0.02
    pB = 2
    pD = 0

    #collecting prices
    price = {'pf': pf, 'pV': pV, 'pB': pB, 'pD': pD}
    params['price'] = price

    #Controller gains
    KcB = 10
    KcD = 10

    #Nominal hold ups
    MDs = 0.5
    MBs = 0.5

    #Nominal flows
    Ds = 0.5
    Bs = 0.5

    #Collecting gains
    gain = {'KcB': KcB, 'KcD': KcD, 'MDs': MDs, 'MBs': MBs, 'Ds': Ds, 'Bs': Bs}
    params['gain'] = gain

    global nx, nu, nk, d, tf, ns

    #problem dimensions
    nx = 84 #number of states (CSTR + Distillation Column)
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
              
