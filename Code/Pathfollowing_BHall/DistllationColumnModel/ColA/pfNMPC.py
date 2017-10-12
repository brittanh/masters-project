#!/opt/local/bin/python
# -*- encoding: ascii -*-
"""
    @purpose: Path- following based Nonlinear Model Predictive Control (pfNMPC)
    @author: Brittany Hall
    @date: 07.10.2017
    @version: 0.1
    @updates:
"""
from numpy import size, zeros, append
from solveOpt import *

import scipy.io as spio

def pfNMPC(optProblem, system, MPCit, N, T, tmeasure, xmeasure, u0, varargin):
    #Dimension of state and input
    nx = size(xmeasure) #Elements in state
    nu = size(u0, axis = 0) #Size of inputs
    #Constructing empty arrays for later use
    Tall = []
    Xall = zeros((MPCit, size(xmeasure, axis = 0)))
    Uall = zeros((MPCit, size(u0, axis = 0)))
    xmeasureAll = []
    uAll = []
    runtime = []
    u_pf_opt = []
    x_pf_opt = []

    #starting NMPC iteration
    iter = 1
    z1 = xmeasure
    #loading in noise data
    data = spio.loadmat('noise1pct.mat', squeeze_me = True)
    noise = data['noise']
    raw_input()
    while (iter <= MPCit):
        print('------------------------------------------------------------\n')
        print('MPC iteration: %d\n', iter)

        #Obtaining new initial value
        def measureInitVal(tmeasure, xmeasure):
            t0 = tmeasure
            x0 = xmeasure
            return t0, x0
        t0,x0 = measureInitVal(tmeasure, xmeasure)

        #adding measurement noise
        n_M = noise[:,iter]                                       #Holdup noise
        n_X = zeros((NT+1,1))                              #Concentration noise
        measure_noise = append(n_M, n)
        x0_measure = x0 + measure_noise     #Adding measurement noise to states

        #advanced-step NMPC
        primalNLP, dualNLP, lb, ub, _, params = solveOpt(optProblem, system, N, t0, x0, u0, T, iter, u_pf_opt, x_pf_opt, z1)

        #re-arrange NLP solutions
        _, x_nlp_opt = plotStates(primalNLP,lb,ub,N)

        p_init = primalNLP[0:nx]
        p_final = x0_measure
        x_start = primalNLP
        y_start = dualNLP

        delta_t = 0.5 #number of path following steps
        lb_init = lb
        ub_init = ub

        #NLP sensitivity (predictor-corrector)
        
